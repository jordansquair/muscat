x <- x1
ng <- 500
nc <- 2e3
nk <- ns <- 3
kp <- sp <- gp <- NULL
pdd <- c(0.2, 0.2, rep(0.15, 4))
lfc <- list(dist = "gamma", shape = 4, rate = 2)
paired <- FALSE

# simData <- function(x, ng, nc, nk, ns, 
#     kp = NULL, sp = NULL, gp = NULL,
#     pdd = diag(6)[1, ], lfc = list("gamma", shape = 4, rate = 2)) {

    # throughout this code...
    # k = cluster ID
    # s = sample ID
    # g = group ID
    # c = DD category
    # 0 = reference
    
    # get reference IDs
    names(kids0) <- kids0 <- levels(x$cluster_id)
    names(sids0) <- sids0 <- levels(x$sample_id)
    nk0 <- ifelse((nk0 <- length(kids0)) == 0, 1, nk0)
    ns0 <- ifelse((ns0 <- length(sids0)) == 0, 1, ns0)
    
    # set simulation IDs & gene, cell names
    names(kids) <- kids <- paste0("cluster", seq_len(nk))
    names(sids) <- sids <- paste0("sample", seq_len(ns))
    names(gids) <- gids <- c("A", "B")
    gs <- paste0("gene", seq_len(ng))
    cs <- paste0("cell", seq_len(nc))
    
    # default to equal probability for each cluster, sample, group ID
    if (is.null(kp)) kp <- rep(1/nk, nk)
    if (is.null(sp)) sp <- rep(1/ns, ns)
    if (is.null(gp)) gp <- rep(0.5, 2)
    
    # sample cell metadata
    cd <- data.frame(
        cluster_id = sample(kids, nc, TRUE, kp),
        sample_id = sample(sids, nc, TRUE, sp),
        group_id = sample(gids, nc, TRUE, gp),
        row.names = cs, stringsAsFactors = FALSE)
    
    # get cell indices for each cluster-sample-group
    dt <- data.table(cell_id = rownames(cd), cd)
    ci <- split(dt, by = colnames(cd), sorted = TRUE, flatten = FALSE)
    ci <- map_depth(ci, ncol(cd), "cell_id")
    
    # split reference cells by cluster-sample
    cd0 <- as.data.frame(colData(x))
    dt0 <- data.table(cell_id = colnames(x), cd0)
    dt0 <- split(dt0, by = c("cluster_id", "sample_id"), sorted = TRUE, flatten = FALSE)
    cs0 <- map_depth(dt0, 2, "cell_id")
    
    # sample number of genes to simulate per DD category
    ndd <- vapply(kids, function(k) 
        table(sample(cats, ng, TRUE, pdd)), 
        numeric(length(cats)))
    
    # get gene indices for each cluster-category
    gi <- vapply(kids, function(k) 
        split(sample(gs), rep.int(cats, ndd[, k])), 
        vector("list", length(cats)))
    
    # sample reference genes to simulate from 
    gs0 <- set_names(sample(rownames(x), ng, ng > nrow(x)), gs)
    
    # sample reference clusters & samples
    kids0_use <- setNames(sample(kids0, nk, nk > nk0), kids)
    if (paired) { 
        # use same set of reference samples for all groups
        sids0_use <- sample(sids0, ns, ns > ns0)
        sids0_use <- replicate(length(gids), sids0_use)
    } else {
        # draw reference samples at random for each group
        sids0_use <- replicate(length(gids), sample(sids0, ns, ns > ns0))
    }
    dimnames(sids0_use) <- list(sids, gids)
    
    # sample logFCs for each cluster-category
    fun <- get(paste0("r", lfc$dist))
    lfc <- lfc[-grep("dist", names(lfc))]
    lfc <- vapply(kids, function(k)
        lapply(unfactor(cats), function(c) {
            n <- ndd[c, k]
            if (c == "ee")
                return(rep(NA, n))
            # sample directions for each gene with 
            # equal probability of up- & down-regulation
            dirs <- sample(c(-1, 1), n, TRUE) 
            # sample logFCs from input distribution
            lfcs <- set_names(do.call(fun, c(n, lfc)), gs0[[k]][[c]])
        }), vector("list", length(cats)))
    
    # fit offsets for each cluster-sample
    fits <- map_depth(dt0, 2, function(u) {
        fits <- fitdistr(u$offset, "normal")
        as.list(fits$estimate)
    })
    
    mapply(function(k, s) {
        dt <- data.table(cell_id = rownames(cd), cd)
        ci <- split(dt, by = colnames(cd), sorted = TRUE, flatten = FALSE)
        ci <- map_depth(ci, ncol(cd), "cell_id")
        dt <- data.table(cell_id = colnames(x), as.data.frame(colData(x)))
        cs0 <- split(dt, by = c("cluster_id", "sample_id"), sorted = TRUE, flatten = FALSE)
        cs0 <- map_depth(cs0, 2, "cell_id")
        
        # get reference cluster, samples, cells
        k0 <- kids0_use[[k]]
        s0 <- sids0_use[s, ]
        cs0 <- cs0[[k0]][s0]
        
        # get simulation cell indices
        ci <- ci[[k]][[s]]
        
        # get number of cells to simulate per group
        ncs <- vapply(ci, length, numeric(1))
        
        lapply(cats[ndd[, k] != 0], function(c) {
            # get reference genes
            gs0 <- gs0[[k]][[c]]
            
            # get simulation gene indices
            gi <- gi[[c, k]]
            
            # compute NB parameters
            ds <- rowData(x)[gs0, "dispersion"] # get dispersions
            os <- lapply(gids, function(g) {    # sample offsets
                fit <- fits[[k0]][[s0[g]]]
                do.call(rnorm, c(ncs[[g]], fit))
            })
            ms <- lapply(gids, function(g) { # compute NB means
                ms <- ds %o% os[[g]]
                dimnames(ms) <- list(gi, ci[[g]])
                return(ms)
            })      
            
            # sample counts
            
        })
    })
}

#' @param nc,ng,ns,nk number of cells, genes, samples, clusters to simulate. 
#' @param sp,kp,gp probability of cells to be in each sample, cluster, group.
# simData2 <- function(x, ng, nc, nk, ns, kp = NULL, sp = NULL, gp = NULL,
#     pdd = diag(6)[1, ], lfc = list(dist = "gamma", shape = 4, rate = 2),
#     paired = FALSE) {
  

simPars <- function(x,
    kid = "cluster_id", sid = "sample_id", gid = NULL, guse = NULL,
    min_count = 1, min_cells = 10, min_genes = 100, min_size = 100, 
    verbose = TRUE) {
    
    # store list of input arguments
    args <- as.list(environment())

    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"),
        is.null(kid) | is.character(kid) & length(kid) == 1 & !is.null(x[[kid]]),
        is.null(sid) | is.character(sid) & length(sid) == 1 & !is.null(x[[sid]]),
        is.null(gid) | is.character(gid) & length(gid) == 1 & !is.null(x[[gid]]),
        is.logical(verbose), length(verbose) == 1)

    data(sce); x <- sce
    args <- list(kid = "cluster_id", sid = "sample_id", gid = "group_id")
    for (id in names(args)) assign(id, args[[id]])
    min_count <- 1; min_cells = 10; min_genes = 100; min_size = 100; guse <- NULL
    
    # convert cell metadata ID columns to factors
    ids <- grep("id", names(args), value = TRUE)
    for (id in args[ids]) x[[id]] <- factor(x[[id]])
    
    nk <- length(kids <- levels(x[[kid]]))
    ns <- length(sids <- levels(x[[sid]]))
    
    # default to using all cells
    if (!is.null(gid)) {
        # default to using first group available
        if (is.null(guse)) {
            guse <- levels(x[[gid]])[1] 
        } else {
            stopifnot(is.character(guse), length(guse) == 1, !is.null(x[[guse]]))
        }
        cs <- x[[gid]] == guse
        x <- x[, cs, drop = FALSE]
    }
    # keep cells with ≥ 'min_genes' detected;
    # keep genes with ≥ 'min_count's in ≥ 'min_cells'
    gs <- rowSums(counts(x) > min_count) >= min_cells
    cs <- colSums(counts(x) > 0) >= min_genes 
    x <- x[gs, cs, drop = FALSE]
    
    # keep cluster-samples with ≥ 'min_size' cells
    nc <- table(x[[kid]], x[[sid]])
    nc <- .filter_matrix(nc, n = min_size)
    cs1 <- x[[kid]] %in% rownames(nc)
    cs2 <- x[[sid]] %in% colnames(nc)
    x <- x[, cs1 & cs2]
   
    # filtering done; drop no longer existent factors
    x <- .update_sce(x)
    
    # construct design matrix
    sep <- ifelse(is.null(kid) | is.null(sid), "", ":")
    f <- paste(kid, sid, sep = sep)
    if (length(f) == 0) {
        mm <- NULL 
    } else {
        f <- sprintf("~0+I(%s)", f)
        cd <- as.data.frame(colData(x))
        mm <- model.matrix(as.formula(f), data = cd)
    }
    
    # estimate NB parameters
    y <- DGEList(counts(x))
    y <- estimateDisp(y, mm)
    fit <- glmFit(y, prior.count = 0)
    
    # fit cell offsets
    if (is.null(mm)) {
        by <-  rep(1, ncol(x))
        depth <- 1
    } else {
        by <- c(kid, sid)
        depth <- length(by)
    }
    cs <- as.character(seq_len(ncol(x)))
    os <- set_names(c(fit$offset), cs)
    dt <- data.table(cd, cid = cs)
    dt_split <- split(dt, by, sorted = TRUE, flatten = FALSE)
    cs_split <- map_depth(dt_split, depth, "cid")
    fits <- map_depth(cs_split, depth, function(cs)
        as.list(fitdistr(os[cs], "normal")$estimate))
    
    # write average logCPM, tagwise & trended 
    # dispersion estimates to gene metadata
    vars <- c("AveLogCPM", "tagwise.dispersion", "trended.dispersion")
    pars <- vapply(vars, function(u) y[[u]], numeric(nrow(x)))
    colnames(pars) <- c("mean_logCPM", "tagwise_disp", "trended_disp")
    rowData(x) <- cbind(rowData(x), pars)
    metadata(x)$common_disp <- y$common.dispersion
    
    # write beta estimates to gene metadata
    bs <- fit$coefficients
    bs <- split(bs, row(bs))
    bs <- lapply(bs, matrix, ns, nk)
    bs <- lapply(bs, array, c(ns, nk, 1), list(sids, kids))
    rowData(x)$beta <- bs

    # write offset estimates to cell metadata
    
    return(x)
}
data(sce); x <- sce[seq_len(1e3), ]
colnames(colData(x)) <- c("a", "b", "c")
colData(x) <- colData(x)[, -3]
x <- simPars(x, "a", "b")

data(sce); x <- sce[seq_len(1e3), ]
x <- x[, x$group_id == "ctrl"]
colData(x) <- colData(x)[, -3]
x$sample_id <- droplevels(x$sample_id)

nk <- length(kids <- levels(x$cluster_id))
ns <- length(sids <- levels(x$sample_id))

y <- DGEList(counts(x))
cd <- as.data.frame(colData(x))
mm <- model.matrix(data = cd, ~ 0+I(cluster_id))
colnames(mm)

y <- estimateDisp(y, mm)
y <- glmFit(y, prior.count = 0)
head(y$coefficients)

bs <- y$coefficients
bs <- split(bs, row(bs))
bs <- lapply(bs, matrix, ns, nk)
bs <- lapply(bs, array, c(ns, nk, 1), list(sids, kids))
rowData(x)$beta <- bs
head(rowData(x))
rowData(x)$beta[[1]]

bs <- matrix(bs, nrow(x), nk, ns)
head(bs)
bs <- split(bs, col(bs))
bs <- matrix(bs, nrow(x), nk, ns, , dimnames = list(kids, sids))

names(gs) <- gs <- rownames(x)
bs <- sapply(gs, function(g)
    list(matrix(y$coefficients[g, ], ns, nk, 
        dimnames = list(sids, kids))))
rowData(x)$beta <- bs
rowData(x)$beta[[1]]
head(rowData(x))

x$offset <- c(y$offset)
rowData(x)$dispersion <- y$dispersion
betas <- paste0("beta.", levels(x$sample_id))
colnames(y$coefficients) <- betas
rowData(x)[, betas] <- y$coefficients



bs <- split(bs, col(bs))
matrix(bs, nlevels(x$cluster_id), nlevels(x$sample_id), nrow(x))
x <- ref

ng <- 500
nc <- 2e3
nk <- ns <- 3
kp <- sp <- gp <- NULL
pdd <- c(0.2, 0.2, rep(0.15, 4))
lfc <- list(dist = "gamma", shape = 4, rate = 2)
paired <- FALSE

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

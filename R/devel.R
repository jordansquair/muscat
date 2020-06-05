# source("R/estPars.R")
# source("R/utils.R")
# suppressPackageStartupMessages({
#     library(edgeR)
#     library(MASS)
#     library(SingleCellExperiment)
# })
# data(sce); x <- sce; x0 <- estPars(x)
# x <- y
# ng <- 500
# nc <- 2e3
# nk <- ns <- 3
# kp <- sp <- gp <- NULL
# pdd <- c(0.2, 0.2, rep(0.15, 4))
# lfc <- list(dist = "gamma", shape = 4, rate = 2)
# paired <- FALSE
# cats <- muscat:::cats

simData2 <- function(x, 
    nc = 2e3, ng = 1e3, nk = 4, ns = 3, nb = 2,
    kp = NULL, sp = NULL, gp = NULL, bp = NULL, paired = FALSE,
    pdd = diag(6)[1, ], lfc = list(family = "gamma", shape = 4, rate = 2)) {

    # ng <- nrow(x)
    # nk <- nlevels(x$cluster_id)
    # ns <- nlevels(x$sample_id)
    
    # throughout this code...
    # k = cluster ID
    # s = sample ID
    # g = group ID
    # c = DD category
    # 0 = reference
    
    args <- as.list(environment())

    # get reference IDs
    names(kids0) <- kids0 <- levels(x$cluster_id)
    names(sids0) <- sids0 <- levels(x$sample_id)
    names(bids0) <- bids0 <- levels(x$batch_id)
    nk0 <- ifelse((nk0 <- length(kids0)) == 0, 1, nk0)
    ns0 <- ifelse((ns0 <- length(sids0)) == 0, 1, ns0)
    nb0 <- ifelse((nb0 <- length(bids0)) == 0, 1, nb0)
    
    # set simulation IDs & gene, cell names
    names(kids) <- kids <- paste0("cluster", seq_len(nk))
    names(sids) <- sids <- paste0("sample", seq_len(ns))
    names(bids) <- bids <- paste0("batch", seq_len(nb))
    names(gids) <- gids <- paste0("group", seq_len(2))
    gs <- paste0("gene", seq_len(ng))
    cs <- paste0("cell", seq_len(nc))
    
    # default to equal probability for each cluster, sample, group ID
    if (is.null(kp)) kp <- rep(1/nk, nk)
    if (is.null(sp)) sp <- rep(1/ns, ns)
    if (is.null(bp)) bp <- rep(1/nb, nb)
    if (is.null(gp)) gp <- rep(0.5, 2)
    
    # sample cell metadata
    cd <- data.frame(
        cluster_id = sample(kids, nc, TRUE, kp),
        sample_id = sample(sids, nc, TRUE, sp),
        batch_id = sample(bids, nc, TRUE, bp),
        group_id = sample(gids, nc, TRUE, gp),
        row.names = cs, stringsAsFactors = FALSE)
    
    # sample number of genes to simulate per DD category
    ndd_ck <- vapply(kids, function(k) 
        table(sample(cats, ng, TRUE, pdd)), 
        numeric(length(cats)))

    # get sim_gene indices for each cluster-category
    gi_ck <- vapply(kids, function(k) 
        split(sample(gs), rep.int(cats, ndd_ck[, k])), 
        vector("list", length(cats)))
    
    # sample reference genes to simulate from 
    gs0 <- set_names(sample(rownames(x), ng, ng > nrow(x)), gs)
    
    # split by cluster & DD categroy
    gs0 <- replicate(nk, gs0, simplify = FALSE)
    gs0 <- set_names(gs0, kids)
    gs0_ck <- vapply(kids, function(k) {
        lapply(unfactor(cats), function(c) 
            gs0[[k]][gi_ck[[c, k]]])
    }, vector("list", length(cats)))
   
    # sample reference clusters & samples
    #bids0_use <- set_names(sample(bids0, nb, nb > nb0), bids)
    #kids0_use <- setNames(sample(kids0, nk, nk > nk0), kids)
    bids0_use <- set_names(rep(sample(bids0, 1), nb), bids)
    kids0_use <- set_names(rep(sample(kids0, 1), nk), kids)
    if (paired) { 
        # use same set of reference samples for all groups
        sids0_use <- sample(sids0, ns, ns > ns0)
        sids0_use <- replicate(length(gids), sids0_use)
    } else {
        # draw reference samples at random for each group
        sids0_use <- replicate(length(gids), sample(sids0, ns, ns > ns0))
    }
    sids0_use <- rep(sample(sids0, 1), ns)
    sids0_use <- replicate(2, sids0_use)
    dimnames(sids0_use) <- list(sids, gids)
  
    # sample logFCs for each cluster-category
    dist <- get(paste0("r", lfc$family))
    pars <- lfc[-grep("family", names(lfc))]
    lfcs_ck <- vapply(kids, function(k)
        lapply(unfactor(cats), function(c) {
            n <- ndd_ck[c, k]
            if (c == "ee")
                return(rep(NA, n))
            # sample directions for each gene with 
            # equal probability of up- & down-regulation
            dirs <- sample(c(-1, 1), n, TRUE) 
            # sample logFCs from input distribution
            set_names(do.call(dist, c(n, pars)), gi_ck[[c, k]])
        }), vector("list", length(cats)))
   
    # sample cell indices for each cluster-sample-batch-group
    dt <- data.table(cell_id = rownames(cd), cd)
    ci_ksb <- split(dt, by = colnames(cd), sorted = TRUE, flatten = FALSE)
    ci_ksb <- lapply(kids, function(k)
        vapply(bids, function(b) {
         lapply(sids, function(s) 
             map(ci_ksb[[k]][[s]][[b]], "cell_id"))
        }, vector("list", ns)))

    # # split reference cells by cluster-sample
    # cd0 <- as.data.frame(colData(x))
    # dt0 <- data.table(cell_id = colnames(x), cd0)
    # dt0 <- split(dt0, 
    #     sorted = TRUE, flatten = FALSE,
    #     by = c("cluster_id", "sample_id"))
    # cs0 <- map_depth(dt0, 2, "cell_id")
    
    # get simulation parameters
    os0 <- x$offset
    os_fit <- fitdistr(os0, "normal")
    os_fit <- as.list(os_fit$estimate)
    bs0 <- metadata(x)$betas
    ds0 <- rowData(x)$trended.dispersion
    names(ds0) <- rownames(x)

    # sample counts for each cluster, sample, batch, category ------------------
    res <- mapply(function(k, s, b) {
        # get reference cluster, samples
        k0 <- kids0_use[[k]]
        s0 <- sids0_use[s, ]
        b0 <- bids0_use[b]
        #cs0_ks <- set_names(cs0[[k0]][s0], gids)
        
        # get sim cell indices & 
        # number of cells to simulate per group
        ci <- ci_ksb[[k]][[s, b]]
        ncs <- vapply(ci, length, numeric(1))
        
        # sample offsets
        os <- lapply(ncs, function(n)
            do.call(rnorm, c(n = n, os_fit)))
        
        # get betas
        bs_k <- bs0$k[[ 1, k0]]
        bs_b <- bs0$b[[k0, b0]]
        bs_s <- set_names(bs0$s[k0, s0], gids)
        # sum across cluster-sample-batch
        bs_g <- lapply(gids, function(g)
            rowSums(cbind(bs_k, bs_s[[g]], bs_b)))
        
        res <- lapply(unfactor(cats[ndd_ck[, k] != 0]), function(c) {
            # get ref genes & sim gene indices
            gs0 <- gs0_ck[[c, k]]
            gi <- gi_ck[[c, k]]
            # get dispersion
            ds <- set_names(ds0[gs0], gi)
            # compute NB means
            ms <- lapply(gids, function(g) {
                bs <- bs_g[[g]][gs0]
                ms <- outer(exp(bs), exp(os[[g]])) 
                rownames(ms) <- gi
                colnames(ms) <- ci[[g]]
                return(ms)
            })
            # sample counts
            res <- .sim(
                cat = c, cs = ci, ms, ds, 
                lfcs = lfcs_ck[[c, k]], ep, dp, dm)
            res$cs[gi, unlist(ci)]
        })
        do.call(rbind, res)[gs, ]
    }, 
        k = rep(kids, each = ns*nb), 
        s = rep(sids, nk*nb), 
        b = rep(bids, each = ns), 
        SIMPLIFY = FALSE)
    y <- do.call(cbind, res)[, cs]
    
    # construct gene metadata table storing ------------------------------------
    # gene | cluster_id | category | logFC | ref. gene | disp | gorup means
    md <- data.frame(
        row.names = NULL,
        gene = unlist(gi_ck),
        cluster_id = rep.int(rep(kids, each = length(cats)), c(ndd_ck)),
        category = rep.int(rep(cats, nk), c(ndd_ck)),
        logFC = unlist(lfcs_ck),
        ref_gene = unlist(gs0))
    o <- order(as.numeric(gsub("[a-z]", "", md$gene)))
    md <- md[o, ]; rownames(md) <- NULL
    
    # construct SCE ------------------------------------------------------------
    # cell metadata including group, sample, cluster IDs
    cd <- mutate_all(cd, as.factor)
    cd$group_id <- droplevels(cd$group_id)
    cd <- mutate(cd, unique_id = factor(paste(
        sample_id, batch_id, group_id, sep = ".")), .before = 1)
    ids <- levels(cd$unique_id)
    m <- match(ids, cd$unique_id)
    o <- order(gids <- cd$group_id[m])
    ei <- select(cd[m, ], -cluster_id)
    rownames(ei) <- NULL
    # gene metadata storing gene classes & specificities
    rd <- NULL
    # simulation metadata including reference 
    # clusters/samples, gene metadata, function call
    md <- list(
        experiment_info = ei,
        n_cells = table(cd$unique_id),
        gene_info = md,
        ref_kids = kids0_use,
        ref_sids = sids0_use, 
        args = args)
    # return SCE
    SingleCellExperiment(
        assays = list(counts = as.matrix(y)),
        colData = cd, rowData = rd, metadata = md)
}

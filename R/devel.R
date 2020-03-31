# source("R/estPars.R")
# source("R/utils.R")
# suppressPackageStartupMessages({
#     library(edgeR)
#     library(MASS)
#     library(SingleCellExperiment)
# })
# data(sce); x <- sce; x0 <- estPars(x)
# x <- x0
# ng <- 500
# nc <- 2e3
# nk <- ns <- 3
# kp <- sp <- gp <- NULL
# pdd <- c(0.2, 0.2, rep(0.15, 4))
# lfc <- list(dist = "gamma", shape = 4, rate = 2)
# paired <- FALSE
# cats <- muscat:::cats

simData <- function(x, 
    nc = 2e3, #ng = 1e3, nk = 4, ns = 3,
    kp = NULL, sp = NULL, gp = NULL, paired = TRUE,
    pdd = diag(6)[1, ], lfc = list(family = "gamma", shape = 4, rate = 2)) {

    ng <- nrow(x)
    nk <- nlevels(x$cluster_id)
    ns <- nlevels(x$sample_id)
    
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
    ci <- vapply(sids, function(s) {
        lapply(kids, function(k) 
            map(ci[[k]][[s]], "cell_id"))
    }, vector("list", nk))
    
    # split reference cells by cluster-sample
    cd0 <- as.data.frame(colData(x))
    dt0 <- data.table(cell_id = colnames(x), cd0)
    dt0 <- split(dt0, 
        sorted = TRUE, flatten = FALSE,
        by = c("cluster_id", "sample_id"))
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
    
    # split by cluster & categroy
    gs0 <- replicate(nk, gs0)
    gs0 <- split(gs0, col(gs0))
    gs0 <- set_names(map(gs0, set_names, gs), kids)
    gs0 <- vapply(kids, function(k) {
        lapply(unfactor(cats), function(c) 
            gs0[[k]][gi[[c, k]]])
    }, vector("list", length(cats)))
    
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
    dist <- get(paste0("r", lfc$family))
    pars <- lfc[-grep("family", names(lfc))]
    lfcs <- vapply(kids, function(k)
        lapply(unfactor(cats), function(c) {
            n <- ndd[c, k]
            if (c == "ee")
                return(rep(NA, n))
            # sample directions for each gene with 
            # equal probability of up- & down-regulation
            dirs <- sample(c(-1, 1), n, TRUE) 
            # sample logFCs from input distribution
            set_names(do.call(dist, c(n, pars)), gs0[[c, k]])
        }), vector("list", length(cats)))
    
    # sample counts for each cluster, sample, category -------------------------
    res <- mapply(function(k, s) {
        # get reference cluster, samples, cells
        k0 <- kids0_use[[k]]
        s0 <- sids0_use[s, ]
        cs0_ks <- set_names(cs0[[k0]][s0], gids)
        
        # get simulation cell indices & 
        # number of cells to simulate per group
        ci_ks <- ci[[k, s]]
        ncs <- lapply(ci_ks, length)
        
        res <- lapply(unfactor(cats[ndd[, k] != 0]), function(c) {
            # get reference genes & simulation gene indices
            gs0_kc <- gs0[[c, k]]
            gi_kc <- gi[[c, k]]
            # sample NB parameters
            ds <- set_names(rowData(x)[gs0_kc, "tagwise.dispersion"], gi_kc)
            bs <- set_names(metadata(x)$betas[k0, s0], gids)
            os <- set_names(metadata(x)$offsets[k0, s0], gids)
            mus <- map(os, "mean")
            sds <- map(os, "sd")
            ms <- lapply(gids, function(g) {
                bs <- bs[[g]][gs0_kc]
                os <- rnorm(ncs[[g]], mus[[g]], sds[[g]])
                outer(exp(bs), exp(os)) %>% 
                    set_rownames(gi_kc) %>% 
                    set_colnames(ci_ks[[g]])
            })
            # sample counts
            res <- .sim(
                cat = c, cs = ci_ks, ms, ds, 
                lfcs = lfcs[[c, k]], ep, dp, dm)
        })
        counts <- map_depth(res, 1, "cs") %>% 
            lapply("[", i = TRUE, j = unlist(ci_ks)) %>% 
            do.call(what = "rbind") %>% `[`(i = gs, j = TRUE)
    }, k = rep(kids, each = ns), s = rep(sids, nk), SIMPLIFY = FALSE)
    y <- do.call("cbind", res)[, cs]
    
    # construct gene metadata table storing ------------------------------------
    # gene | cluster_id | category | logFC | ref. gene | disp | gorup means
    md <- data.frame(
        gene = unlist(gi),
        cluster_id = rep.int(rep(kids, each = length(cats)), c(ndd)),
        category = rep.int(rep(cats, nk), c(ndd)),
        logFC = unlist(lfcs),
        ref_gene = unlist(gs0))
    o <- order(as.numeric(gsub("[a-z]", "", md$gene)))
    md <- md[o, ]; rownames(gi) <- NULL
    
    # construct SCE ------------------------------------------------------------
    # cell metadata including group, sample, cluster IDs
    cd <- mutate_all(cd, as.factor)
    cd$group_id <- droplevels(cd$group_id)
    cd$sample_id <- factor(paste(cd$group_id, cd$sample_id, sep = "."))
    sids <- levels(cd$sample_id)
    m <- match(sids, cd$sample_id)
    o <- order(gids <- cd$group_id[m])
    ei <- data.frame(sample_id = sids[o], group_id = gids[o])
    # gene metadata storing gene classes & specificities
    rd <- NULL
    # simulation metadata including reference 
    # clusters/samples, gene metadata, function call
    md <- list(
        experiment_info = ei,
        n_cells = table(cd$sample_id),
        gene_info = gi,
        ref_kids = kids0_use,
        ref_sids = sids0_use, 
        args = args)
    # return SCE
    SingleCellExperiment(
        assays = list(counts = as.matrix(y)),
        colData = cd, rowData = rd, metadata = md)
}

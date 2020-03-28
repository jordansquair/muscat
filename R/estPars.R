#' @param nc,ng,ns,nk number of cells, genes, samples, clusters to simulate. 
#' @param sp,kp,gp probability of cells to be in each sample, cluster, group.
# simData2 <- function(x, ng, nc, nk, ns, kp = NULL, sp = NULL, gp = NULL,
#     pdd = diag(6)[1, ], lfc = list(dist = "gamma", shape = 4, rate = 2),
#     paired = FALSE) {

estPars <- function(x,
    kid = "cluster_id", sid = "sample_id", gid = NULL, guse = NULL,
    min_count = 1, min_cells = 10, min_genes = 100, min_size = 100,
    verbose = TRUE) {

    # store list of input arguments
    args <- as.list(environment())
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"),
        is.null(kid) || is.character(kid) & length(kid) == 1 & !is.null(x[[kid]]),
        is.null(sid) || is.character(sid) & length(sid) == 1 & !is.null(x[[sid]]),
        is.null(gid) || is.character(gid) & length(gid) == 1 & !is.null(x[[gid]]),
        is.logical(verbose), length(verbose) == 1)

    # data(sce); x <- sce
    # args <- list(kid = NULL, NULL)
    # for (id in names(args)) assign(id, args[[id]])
    # min_count <- 1; min_cells = 10; min_genes = 100; min_size = 100; guse <- NULL

    # convert cell metadata ID columns to factors
    ids <- grep("id", names(args), value = TRUE)
    rmv <- vapply(args[ids], is.null, logical(1))
    for (id in args[ids[!rmv]]) 
        x[[id]] <- factor(x[[id]])

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
    if (!is.null(c(kid, sid))) {
        id <- which(c(!is.null(kid), !is.null(sid)))
        if (length(id) == 1) {
            id <- switch(paste(id), "1" = kid, "2" = sid) 
            ns <- table(x[[id]])
            ids <- levels(x[[id]])
            ids <- ids[ns > min_cells]
            cs <- x[[id]] %in% ids
        } else {
            nc <- table(x[[kid]], x[[sid]])
            nc <- .filter_matrix(nc, n = min_size)
            cs1 <- x[[kid]] %in% rownames(nc)
            cs2 <- x[[sid]] %in% colnames(nc)
            cs <- cs1 & cs2
           
        }
        x <- x[, cs]
    }

    # filtering done; drop no longer existent factors
    x <- .update_sce(x)
    
    # get final cell metadata
    cd <- as.data.frame(colData(x))
    if (!is.null(kid)) nk <- length(kids <- set_names(levels(x[[kid]])))
    if (!is.null(sid)) ns <- length(sids <- set_names(levels(x[[sid]])))
    
    # construct design matrix
    sep <- ifelse(is.null(kid) | is.null(sid), "", ":")
    f <- paste(kid, sid, sep = sep)
    if (length(f) == 0) {
        mm <- NULL 
    } else {
        f <- sprintf("~0+I(%s)", f)
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
    cs <- colnames(x)
    os <- set_names(c(fit$offset), cs)
    if (!is.null(c(kid, sid))) {
        dt <- data.table(cd, cid = cs)
        dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
        cs_split <- map_depth(dt_split, depth, "cid")
        musd <- map_depth(cs_split, depth, function(cs) 
            as.list(fitdistr(os[cs], "normal")$estimate))
    } else {
        musd <- set_names(as.list(fitdistr(os, "normal")$estimate), cs)
    }
    
    # write average logCPM, tagwise & trended 
    # dispersion estimates to gene metadata
    vars <- c("AveLogCPM", "tagwise.dispersion", "trended.dispersion")
    pars <- vapply(set_names(vars), function(u) y[[u]], numeric(nrow(x)))
    rowData(x) <- cbind(rowData(x), pars)
    
    # write beta estimates to gene metadata
    bs <- fit$coefficients
    bs <- split(bs, row(bs))
    bs <- lapply(bs, matrix, ns, nk)
    bs <- lapply(bs, array, c(ns, nk, 1), list(sids, kids))
    rowData(x)$beta <- bs
    
    # write common dispersion & offset parameters to metadata
    metadata(x)$common_disp <- y$common.dispersion
    metadata(x)$offsets <- musd
    metadata(x)$args <- args
    
    # return SCE
    return(x)
}
data(sce); x0 <- sce
x1 <- estPars(x0, sid = NULL, gid = "group_id")
plotMeanDisp(x1)

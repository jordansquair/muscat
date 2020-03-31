#' @param nc,ng,ns,nk number of cells, genes, samples, clusters to simulate. 
#' @param sp,kp,gp probability of cells to be in each sample, cluster, group.
# simData2 <- function(x, ng, nc, nk, ns, kp = NULL, sp = NULL, gp = NULL,
#     pdd = diag(6)[1, ], lfc = list(dist = "gamma", shape = 4, rate = 2),
#     paired = FALSE) {

estPars <- function(x,
    kid = "cluster_id", sid = "sample_id", gid = "group_id", guse = NULL,
    min_count = 1, min_cells = 10, min_genes = 100, min_size = 100, verbose = TRUE) {

    # store list of input arguments
    args <- as.list(environment())
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"),
        length(grep("^counts$", assayNames(x))) == 1, 
        is.null(kid) || is.character(kid) & length(kid) == 1 & !is.null(x[[kid]]),
        is.null(sid) || is.character(sid) & length(sid) == 1 & !is.null(x[[sid]]),
        is.null(gid) || is.character(gid) & length(gid) == 1 & !is.null(x[[gid]]),
        is.null(guse) || is.character(guse) & length(guse) == 1 & !is.null(x[[guse]]),
        is.numeric(min_count), length(min_count) == 1, 
        is.numeric(min_cells), length(min_cells) == 1, 
        is.numeric(min_genes), length(min_genes) == 1, 
        is.numeric(min_size), length(min_size) == 1, 
        is.logical(verbose), length(verbose) == 1)

    # data(sce); x <- sce
    # names(colData(x)) <- c("abc", "x", "y")
    # args <- list(kid = "abc", sid = NULL, gid = NULL)
    # for (id in names(args)) assign(id, args[[id]])
    # min_count <- 1; min_cells = 10; min_genes = 100; min_size = 100; guse <- NULL

    if (is.null(rownames(x))) rownames(x) <- paste0("ref_gene", seq_len(nrow(x)))
    if (is.null(colnames(x))) colnames(x) <- paste0("ref_cell", seq_len(ncol(x)))
    
    # convert cell metadata ID columns to factors
    ids <- grep("id", names(args), value = TRUE)
    defs <- unlist(formals("estPars")[ids])
    for (id in ids) {
        if (is.null(args[[id]])) {
            def <- defs[id]
            x[[def]] <- NULL
            x[[def]] <- paste0(gsub("_id$", "", def), "1")
        } else {
            m <- match(args[[id]], names(colData(x)))
            def <- defs[id]
            names(colData(x))[m] <- def
            x[[def]] <- factor(x[[def]])
        }
    }
    # drop all other cell metadata
    colData(x) <- colData(x)[, defs]
    for (id in defs) x[[id]] <- factor(x[[id]])

    # default to using first group available
    if (is.null(guse)) 
        guse <- levels(x[[gid]])[1] 
    cs <- x[[gid]] == guse
    x <- x[, cs]
    
    # keep cells with ≥ 'min_genes' detected;
    # keep genes with ≥ 'min_count's in ≥ 'min_cells'
    gs <- rowSums(counts(x) > min_count) >= min_cells
    cs <- colSums(counts(x) > 0) >= min_genes 
    x <- x[gs, cs]

    # keep cluster-samples with ≥ 'min_size' cells
    if (!is.null(c(kid, sid))) {
        nc <- table(x$cluster_id, x$sample_id)
        nc <- .filter_matrix(nc, n = min_size)
        cs1 <- x$cluster_id %in% rownames(nc)
        cs2 <- x$sample_id %in% colnames(nc)
        x <- x[, cs1 & cs2]
    }

    # filtering done; drop no longer existent factors
    x <- .update_sce(x)
    
    # get final cell metadata
    cd <- as.data.frame(colData(x))
    nk <- length(kids <- purrr::set_names(levels(x$cluster_id)))
    ns <- length(sids <- purrr::set_names(levels(x$sample_id)))
    
    # construct design matrix
    f <- paste( 
        ifelse(nk == 1, "", "cluster_id"), 
        ifelse(ns == 1, "", "sample_id"),
        sep = ifelse(nk == 1 | ns == 1, "", ":"))
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
    # format into 'nk' x 'ns' matrix
    bs <- fit$coefficients
    bs <- split(bs, col(bs))
    bs <- map(bs, set_names, rownames(x))
    bs <- matrix(bs, nk, ns, TRUE, list(kids, sids))
    
    # fit cell offsets by cluster-sample
    cs <- colnames(x)
    os <- set_names(c(fit$offset), cs)
    cs_by_ks <- data.table(cd, cid = cs) %>% 
        split(by = c("cluster_id", "sample_id"), 
            sorted = TRUE, flatten = FALSE) %>% 
        map_depth(2, "cid")
    os <- vapply(sids, function(s) 
        lapply(kids, function(k) {
            cs <- cs_by_ks[[k]][[s]]
            fit <- fitdistr(os[cs], "normal")
            as.list(fit$estimate)
        }), vector("list", nk))
    
    # write average logCPM, tagwise & trended 
    # dispersion estimates to gene metadata
    vars <- c("AveLogCPM", "tagwise.dispersion", "trended.dispersion")
    pars <- vapply(purrr::set_names(vars), function(u) y[[u]], numeric(nrow(x)))
    rowData(x) <- cbind(rowData(x), pars)
    
    # store common dispersion, betas, 
    # offset fits & function in metadata
    metadata(x) <- list(
        common.dispersion = y$common.dispersion,
        betas = bs, offsets = os, args = args[-1])
    
    # return SCE
    return(x)
}
# data(sce); x0 <- sce
# x1 <- estPars(x0)
# plotMeanDisp(x1)
# 
# metadata(x1)
# metadata(x1)$betas
# metadata(x1)$offsets

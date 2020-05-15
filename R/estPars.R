#' @param nc,ng,ns,nk number of cells, genes, samples, clusters to simulate. 
#' @param sp,kp,gp probability of cells to be in each sample, cluster, group.
# simData2 <- function(x, ng, nc, nk, ns, kp = NULL, sp = NULL, gp = NULL,
#     pdd = diag(6)[1, ], lfc = list(dist = "gamma", shape = 4, rate = 2),
#     paired = FALSE) {

data(sce)
x <- sce
x$batch_id <- sample(x$group_id)
x <- x[, x$group_id == "ctrl"]
x$group_id <- factor("a")
kid <- "cluster_id"
sid <- "sample_id"
bid <- "batch_id"
gid <- "group_id"
guse <- NULL
min_count <- 1
min_cells <- 20
min_genes <- 100
min_size <- 100
verbose <- TRUE

estPars <- function(x,
    kid = "cluster_id", sid = "sample_id", bid = "batch_id", gid = "group_id", guse = NULL,
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

    #if (is.null(rownames(x))) rownames(x) <- paste0("ref_gene", seq_len(nrow(x)))
    #if (is.null(colnames(x))) colnames(x) <- paste0("ref_cell", seq_len(ncol(x)))
    
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
    x <- x[, x[[gid]] == guse]
    
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
    nk <- length(names(kids) <- kids <- levels(x[[kid]]))
    ns <- length(names(sids) <- sids <- levels(x[[sid]]))
    nb <- length(names(bids) <- bids <- levels(x[[bid]]))
    
    # construct design matrix
    # f <- paste(
    #     ifelse(nk == 1, "", "cluster_id"),
    #     ifelse(ns == 1, "", "sample_id"),
    #     sep = ifelse(nk == 1 || ns == 1, "", ":"))
    # if (length(f) == 0) {
    #     mm <- NULL
    # } else {
    #     df <- mutate(cd, foo = paste(cluster_id, sample_id, sep = "."))
    #     f <- as.formula(paste("~0+foo"))
    #     mm <- model.matrix(f, data = df)
    # }
    
    #df <- mutate(cd, foo = paste(cluster_id, sample_id, sep = "."))
    f <- as.formula(paste("~0+cluster_id+cluster_id:sample_id+cluster_id:batch_id"))
    mm <- model.matrix(f, data = cd)
    
    # estimate NB parameters
    y <- DGEList(counts(x))
    y <- estimateDisp(y, mm)
    fit <- glmFit(y, prior.count = 0)
    x$offset <- c(fit$offset)
    
    # format betas into 'nk' x 'ns' matrix
    bs <- fit$coefficients
    bs <- split(bs, col(bs))
    bs <- map(bs, set_names, rownames(x))
    names(bs) <- colnames(mm)
    
    # split by cluster, sample, batch effect
    s <- grep("sample_id", names(bs))
    b <- grep("batch_id", names(bs))
    k <- seq_len(ncol(mm))[-c(s, b)]
    bs_k <- matrix(bs[k],  1, nk, TRUE, list(NULL, kids))
    bs_s <- matrix(bs[s], nk, ns, TRUE, list(kids, sids))
    bs_b <- matrix(bs[b], nk, ns, TRUE, list(kids, bids))

    # write average logCPM, tagwise & trended 
    # dispersion estimates to gene metadata
    names(vars) <- vars <- c("AveLogCPM", "tagwise.dispersion", "trended.dispersion")
    pars <- vapply(vars, function(u) y[[u]], numeric(nrow(x)))
    rowData(x) <- cbind(rowData(x), pars)
    
    # store common dispersion, betas, 
    # offset fits & function in metadata
    metadata(x) <- list(args = args[-1],
        common.dispersion = y$common.dispersion,
        betas = list(k = bs_k, s = bs_s, b = bs_b))
    
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

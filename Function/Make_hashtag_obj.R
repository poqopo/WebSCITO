make_hashtag_obj <-function(raw, abNum, batchNum) {
  tx <- raw$`Gene Expression`[, which(Matrix::colSums(raw$`Gene Expression`) > 0)]
  ab <- raw$`Antibody Capture`[, which(Matrix::colSums(raw$`Antibody Capture`) > 0)]
  ccd_indices <- intersect(which(Matrix::colSums(tx) > 0), which(Matrix::colSums(ab) > 0))
  
  ab <- ab[, ccd_indices]
  tx <- tx[, ccd_indices]
  
  colnames(ab) <- make.names(colnames(ab), unique = TRUE)
  colnames(tx) <- make.names(colnames(tx), unique = TRUE)
  
  ab_hashtag <- data.frame()
  
  # 반복문 실행
  for (i in 1:as.integer(batchNum)) {
    start <- (1 + abNum*(i-1))
    end <- i*abNum
    result <- round(colSums(t(t(ab[start:end,])/mean(ab[start,]))*100),0)
    ab_hashtag <- rbind(ab_hashtag, result)
  }
  
  # ab_hashtag의 열 이름 설정
  colnames(ab_hashtag) <- colnames(ab)
  rownames(ab_hashtag) <- paste("Batch", seq(1:as.integer(batchNum)), sep = "")
  ab_hashtag <- as.matrix(ab_hashtag)
  
  MaxN <- function(x, N = 2){
    len <- length(x)
    if (N > len) {
      warning('N greater than length(x).  Setting N=length(x)')
      N <- length(x)
    }
    sort(x, partial = len - N + 1)[len - N + 1]
  }
  
  HTODemux <- function (object, assay = "HTO", positive.quantile = 0.99, init = NULL, 
                        nstarts = 100, kfunc = "clara", nsamples = 100, seed = 42, 
                        verbose = TRUE) {
    if (!is.null(x = seed)) {
      set.seed(seed = seed)
    }
    assay <- assay %||% DefaultAssay(object = object)
    data <- GetAssayData(object = object, assay = assay)
    counts <- GetAssayData(object = object, assay = assay, slot = "counts")[, 
                                                                            colnames(x = object)]
    counts <- as.matrix(x = counts)
    ncenters <- init %||% (nrow(x = data) + 1)
    switch(EXPR = kfunc, kmeans = {
      init.clusters <- kmeans(x = t(x = GetAssayData(object = object, 
                                                     assay = assay)), centers = ncenters, nstart = nstarts)
      Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
    }, clara = {
      init.clusters <- clara(x = t(x = GetAssayData(object = object, 
                                                    assay = assay)), k = ncenters, samples = nsamples)
      Idents(object = object, cells = names(x = init.clusters$clustering), 
             drop = TRUE) <- init.clusters$clustering
    }, stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'"))
    average.expression <- AverageExpression(object = object, 
                                            assays = assay, verbose = FALSE)[[assay]]
    if (sum(average.expression == 0) > 0) {
      stop("Cells with zero counts exist as a cluster.")
    }
    discrete <- GetAssayData(object = object, assay = assay)
    discrete[discrete > 0] <- 0
    for (iter in rownames(x = data)) {
      values <- counts[iter, colnames(object)]
      values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, 
      ])]])]
      fit <- suppressWarnings(expr = fitdist(data = values.use, 
                                             distr = "nbinom"))
      cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
      discrete[iter, names(x = which(x = values > cutoff))] <- 1
      if (verbose) {
        message(paste0("Cutoff for ", iter, " : ", cutoff, 
                       " reads"))
      }
    }
    npositive <- Matrix::colSums(x = discrete)
    classification.global <- npositive
    classification.global[npositive == 0] <- "Negative"
    classification.global[npositive == 1] <- "Singlet"
    classification.global[npositive == 2] <- "Doublet"
    classification.global[npositive == 3] <- "Triplet"
    classification.global[npositive == 4] <- "Quadruplet"
    classification.global[npositive == 5] <- "Quintuplet"
    classification.global[npositive == 6] <- "Sextuplet"
    classification.global[npositive == 7] <- "Septuplet"
    donor.id = rownames(x = data)
    hash.max <- apply(X = data, MARGIN = 2, FUN = max)
    hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
    hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
    hash.maxID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data), 
                                                   FUN = function(x) {
                                                     return(which(x = data[, x] == hash.max[x])[1])
                                                   })])
    hash.secondID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data), 
                                                      FUN = function(x) {
                                                        return(which(x = data[, x] == hash.second[x])[1])
                                                      })])
    ##browser();
    hash.margin <- hash.max - hash.second
    doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
      return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                   collapse = "_"))
    })
    classification <- classification.global
    classification[classification.global == "Negative"] <- "Negative"
    classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == 
                                                                             "Singlet")]
    classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == 
                                                                             "Doublet")]
    classification.metadata <- data.frame(hash.maxID, hash.secondID, 
                                          hash.margin, classification, classification.global)
    colnames(x = classification.metadata) <- paste(assay, c("maxID", 
                                                            "secondID", "margin", "classification", "classification.global"), 
                                                   sep = "_")
    object <- AddMetaData(object = object, metadata = classification.metadata)
    Idents(object) <- paste0(assay, "_classification")
    doublets <- rownames(x = object[[]])[which(object[[paste0(assay, 
                                                              "_classification.global")]] == "Doublet")]
    Idents(object = object, cells = doublets) <- "Doublet"
    object$hash.ID <- Idents(object = object)
    ##object$discrete <- discrete;
    return(list(object,discrete))
  }
  
  
  ab_hashtag_obj <- CreateSeuratObject(counts = ab)
  ab_hashtag_obj[["HTO"]] <- CreateAssayObject(counts = ab_hashtag)
  ab_hashtag_obj <- NormalizeData(ab_hashtag_obj, assay = "HTO", normalization.method = "CLR")
  
  rst <- HTODemux(ab_hashtag_obj, assay = "HTO", positive.quantile = 0.99)
  
  ab_hashtag_obj <- rst[[1]]
  discrete <- rst[[2]]
  
  Idents(ab_hashtag_obj) <- "HTO_maxID"
  
  return(list(ab = ab,
              tx= tx,
              discrete = discrete,
              ab_hashtag_obj = ab_hashtag_obj))
}

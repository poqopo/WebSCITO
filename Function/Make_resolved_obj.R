make_resolved_obj <- function(ab, tx, aux, discrete, ab_hashtag_obj, abNum, batchNum) {
  ab_resolved_list <- list()
  tx_resolved_list <- list()
  
  for (i in 1:as.integer(batchNum)) {
    indices <- which(discrete[i, ] == 1)
    start <- (1 + abNum*(i-1))
    end <- i*abNum
    
    ab_subset <- ab[start:end, indices]
    ab_resolved_list[[i]] <- ab_subset
    
    tx_subset <- tx[, indices]
    tx_resolved_list[[i]] <- tx_subset
  }
  
  ab_resolved <- do.call(cbind, ab_resolved_list)
  tx_resolved <- do.call(cbind, tx_resolved_list)
  colnames(ab_resolved) <- make.names(colnames(ab_resolved), unique = TRUE)
  colnames(tx_resolved) <- make.names(colnames(tx_resolved), unique = TRUE)
  
  classification <- as.character(ab_hashtag_obj$HTO_classification.global)
  meta <- data.frame()
  aux <- NULL
  if (is.null(aux) == FALSE) {
    for (i in 1:as.integer(batchNum)) {
      indices <- which(discrete[i, ] == 1)
      
      batch_aux <- aux[indices, ]
      batch_meta <- data.frame(
        Batch = rep(paste0("Batch", i), sum(discrete[i, ])),
        Type = classification[indices]
      )
      batch_meta <- cbind(batch_meta, batch_aux)
      
      meta <- rbind(meta, batch_meta)
    }
  } else {
    for (i in 1:as.integer(batchNum)) {
      indices <- which(discrete[i, ] == 1)
      
      batch_meta <- data.frame(
        Batch = rep(paste0("Batch", i), sum(discrete[i, ])),
        Type = classification[indices]
      )
      
      meta <- rbind(meta, batch_meta)
    }
  }
  rownames(meta) <- colnames(ab_resolved)
  ab_resolved_obj <- CreateSeuratObject(counts = ab_resolved, meta.data = meta)
  ab_resolved_obj[["ADT"]] <- CreateAssay5Object(counts = ab_resolved)
  DefaultAssay(ab_resolved_obj) <- "ADT"
  Idents(ab_resolved_obj) <- "Batch"
  return(list(selected_ab_cell_matrix = as.data.frame(ab_resolved_obj[["ADT"]]$counts), meta_ab_cell_matrix = as.data.frame(ab_resolved_obj@meta.data)))
}

after_processing <-function(selected_ab_cell_matrix) {
  selected_ab_cell_matrix['CD235/61',] <- selected_ab_cell_matrix['CD235-barcode1',] + selected_ab_cell_matrix['CD61-barcode1',]
  selected_ab_cell_matrix <- selected_ab_cell_matrix[-c(2,3),]
  return(selected_ab_cell_matrix)
}
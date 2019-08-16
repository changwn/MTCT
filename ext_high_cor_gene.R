

ext_high_cor_gene <- function(pseudo_BULK, verbose=T)
{

  aaa <- cor(t(pseudo_BULK[[1]]), t(pseudo_BULK[[2]]))
  dim(aaa)
  high_cor_gene <- c()
  all_detail <- list()
  cell_type_cor <- c()
  for (i in 1:ncol(aaa))
  {
    if(verbose == TRUE)
    {
      print("another cell type---------")
      print(i)
      print(colnames(aaa)[i])
      print(sort(aaa[, i], decreasing = TRUE)[1:20])
    }
    high_cor_gene <- c(high_cor_gene, names(sort(aaa[, i], decreasing = TRUE)[1:20]))
    all_detail[[i]] <- sort(aaa[, i], decreasing = TRUE)[1:20]
    cell_type_cor<- c(cell_type_cor, mean(sort(aaa[, i], decreasing = TRUE)[1:20]))
  }
  names(all_detail) <- colnames(aaa)
  ave_cor <- mean(cell_type_cor)
  
  return(list(high_cor_gene, all_detail, cell_type_cor, ave_cor))
}




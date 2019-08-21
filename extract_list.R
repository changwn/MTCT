
extract_list <- function(mylist){
  nlength <- length(mylist)
  total <- c()
  for(i in 1:nlength){
    total <- c(total, mylist[[i]])
  }
  
  return(total)
}


check_pc <- function(my.df)
{
	gene_count <- nrow(my.df)
	independent_satisfy_count <- 0
	independent_flag <- TRUE
	for(i in 1:nrow(my.df))
	{
		if(my.df[i, 3] > 0.3 && my.df[i, 4] < 0.3){
			independent_satisfy_count <- independent_satisfy_count + 1 
		}
	}

	if(independent_satisfy_count < gene_count){
		independent_flag <- FALSE
	}

	return(independent_flag)

}


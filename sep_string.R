sep_string <- function(my.str, length_cut)
{
	tmp <- c()
	tmp <- unlist(strsplit(my.str, split=', '))
	if(length(tmp) > length_cut){
		tmp <- NULL
	}

	return(tmp)
}


pop_no_cor_pc <- function(CTES_row_bases, pass_id)
{
	pop_id <- c()
	remain_id <- setdiff(c(1:nrow(CTES_row_bases)), pass_id)
	for(i in remain_id)
	{
		count_test_independent <- 0
		for(j in pass_id)
		{
			if(cor(CTES_row_bases[i, ], CTES_row_bases[j, ]) < 0.1)
				count_test_independent <- count_test_independent + 1
				#print(cor(CTES_row_bases[i, ], CTES_row_bases[j, ]))
		}
		if(count_test_independent >= length(pass_id))
			pop_id<- c(pop_id, i)

	}

	return(pop_id)

}


pop_no_cor_pc_dynamic <- function(CTES_row_bases, pass_id)
{
	pop_id <- c()
	remain_id <- setdiff(c(1:nrow(CTES_row_bases)), pass_id)
	for(i in remain_id)
	{
		count_test_independent <- 0
		for(j in pass_id)
		{
			if(cor(CTES_row_bases[i, ], CTES_row_bases[j, ]) < 0.1)
				count_test_independent <- count_test_independent + 1
				#print(cor(CTES_row_bases[i, ], CTES_row_bases[j, ]))
		}
		if(count_test_independent >= length(pass_id))
			pass_id <- c(pass_id, i)
			#order is important! We just the first element once we detect the independence!

	}

	return(pass_id)

}

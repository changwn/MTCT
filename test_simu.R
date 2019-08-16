
setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_06_11_simulator')


#source 
source('scsimu_old_dirichlet.R')
source('scsimu_sep.R')
#load GSE720560 for testing
#load('C:/Users/wnchang/Documents/F/PhD_Research/2018_11_26_simulate_bulk/GSE72056_scRNA.RData')


ttt <- Bulk_Simu_no_coinfil_v1(GSE72056_scRNA, GSE72056_cell_key, cellNumber=1000, sampleNumber=50)
apply(ttt[[2]], 2,sum)

par(mfrow=c(2,4))
for(i in 1:8)
{
	hist(ttt[[2]][i, ], breaks=20, col='lightblue')
}

aaa <- rdirichlet(50, rep(1, 50))

ttt2 <- Bulk_Simu_coinfil(GSE72056_scRNA, GSE72056_cell_key,cellnumber=1000,'GSE72056',row1=1,row2=2,cor=0.3)


hist(ttt2[[3]][which(GSE72056_cell_key=='B_cell'),1], col='lightblue')
dim(ttt2[[2]])
ttt2[[2]][,1]
table(GSE72056_cell_key)

sum(ttt2[[3]][which(GSE72056_cell_key=='B_cell'),1])

library(ICTD)
prop <- ICTD(ttt2[[1]])[[1]]
cor(prop, ttt2[[2]])


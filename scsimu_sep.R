

#--------------------------
# 2019-08-09
# 
#
#--------------------------
# content
#Bulk_Simu_no_coinfil_v1
#Bulk_Simu_coinfil_v1
#mat_dirichlet












mat_dirichlet <- function(celltype_number, sample_number, centre){  
  if(!require('DirichletReg')) install.packages('DirichletReg')
  
  aaa <- rdirichlet(sample_number, rep(1, celltype_number))  
  aaa1 <- aaa * centre   #centre equal to 1000,which means sum of proportion of cell types equal to 100
  bbb <- t(aaa1)
  
  ccc <- rdirichlet(sample_number, rep(1, celltype_number))
  ccc1 <- ccc * centre 
  ddd <- t(ccc1)
  
  mat <- bbb + ddd
  mat <- mat * 0.5
  
  return(mat)
}

Bulk_Simu_no_coinfil_v1 <- function(sc_RNA, Cell_Type, cellNumber=1000, sampleNumber=100)
{
  print(Sys.time())
  print('Start simulation without co-infiltration!')
  BULK <- list()
  BULK[[1]] <- matrix(0, nrow = nrow(sc_RNA), ncol = sampleNumber)
  #BULK[[2]]<-matrix(0,ncol = sampleNumbe,nrow = length(unique(Cell_Type)))
  BULK[[2]] <- mat_dirichlet( length(unique(Cell_Type)),  sampleNumber, centre = cellNumber)
  BULK[[2]] <- round(BULK[[2]])
  BULK[[3]]<-matrix(0,ncol = sampleNumber,nrow = length(Cell_Type))
  #BULK[[4]]<-list()
  #names(BULK)<-c("Simu_Bulk","Cell_Type_Prop","Cell_Mat","Cell_Exp_Prop")
  names(BULK)<-c("Simu_Bulk","Cell_Type_Prop","Cell_Mat")
  rownames(BULK[[1]])<-rownames(sc_RNA)
  rownames(BULK[[2]])<-unique(Cell_Type)
  rownames(BULK[[3]])<-names(Cell_Type)
  
  for (i in 1:sampleNumber) {
    print(i)
    Cell_use<-NULL
    for (j in 1:nrow(BULK[[2]])) {
      #print(j)
      Cell_Cell<-sample(names(Cell_Type)[which(Cell_Type==rownames(BULK[[2]])[j])],BULK[[2]][j,i],replace = TRUE)
      Cell_use<-c(Cell_use,Cell_Cell)
      #TEMP<-sc_RNA[,Cell_Cell,drop=F]  #?? replicate
      #BULK[[4]][[j]][,i]<-rowSums(TEMP)
    }
    TEMP<-sc_RNA[,Cell_use]
    BULK[[1]][,i]<-rowSums(TEMP)/length(Cell_use)
    BULK[[3]][,i]<-table(Cell_use)[rownames(BULK[[3]])]
  }

  #BULK[[2]]<-sweep(BULK[[2]],2,colSums(BULK[[2]]),FUN = "/") #normalization!
  BULK[[3]][is.na(BULK[[3]])]<-0

  print("Done")
  print(Sys.time())
  return(BULK)
  
}

Bulk_Simu_coinfil_v1 <- function(sc_RNA, Cell_Type, cellNumber=1000, sampleNumber=100, row1, row2, cor)
{
  print(Sys.time())
  print('Start simulation without co-infiltration!')
  BULK <- list()
  BULK[[1]] <- matrix(0, nrow = nrow(sc_RNA), ncol = sampleNumber)
  print(cor)
  nnn <- mat_dirichlet_coinfil_colsum(length(unique(Cell_Type)), Sample,centre=cellnumber,GSEname_short,low=cor-0.05,up=cor+0.05, row1, row2,cor) 
  BULK[[2]] <- nnn
  BULK[[3]]<-matrix(0,ncol = sampleNumber,nrow = length(Cell_Type))
  BULK[[4]]<-list()
  names(BULK)<-c("Simu_Bulk","Cell_Type_Prop","Cell_Mat","Cell_Exp_Prop")
  rownames(BULK[[1]])<-rownames(sc_RNA)
  rownames(BULK[[2]])<-unique(Cell_Type)
  rownames(BULK[[3]])<-names(Cell_Type)
  
  for (i in 1:sampleNumber) {
    print(i)
    #BULK[[2]][,i]<-round(runif(nrow(BULK[[2]]),500,1000))  #already calculated dirichlet distribution 
    #BULK[[2]][,i]<-round(runif(nrow(BULK[[2]]),2000,4000))
    Cell_use<-NULL
    for (j in 1:nrow(BULK[[2]])) {
      #print(j)
      Cell_Cell<-sample(names(Cell_Type)[which(Cell_Type==rownames(BULK[[2]])[j])],BULK[[2]][j,i],replace = TRUE)
      Cell_use<-c(Cell_use,Cell_Cell)
      #TEMP<-sc_RNA[,Cell_Cell,drop=F]  #?? replicate
      #BULK[[4]][[j]][,i]<-rowSums(TEMP)
    }
    TEMP<-sc_RNA[,Cell_use]
    BULK[[1]][,i]<-rowSums(TEMP)/length(Cell_use)
    BULK[[3]][,i]<-table(Cell_use)[rownames(BULK[[3]])]
  }

  #BULK[[2]]<-sweep(BULK[[2]],2,colSums(BULK[[2]]),FUN = "/")
  BULK[[3]][is.na(BULK[[3]])]<-0

  print("Done")
  print(Sys.time())
  return(BULK)

}



mat_dirichlet <- function(celltype_number, sample_number,centre){  
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

replace_coinfil_two_vector <- function(ttt,loca_b=3,loca_t=9){
        aa <- ttt[loca_b,] #Bcell
        bb <- ttt[loca_t,] #Tcell

        aa1 <- aa[1:30];  aa2 <- aa[31:60] ;  aa3 <- aa[61:100];  
        bb1 <- bb[1:30];  bb2 <- bb[31:60];   bb3 <- bb[61:100];
        mean_aa1 <- mean(aa1); mean_aa2 <- mean(aa2); mean_aa3 <- mean(aa3); mean_bb1<-mean(bb1); mean_bb2<-mean(bb2);mean_bb3<-mean(bb3);
        scale <- 0.05
        sd_aa1<-scale*mean_aa1; sd_aa2<-scale*mean_aa2; sd_aa3<-scale*mean_aa3;
        sd_bb1<-scale*mean_bb1; sd_bb2<-scale*mean_bb2; sd_bb3<-scale*mean_bb3;

        twoV1 <- simcor(30,mean_aa1,sd_aa1,mean_bb1,sd_bb1, 0.4)
        twoV2 <- simcor(30,mean_aa2,sd_aa2,mean_bb2,sd_bb2, 0.4)
        twoV3 <- simcor(40,mean_aa3,sd_aa3,mean_bb3,sd_bb3, 0.4)
        xx <- c(twoV1[,1],twoV2[,1],twoV3[,1])
        yy <- c(twoV1[,2],twoV2[,2],twoV3[,2])

        ttt[loca_b,] <- xx
        ttt[loca_t,] <- yy

        #scale will change the correlation
        #sum_count <- apply(ttt,2,sum)
        #scale_factor <- rep(centre,100) / sum_count 
        #for(i in 1:100){
        #    ttt[,i] <- ttt[,i] * scale_factor[i]
        #}
        return(ttt) 
}

split_two_vector_predefined_cor <- function(ttt, loca_a=row1, loca_b=row2,ave_my,vari){
        #ttt_original <- ttt
        sum_ab <- ttt[loca_a,] + ttt[loca_b,]
        factor <- rnorm(100, mean=ave_my, sd=vari)
           factor[which(factor<0)] = 0.5
        aa <- factor * sum_ab
        bb <- sum_ab - aa
        co_infiltrate <- cor(aa,bb)
        print(paste("make co-infiltration: cor= ", co_infiltrate, sep=""))

        #ttt[loca_a,] <- aa
        #ttt[loca_b,] <- bb

        return(list(co_infiltrate=co_infiltrate,aa=aa,bb=bb) )
}


replace_two_vector_predefined_cor <- function(ttt, loca_a, loca_b,ave_my,vari){
        #sum_ab <- ttt[loca_a,] + ttt[loca_b,]
        scale=0.05
        dframe  <- simcor(100, mean( ttt[loca_a,]), scale*mean( ttt[loca_a,]), mean(ttt[loca_b,]), scale*mean(ttt[loca_b,]),ave_my )

        ll <- list()
        ll[[1]] <- cor(dframe[,1],dframe[,2])
        print(paste("cor=",ll[[1]],sep="") )
        ll[[2]] <- dframe[,1]
        ll[[3]] <- dframe[,2]

        return(ll)

}

simcor <- function (n, xmean, xsd, ymean, ysd, correlation) {
    x <- rnorm(n)
    y <- rnorm(n)
    z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) * 
             scale(resid(lm(y ~ x)))[,1]
    xresult <- xmean + xsd * scale(x)[,1]
    yresult <- ymean + ysd * z
    data.frame(x=xresult,y=yresult)
}

mat_dirichlet_coinfil <- function(celltype_number, sample_number,centre,GSEname){  
    ttt <- mat_dirichlet(celltype_number, sample_number,centre)
    vvlibrary <- unlist(strsplit(GSEname, "_", fixed=T))
    if("GSE103322" %in% vvlibrary){
        print("add co-infiltration to GSE103322")
       ttt1 <- replace_coinfil_two_vector(ttt,loca_b=3,loca_t=9)
       #cor(ttt[3,1:30],ttt[9,1:30]); cor(ttt[3,31:60],ttt[9,31:60]);cor(ttt[3,61:100],ttt[9,61:100]);
       #cor(ttt1[3,1:30],ttt1[9,1:30]); cor(ttt1[3,31:60],ttt1[9,31:60]);cor(ttt1[3,61:100],ttt1[9,61:100]);
       
    }
    return(ttt1)
}

mat_dirichlet_coinfil_colsum <- function(celltype_number, sample_number,centre,GSEname,low,up,row1,row2,cor){  
    print(paste("lower:",low,sep=""))
    print(paste("upper:",up,sep=""))

    ttt <- mat_dirichlet(celltype_number, sample_number,centre)
    vvlibrary <- unlist(strsplit(GSEname, "_", fixed=T))
    if(cor < 0.4) vari=0.12
    if(cor >=0.4 & cor < 0.8) vari=0.1
    if(cor >= 0.8 & cor < 0.9) vari=0.04
    if(cor >= 0.9) vari=0.02

    

    if("GSE103322" %in% vvlibrary){
        print("add co-infiltration to GSE103322")
        co_infiltrate <- 0
        lower = low
        upper = up
        while(T){
            ttt1 <- split_two_vector_predefined_cor(ttt,loca_a=row1,loca_b=row2, ave_my=cor, vari)
            #ttt1 <- replace_two_vector_predefined_cor(ttt, loca_a=row1, loca_b=row2,ave_my,vari)
            # 3,9 bcell,tcell
            #4,6 CAF,endothelial
            # 2,4 myofibro, CAF
            co_infiltrate <- ttt1[[1]]
            if(co_infiltrate>lower & co_infiltrate<upper)
            {
                 ttt[row1,] <- ttt1[[2]]#aa
                 ttt[row2,] <- ttt1[[3]]#bb
                 break
            }
        }
    }
    if("GSE81861" %in% vvlibrary){
        print("add co-infiltration to GSE81861")
        co_infiltrate <- 0
        lower = low
        upper = up
        while(T){
            ttt1 <- split_two_vector_predefined_cor(ttt,loca_a=row1,loca_b=row2, ave_my=cor, vari)
            #ttt1 <- replace_two_vector_predefined_cor(ttt, loca_a=row1, loca_b=row2,ave_my,vari)
            co_infiltrate <- ttt1[[1]]
            if(co_infiltrate>lower & co_infiltrate<upper)
            {
                 ttt[row1,] <- ttt1[[2]]#aa
                 ttt[row2,] <- ttt1[[3]]#bb
                 break
            }
        }
    }
      if("GSE75688" %in% vvlibrary){
        print("add co-infiltration to GSE75688")
        co_infiltrate <- 0
        lower = low
        upper = up
        while(T){
            ttt1 <- split_two_vector_predefined_cor(ttt,loca_a=row1,loca_b=row2, ave_my=cor, vari)
            #ttt1 <- replace_two_vector_predefined_cor(ttt, loca_a=row1, loca_b=row2,ave_my,vari)
            co_infiltrate <- ttt1[[1]]
            if(co_infiltrate>lower & co_infiltrate<upper)
            {
                 ttt[row1,] <- ttt1[[2]]#aa
                 ttt[row2,] <- ttt1[[3]]#bb
                 break
            }
        }
    }
    if("GSE72056" %in% vvlibrary){
        print("add co-infiltration to GSE72056")
        co_infiltrate <- 0
        lower = low
        upper = up
        while(T){
            ttt1 <- split_two_vector_predefined_cor(ttt,loca_a=row1,loca_b=row2, ave_my=cor, vari)
            #ttt1 <- replace_two_vector_predefined_cor(ttt, loca_a=row1, loca_b=row2,ave_my,vari)
            co_infiltrate <- ttt1[[1]]
            if(co_infiltrate>lower & co_infiltrate<upper)
            {
                 ttt[row1,] <- ttt1[[2]]#aa
                 ttt[row2,] <- ttt1[[3]]#bb
                 break
            }
        }
    }
  
    return(ttt)

}



Bulk_Simu_coinfil<-function(sc_RNA,Cell_Type,cellnumber=1000,GSEname_short,row1,row2,cor){
  Sample<-100
  print(Sys.time())
  print("Start Simulation")
  BULK<-list()
  BULK[[1]]<-matrix(0,ncol = Sample,nrow = nrow(sc_RNA))
  #BULK[[2]]<-matrix(0,ncol = Sample,nrow = length(unique(Cell_Type)))
  #prepare dirichlet distribute matrix for BULK[[2]]
  #BULK[[2]] <- mat_dirichlet( length(unique(Cell_Type)),  Sample， centre = 10000)
  #BULK[[2]] <- mat_dirichlet( length(unique(Cell_Type)),  Sample， centre = cellnumber)
  #BULK[[2]]<- mat_dirichlet_coinfil(length(unique(Cell_Type)), Sample,centre = cellnumber,GSEname_short)
  print(cor)
  nnn <- mat_dirichlet_coinfil_colsum(length(unique(Cell_Type)), Sample,centre=cellnumber,GSEname_short,low=cor-0.05,up=cor+0.05, row1, row2,cor) 
  BULK[[2]] <- nnn
#BULK[[2]] <- mat_dirichlet_coinfil(length(unique(Cell_Type)), Sample,centre = cellnumber,"GSE103322")
  BULK[[2]] <- round(BULK[[2]])
  BULK[[3]]<-matrix(0,ncol = Sample,nrow = length(Cell_Type))
  BULK[[4]]<-list()
  names(BULK)<-c("Simu_Bulk","Cell_Type_Prop","Cell_Mat","Cell_Exp_Prop")
  rownames(BULK[[1]])<-rownames(sc_RNA)
  rownames(BULK[[2]])<-unique(Cell_Type)
  rownames(BULK[[3]])<-names(Cell_Type)
  for (i in 1:length(unique(Cell_Type))) {
    BULK[[4]][[i]]<-BULK[[1]]*0
  }
  names(BULK[[4]])<-unique(Cell_Type)
  

  for (i in 1:Sample) {
    print(i)
    #BULK[[2]][,i]<-round(runif(nrow(BULK[[2]]),500,1000))  #already calculated dirichlet distribution 
    #BULK[[2]][,i]<-round(runif(nrow(BULK[[2]]),2000,4000))
    Cell_use<-NULL
    for (j in 1:nrow(BULK[[2]])) {
      #print(j)
      Cell_Cell<-sample(names(Cell_Type)[which(Cell_Type==rownames(BULK[[2]])[j])],BULK[[2]][j,i],replace = TRUE)
      Cell_use<-c(Cell_use,Cell_Cell)
      TEMP<-sc_RNA[,Cell_Cell,drop=F]  #?? replicate
      BULK[[4]][[j]][,i]<-rowSums(TEMP)
    }
    TEMP<-sc_RNA[,Cell_use]
    BULK[[1]][,i]<-rowSums(TEMP)/length(Cell_use)
    BULK[[3]][,i]<-table(Cell_use)[rownames(BULK[[3]])]
  }

  #BULK[[2]]<-sweep(BULK[[2]],2,colSums(BULK[[2]]),FUN = "/")
  BULK[[3]][is.na(BULK[[3]])]<-0
  TEMP<-sweep(BULK[[1]],2,colSums(BULK[[3]]),FUN = "*")
  for (j in 1:length(BULK[[4]])) {
    BULK[[4]][[j]]<-BULK[[4]][[j]]/TEMP
  }
  print("Done")
  print(Sys.time())
  return(BULK)
}


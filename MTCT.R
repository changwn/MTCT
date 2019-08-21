setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_06_11_simulator')

source('scsimu_sep.R')
source('extract_list.R')
source('ext_high_cor_gene.R')
source('check_pc.R')
source('sep_string.R')
source('pop_no_cor_pc.R')

load('C:/Users/wnchang/Documents/F/PhD_Research/2019_03_20_simulation_pipeline/E_MTAB_6149_selected_cells.RData')
label_name <- as.character(selected_label[,4])
names(label_name) <- colnames(scData)	
table(label_name)

#------------------------------------
#
# 1. using original label 
#
#------------------------------------
pseudo_BULK <- Bulk_Simu_no_coinfil_v1(scData, label_name, cellNumber=10000, sampleNumber=50)
data.matrix <- pseudo_BULK[[1]]
tProp <- pseudo_BULK[[2]]

pseudo_BULK_orig <- pseudo_BULK
high_gene_orig <- ext_high_cor_gene(pseudo_BULK_orig)

#---seurat
marker_gene <- high_gene_orig[[1]]
MAT <- scData
library(Seurat)
MAT_seurat_orig<-CreateSeuratObject(counts = MAT, project = "MAT", min.cells = 3, min.features = 200)
MAT_seurat_orig[["percent.mt"]] <- PercentageFeatureSet(MAT_seurat_orig, pattern = "^mt-")
#VlnPlot(MAT_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MAT_seurat_orig <- NormalizeData(MAT_seurat_orig)
MAT_seurat_orig<-FindVariableFeatures(MAT_seurat_orig, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MAT_seurat_orig)
MAT_seurat_orig <- ScaleData(MAT_seurat_orig, features = all.genes)
#MAT_seurat <- RunPCA(MAT_seurat, features = VariableFeatures(object = MAT_seurat))
MAT_seurat_orig <- RunPCA(MAT_seurat_orig, features = marker_gene)
MAT_seurat_orig <- JackStraw(MAT_seurat_orig, num.replicate = 100)
MAT_seurat_orig <- ScoreJackStraw(MAT_seurat_orig, dims = 1:20)
#JackStrawPlot(MAT_seurat, dims = 1:20)
#ElbowPlot(MAT_seurat)
MAT_seurat_orig <- FindNeighbors(MAT_seurat_orig, dims = 1:20)
MAT_seurat_orig <- FindClusters(MAT_seurat_orig, resolution = 0.5)
MAT_seurat_orig <- RunTSNE(MAT_seurat_orig, dims = 1:20)
DimPlot(MAT_seurat_orig, reduction = "tsne", label=T)

MAT_ident_orig<-as.vector(MAT_seurat_orig@meta.data[,6])
names(MAT_ident_orig) <- rownames(MAT_seurat_orig@meta.data)

library(cluster)
tsne_position <- MAT_seurat_orig@reductions$tsne@cell.embeddings
test_dist<-dist(tsne_position,method = "canberra")
cluster_label <- strtoi(MAT_ident_orig, base = 0L)
sil<-silhouette(cluster_label,test_dist)
sum(sil[,3])

# library(cluster)
# tsne_position <- MAT_seurat_orig@reductions$tsne@cell.embeddings
# test_dist<-dist(tsne_position,method = "canberra")
# #cluster_label <- strtoi(label_name, base = 0L)
# label_orig <- as.factor(label_name)
# levels(label_orig) <- 1:length(levels(label_orig))
# label_orig <- as.numeric(label_orig)
# sil<-silhouette(label_orig,test_dist)
# sum(sil[,3])
#---result

#-------------------------------
#
# 2. seurat clustering, plot tSNE
#
#-------------------------------

#---seurat
MAT <- scData
library(Seurat)
MAT_seurat<-CreateSeuratObject(counts = MAT, project = "MAT", min.cells = 3, min.features = 200)
MAT_seurat[["percent.mt"]] <- PercentageFeatureSet(MAT_seurat, pattern = "^mt-")
#VlnPlot(MAT_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MAT_seurat <- NormalizeData(MAT_seurat)
MAT_seurat<-FindVariableFeatures(MAT_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MAT_seurat)
MAT_seurat <- ScaleData(MAT_seurat, features = all.genes)
MAT_seurat <- RunPCA(MAT_seurat, features = VariableFeatures(object = MAT_seurat))
#MAT_seurat <- RunPCA(MAT_seurat, features = marker_gene) #use default seurat, no marker gene
MAT_seurat <- JackStraw(MAT_seurat, num.replicate = 100)
MAT_seurat <- ScoreJackStraw(MAT_seurat, dims = 1:20)
#JackStrawPlot(MAT_seurat, dims = 1:20)
#ElbowPlot(MAT_seurat)
MAT_seurat <- FindNeighbors(MAT_seurat, dims = 1:20)
MAT_seurat <- FindClusters(MAT_seurat, resolution = 0.5)
MAT_seurat <- RunTSNE(MAT_seurat, dims = 1:20)
DimPlot(MAT_seurat, reduction = "tsne", label=T)
MAT_ident<-as.vector(MAT_seurat@meta.data[,6])
names(MAT_ident) <- rownames(MAT_seurat@meta.data)

library(cluster)
tsne_position <- MAT_seurat@reductions$tsne@cell.embeddings
test_dist<-dist(tsne_position,method = "canberra")
cluster_label <- strtoi(MAT_ident, base = 0L)
sil<-silhouette(cluster_label,test_dist)
sum(sil[,3])

#----------------------
#
# 3. MTCT: simulate using seurat cluster label, run ictd, clustering using ctes3 marker
#
#----------------------

pseudo_BULK <- Bulk_Simu_no_coinfil_v1(scData, MAT_ident, cellNumber=10000, sampleNumber=50)
data.matrix <- pseudo_BULK[[1]]
tProp <- pseudo_BULK[[2]]

#-----------------
mainDir <- "C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline"
setwd(mainDir)

#source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Loading_header.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Loading_header_new_IM.r")
#source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Loading_header_new_Brain.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/10_18_update/All_new_functions.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/10_18_update/Step_3_linking_graph_function.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_10_08_extract_data/R2_all_function.R")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/temp_func.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/temp_func3.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/temp_func3_new.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/01-23NMF/NMF_functions_new.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Step2plus_Celltype_marker_Hierarchical_CTES.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Step2plus_Celltype_marker_Hierarchical_CTES2.r")

library(NMF)
set.seed(123456)
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/nmf.library.R")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/ini.R")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/Step_3_linking_graph_function.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/R4RR_selection_functions.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/NMF_functions_new.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/NMF_method1_test_version.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/NMF_method1_test_version2.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/NMF_method1_test_version3.r")


source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/debug_20190217.R")  #debug "compute_CompRowspace_NN_uni" function
#-----------------
###########################################################
if(length(colnames(data.matrix)) == 0) {
  warning("input data do NOT have colnames")
  colnames(data.matrix) <- paste( "Setsample", 1:ncol(data.matrix), sep="")    
}
data.matrix <- rm_zero_row(data.matrix)

###########  RNAseq log(X+1)
d.matrix<-log(data.matrix + 1)
d.matrix<-as.matrix(d.matrix)
#d.matrix <- data.matrix #log will cause bad performance sometimes!!!!!!!!!72056,81861
########### normalize data2, same for RNA-seq data
data01<-normalize_data2(d.matrix)
data0<-d.matrix

########### prepare input data for cancer module identification
data2<-data0

###########take the highest expressed probe for the genes with duplicated probes
data21<-data2[order(-apply(data2,1,mean)),]
data22<-data21[unique(rownames(data21)),]
data22_code<-data22[intersect(rownames(data22),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]
data23<-data22[intersect(rownames(data22),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]
data23<-normalize_data2(data23)		#BUG source: get NaN value after "normalize_data2" because some zero rows
##############

#data_ccc<-data_CORS_cancer
data_CORS_cancer <- data23
data_ccc <- data23
list_c1<-MRHCA_IM_compute_MR(data_CORS_cancer=data_ccc,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)	#MR_1_rank_test_function_v1.1
MR_IM_result_new_c<-MRHCA_IM_compute_full_pub_new(data_CORS_cancer=data_ccc,list_c=list_c1,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
tg_key="nonono"
#list_new_c2<-Process_MR_IM_result_new(MR_IM_result_new_c,tg_key_c=tg_key,cell_type_enrich_cut=0.4,cor_cut0=0.8,num_cut=10,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
#R1_filter_step1_results_new<-R1_list_filtering_step1_new(list_new_c2,data_CORS_cancer=data_ccc,max_cut=20,cutn0=10,cut10=0.8,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)

list_new_c2<-Process_MR_IM_result_new(MR_IM_result_new_c,tg_key_c=tg_key,cell_type_enrich_cut=0.4,cor_cut0=0.7,num_cut=7,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
R1_filter_step1_results_new<-R1_list_filtering_step1_new(list_new_c2,data_CORS_cancer=data_ccc,max_cut=7,cutn0=7,cut10=0.7,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)

tg_R1_lists<-R1_filter_step1_results_new[[4]]
print(length(tg_R1_lists))
#R1_selectedCM_step2_results_new<-Step2plus_Celltype_marker_inference(tg_R1_lists,data_CORS_cancer,tg_R1_cut=6,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],resolution_level=resolution_level0)
#R1_selectedCM_step2_results_new <- Step2plus_Celltype_marker_inference_new(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=20)
#comb <- c(R1_selectedCM_step2_results_new[[1]],R1_selectedCM_step2_results_new[[2]])

# R1_selectedCM_step2_results_Hclust<-Step2plus_Celltype_marker_inference_Hclust(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=40)
# R1_selectedCM_step2_results_fixedCT<-Step2plus_Celltype_marker_inference_fixedCT(tg_R1_lists,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]])
# R1_selectedCM_step2_results_HP<-Step2plus_Celltype_marker_inference_Hclust_plus(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.5,resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.85,IM_reso_level=IM_reso_level)
# R1_selectedCM_step2_results_HCTES<-Step2plus_Celltype_marker_inference_Hclust_CTES(tg_R1_lists,data_CORS_cancer,tg_R1_cut=5,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.95,IM_reso_level=IM_reso_level,extra_ctn=3)
# CTES<-c(R1_selectedCM_step2_results_HCTES[[1]],R1_selectedCM_step2_results_HCTES[[2]])

#R1_selectedCM_step2_results_Hierarchical_CTES<-Step2plus_Celltype_marker_Hierarchical_CTES(tg_R1_lists=tg_R1_lists,data_CORS_cancer=data_CORS_cancer,data.matrix=data.matrix,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.45,tg_R1_cut=6)
R1_selectedCM_step2_results_Hierarchical_CTES<- Step2plus_Celltype_marker_Hierarchical_CTES3(tg_R1_lists,data_CORS_cancer,data.matrix,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.45,tg_R1_cut=6,ext_cn=10,CT_balance=1)  #72056_ext=6, ctBALANCE=1
CTES3<-R1_selectedCM_step2_results_Hierarchical_CTES[[1]]
c9_Hclust_CTES<-Compute_Rbase_SVD(data.matrix,CTES3)
d9<-cor(t(c9_Hclust_CTES),t(tProp))
apply(d9,2,max)


#--------------------------

MAT <- scData
ictd_genelist <- extract_list(CTES3)
#ictd_genelist <- extract_list(tg_R1_lists)
library(Seurat)
MAT_seurat_mtct<-CreateSeuratObject(counts = MAT, project = "MAT", min.cells = 3, min.features = 200)
MAT_seurat_mtct[["percent.mt"]] <- PercentageFeatureSet(MAT_seurat_mtct, pattern = "^mt-")
#VlnPlot(MAT_seurat_mtct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MAT_seurat_mtct <- NormalizeData(MAT_seurat_mtct)
MAT_seurat_mtct<-FindVariableFeatures(MAT_seurat_mtct, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MAT_seurat_mtct)
MAT_seurat_mtct <- ScaleData(MAT_seurat_mtct, features = all.genes)
#MAT_seurat <- RunPCA(MAT_seurat, features = VariableFeatures(object = MAT_seurat))
MAT_seurat_mtct <- RunPCA(MAT_seurat_mtct, features = ictd_genelist)
MAT_seurat_mtct <- JackStraw(MAT_seurat_mtct, num.replicate = 100)
MAT_seurat_mtct <- ScoreJackStraw(MAT_seurat_mtct, dims = 1:20)
#JackStrawPlot(MAT_seurat_mtct, dims = 1:20)
#ElbowPlot(MAT_seurat_mtct)
MAT_seurat_mtct <- FindNeighbors(MAT_seurat_mtct, dims = 1:20)
MAT_seurat_mtct <- FindClusters(MAT_seurat_mtct, resolution = 0.5)
MAT_seurat_mtct <- RunTSNE(MAT_seurat_mtct, dims = 1:20)
DimPlot(MAT_seurat_mtct, reduction = "tsne", label=T)

MAT_ident_mtct<-as.vector(MAT_seurat_mtct@meta.data[,6])
names(MAT_ident_mtct) <- rownames(MAT_seurat_mtct@meta.data)

# #calculate silhouette
# sc_data_ori <- t(scData)
# library(clusterCrit)
# label_CTES3 <- strtoi(MAT_ident_CTES3, base = 0L)
# intCriteria(sc_data_ori,label_CTES3,"Silhouette")
# label_orig <- as.numeric(MAT_ident_orig)
# intCriteria(sc_data_ori,label_orig,"Silhouette")

library(cluster)
tsne_position <- MAT_seurat_mtct@reductions$tsne@cell.embeddings
test_dist<-dist(tsne_position,method = "canberra")
cluster_label <- strtoi(MAT_ident_mtct, base = 0L)
sil<-silhouette(cluster_label,test_dist)
sum(sil[,3])

heatmap(table(MAT_ident_mtct, MAT_ident))





# #-----------------
# # find seurat marker gene
# MAT.rna.markers <- FindAllMarkers(MAT_seurat, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
# # find mtct marker gene
# MAT.rna.markers_mtct <- FindAllMarkers(MAT_seurat_mtct, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)


# #1. plot mtct T marker on mtct/seurat plot
# MAT.rna.markers[which(MAT.rna.markers[, 6] == 3), ]
# rownames(MAT.rna.markers[which(MAT.rna.markers[, 6] == 3), ])[1:6]
# FeaturePlot(MAT_seurat_mtct, features = c("CD3D","LCK","SH2D1A","CD96","CD8A","TIGIT" ), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)
# # MAT.rna.markers_mtct[which(MAT.rna.markers_mtct[, 6] == 4), ][,7][1:6]
# # FeaturePlot(MAT_seurat_mtct, features = c("CD2","CD3D","CD3E","ITM2A","CORO1A","LCK"), 
# #             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)
# FeaturePlot(MAT_seurat_mtct, features = c("CD3D" ,"CD3E" ,"CD2" ,"IL32" ,"CORO1A" ,"LCK"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)

# FeaturePlot(MAT_seurat, features = c("CD3D","LCK","SH2D1A","CD96","CD8A","TIGIT" ), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)
# FeaturePlot(MAT_seurat, features = c("CD3D" ,"CD3E" ,"CD2" ,"IL32" ,"CORO1A" ,"LCK"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)

# FeaturePlot(MAT_seurat_orig, features = c("CD3D","LCK","SH2D1A","CD96","CD8A","TIGIT" ), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)
# FeaturePlot(MAT_seurat_orig, features = c("CD3D" ,"CD3E" ,"CD2" ,"IL32" ,"CORO1A" ,"LCK"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)


# #2. B marker
# MAT.rna.markers[which(MAT.rna.markers[, 6] == 11), ]
# MAT.rna.markers[which(MAT.rna.markers[, 6] == 11), ][,7][1:6]
# FeaturePlot(MAT_seurat_mtct, features = c("CD79A","FCRL5","TNFRSF17","PNOC","MZB1","CD19"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)
# FeaturePlot(MAT_seurat_mtct, features = c("CD79A","MZB1","ISG20","DERL3","PIM2","ITM2C"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)

# FeaturePlot(MAT_seurat, features = c("CD79A","FCRL5","TNFRSF17","PNOC","MZB1","CD19"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)
# FeaturePlot(MAT_seurat, features = c("CD79A","MZB1","ISG20","DERL3","PIM2","ITM2C"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)

# FeaturePlot(MAT_seurat_orig, features = c("CD79A","FCRL5","TNFRSF17","PNOC","MZB1","CD19"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)
# FeaturePlot(MAT_seurat_orig, features = c("CD79A","MZB1","ISG20","DERL3","PIM2","ITM2C"), 
#             min.cutoff = "q05", max.cutoff = "q95", ncol = 3)


# # automatically
# for(i in 1:4)
# {
#  	i <- i + 1
# 	file_name <- paste(i,'_',names(CTES3)[i],'_mtct.On.mtct--seu.On.mtct--mtct.on.seu--seu.On.seu','.pdf',sep='')
# 	pdf(file_name)
# 	#i = 2
# 	mtct_gene_tmp <- CTES3[[i]]
# 	cell_type_tmp <- names(CTES3)[i]
# 	#return cluster id
# 	print(cell_type_tmp)
# 	print(mtct_gene_tmp)
# 	print(MAT.rna.markers[match(intersect(mtct_gene_tmp,MAT.rna.markers[, 7]), MAT.rna.markers[, 7] ),])
# 	cluster_id_tmp <- names(table(as.numeric(MAT.rna.markers[match(intersect(mtct_gene_tmp, MAT.rna.markers[, 7]), MAT.rna.markers[, 7] ),][, 6])-1))[1]	
# 	cluster_id_tmp <-  strtoi(cluster_id_tmp, base = 0L)
# 	seurat_gene_tmp <- MAT.rna.markers[which(MAT.rna.markers[, 6] == cluster_id_tmp), ][,7][1:6]


# 	FeaturePlot(MAT_seurat_mtct, features = mtct_gene_tmp, 
# 	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
# 	FeaturePlot(MAT_seurat_mtct, features = seurat_gene_tmp, 
# 	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)

# 	FeaturePlot(MAT_seurat, features = mtct_gene_tmp, 
# 	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
# 	FeaturePlot(MAT_seurat, features = seurat_gene_tmp, 
# 	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
# 	dev.off()

# }


###MAT.rna.markers_mtct[which(MAT.rna.markers_mtct[, 6] == 4), ][,7][1:6]

# # auto: directly use mtct cluster
# label_relation <- table(MAT_ident, MAT_ident_mtct)
# cluster_mtct <- colnames(label_relation)
# for(i in 1:length(cluster_mtct))
# {
# 	i <- i+1
# 	cluster_mtct_tmp <- cluster_mtct[i]
# 	cluster_seurat_tmp <- rownames(label_relation)[which(label_relation[,i]==max(label_relation[,i]))]

# 	mtct_gene_tmp <- MAT.rna.markers_mtct[which(MAT.rna.markers_mtct[, 6] == cluster_mtct_tmp), ][,7][1:6]
# 	seurat_gene_tmp <- MAT.rna.markers[which(MAT.rna.markers[, 6] == cluster_seurat_tmp), ][,7][1:6]

# 	file_name <- paste('cluster_', i, '_mtct.On.mtct--seu.On.mtct--mtct.on.seu--seu.On.seu','.pdf',sep='')
# 	pdf(file_name)

# 	FeaturePlot(MAT_seurat_mtct, features = mtct_gene_tmp, 
# 	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
# 	FeaturePlot(MAT_seurat_mtct, features = seurat_gene_tmp, 
# 	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)

# 	FeaturePlot(MAT_seurat, features = mtct_gene_tmp, 
# 	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
# 	FeaturePlot(MAT_seurat, features = seurat_gene_tmp, 
# 	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
# 	dev.off()

# }

# use cluster label to filter CTES3: doesn't work
# use pc2 component to filter!
pass_id <- c()
for(i in 1:length(CTES3))
{
	#i <- i+1
	mtct_gene_tmp <- CTES3[[i]]
	cell_type_tmp <- names(CTES3)[i]
	print(cell_type_tmp)
	print(mtct_gene_tmp)
	print(MAT.rna.markers_mtct[match(intersect(mtct_gene_tmp,MAT.rna.markers_mtct[, 7]), MAT.rna.markers_mtct[, 7] ),])
	ttt <- MAT.rna.markers_mtct[match(intersect(mtct_gene_tmp,MAT.rna.markers_mtct[, 7]), MAT.rna.markers_mtct[, 7] ),]
	if(nrow(ttt) == 0) next
	if(check_pc(ttt) && nrow(ttt)>4)
	{
		print(paste('CTES...', i, '_pass independent check!',sep=''))
		#FeaturePlot(MAT_seurat_mtct, features = mtct_gene_tmp, 
	    #        min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
		pass_id <- c(pass_id, i)
	}

			 # FeaturePlot(MAT_seurat_mtct, features = mtct_gene_tmp, 
	   #           min.cutoff = "q05", max.cutoff = "q95", ncol = 2)

}
CTES_pass <- CTES3[pass_id]
# automatically
for(i in 1:length(CTES_pass))
{
 	i <- i + 1
	file_name <- paste(i,'_',names(CTES_pass)[i],'_newCTES3_tsne','.pdf',sep='')
	pdf(file_name)
	#i = 2
	mtct_gene_tmp <- CTES_pass[[i]]
	cell_type_tmp <- names(CTES_pass)[i]
	#return cluster id
	print(cell_type_tmp)
	print(mtct_gene_tmp)
	print(MAT.rna.markers[match(intersect(mtct_gene_tmp,MAT.rna.markers[, 7]), MAT.rna.markers[, 7] ),])
	cluster_id_tmp <- names(table(as.numeric(MAT.rna.markers[match(intersect(mtct_gene_tmp, MAT.rna.markers[, 7]), MAT.rna.markers[, 7] ),][, 6])-1))[1]	
	cluster_id_tmp <-  strtoi(cluster_id_tmp, base = 0L)
	seurat_gene_tmp <- MAT.rna.markers[which(MAT.rna.markers[, 6] == cluster_id_tmp), ][,7][1:6]


	FeaturePlot(MAT_seurat_mtct, features = mtct_gene_tmp, 
	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
	FeaturePlot(MAT_seurat_mtct, features = seurat_gene_tmp, 
	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)

	FeaturePlot(MAT_seurat, features = mtct_gene_tmp, 
	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
	FeaturePlot(MAT_seurat, features = seurat_gene_tmp, 
	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
	dev.off()

}


#------2. check with published marker library  XXXXXXXXXXXXXXXXX
human_marker <- read.table('Single_cell_markers.txt', header=T, sep='\t')
#rownames(human_marker) <- human_marker[,6]
length(table(human_marker[,6]))
human_marker_B <- human_marker[grep('B cell', human_marker[, 6]), ]
B_cell_lib <- c()
human_marker_B[, 8] <- as.character(human_marker_B[, 8])
for(i in 1:nrow(human_marker_B))
{
	B_cell_lib <- sort(unique(c(B_cell_lib, sep_string(human_marker_B[i, 8], 10))))
	#print(human_marker_B[i, 8])
}
B_cell_lib <- toupper(B_cell_lib)
intersect(B_cell_lib, CTES_pass[[6]])
## few overlap!


#----3. draw heatmap
CTES_pass_gene <- extract_list(CTES_pass)
library(gplots)
library("RColorBrewer")
#pheatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(299)
pheatmap_color <- colorRampPalette(c("blue", "white","red"))(100)
ht_dt <- scData[CTES_pass_gene, ]
# pdf('sc_heatmap_1.pdf')
# heatmap.2(ht_dt, col=pheatmap_color,scale="row", trace='none',Rowv = F, Colv = T)
# dev.off()

DoHeatmap(MAT_seurat_mtct,features = CTES_pass_gene)

# add other cell for CTES_pass
# heatmap for all CTES3
#DoHeatmap(MAT_seurat_mtct,features = ictd_genelist)


#------------------------
#
#make up: add more bases 
#
#------------------------
library(factoextra)
# example for PCA
# data(decathlon2)
# decathlon2.active <- decathlon2[1:23, 1:10]
# head(decathlon2.active[, 1:6])
# res.pca <- prcomp(decathlon2.active, scale = TRUE)
# fviz_eig(res.pca)
# #Graph of individuals. Individuals with a similar profile are grouped together.
# fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
#              )
# #Graph of variables. Positive correlated variables point to the same side of the plot.
# # Negative correlated variables point to opposite sides of the graph.
# fviz_pca_var(res.pca,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
#              )
# #Biplot of individuals and variables
# fviz_pca_biplot(res.pca, repel = TRUE,
#                 col.var = "#2E9FDF", # Variables color
#                 col.ind = "#696969"  # Individuals color
#                 )

CTES_row_bases <- Compute_Rbase_SVD(scData, CTES3)
colnames(CTES_row_bases) <- colnames(scData)
rownames(CTES_row_bases) <- paste(c(1:nrow(CTES_row_bases)), '_', rownames(CTES_row_bases), sep='')
# length(ictd_genelist)
# length(unique(ictd_genelist))
res.pca <- prcomp(t(CTES_row_bases), scale = TRUE)
plot.new()
fviz_eig(res.pca)
# plot dot of the each player (each single cell)
# plot.new()
# fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
#              )
# plot arrow of the cell type
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F     # Avoid text overlapping
             )
#Biplot of individuals and variables
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )


res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 




cor(CTES_row_bases[2,], CTES_row_bases[8,])

#add_id_candidate <- pop_no_cor_pc(CTES_row_bases, pass_id)
pass_and_add <- pop_no_cor_pc_dynamic(CTES_row_bases, pass_id)


pass_add_gene <- extract_list(CTES3[pass_and_add])
names(CTES3)[pass_and_add]

DoHeatmap(MAT_seurat_mtct,features = pass_add_gene)

for(i in 1:length(pass_and_add))
{
	i <- i + 1
	vl_file <- paste(pass_and_add[i], '_', names(CTES3)[pass_and_add][i], '_vlplot.pdf', sep='')
	pdf(vl_file)
	VlnPlot(object = MAT_seurat_mtct, features = CTES3[pass_and_add][[i]], pt.size = 0.1)
	dev.off()


}

RidgePlot(MAT_seurat_mtct, features = CTES3[pass_and_add][[1]], ncol = 3)



#---------------------
#
# CITE-seq validation
#
#---------------------









# setwd('C:/Users/wnchang/Documents/F/PhD_Research/2019_06_11_simulator')
# save(list=c('MAT_ident', 'scData','Bulk_Simu_no_coinfil_v1'), file="lung_sc_c1.RData")

# save(MAT_ident, file='cluster_origi.RData')
# 
# heatmap(table(MAT_ident_orig, MAT_ident_CTES3))
# 
# heatmap(table(MAT_ident_orig, MAT_ident_R1))


#--------------------------
#
#simulation 2nd time using MAT_ident_mtct
#
#--------------------------

pseudo_BULK <- Bulk_Simu_no_coinfil_v1(scData, MAT_ident_mtct, cellNumber=10000, sampleNumber=50)
data.matrix <- pseudo_BULK[[1]]
tProp <- pseudo_BULK[[2]]

#-----------------
mainDir <- "C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline"
setwd(mainDir)

#source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Loading_header.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Loading_header_new_IM.r")
#source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Loading_header_new_Brain.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/10_18_update/All_new_functions.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/10_18_update/Step_3_linking_graph_function.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_10_08_extract_data/R2_all_function.R")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/temp_func.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/temp_func3.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/temp_func3_new.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/01-23NMF/NMF_functions_new.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Step2plus_Celltype_marker_Hierarchical_CTES.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/Step2plus_Celltype_marker_Hierarchical_CTES2.r")

library(NMF)
set.seed(123456)
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/nmf.library.R")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/ini.R")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/Step_3_linking_graph_function.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/R4RR_selection_functions.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/NMF_functions_new.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/NMF_method1_test_version.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/NMF_method1_test_version2.r")
source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/02-01-NMF_function/NMF_method1_test_version3.r")


source("C:/Users/wnchang/Documents/F/PhD_Research/2018_09_11_ICAD_pipeline/debug_20190217.R")  #debug "compute_CompRowspace_NN_uni" function
#-----------------
###########################################################
if(length(colnames(data.matrix)) == 0) {
  warning("input data do NOT have colnames")
  colnames(data.matrix) <- paste( "Setsample", 1:ncol(data.matrix), sep="")    
}
data.matrix <- rm_zero_row(data.matrix)

###########  RNAseq log(X+1)
d.matrix<-log(data.matrix + 1)
d.matrix<-as.matrix(d.matrix)
#d.matrix <- data.matrix #log will cause bad performance sometimes!!!!!!!!!72056,81861
########### normalize data2, same for RNA-seq data
data01<-normalize_data2(d.matrix)
data0<-d.matrix

########### prepare input data for cancer module identification
data2<-data0

###########take the highest expressed probe for the genes with duplicated probes
data21<-data2[order(-apply(data2,1,mean)),]
data22<-data21[unique(rownames(data21)),]
data22_code<-data22[intersect(rownames(data22),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]
data23<-data22[intersect(rownames(data22),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]
data23<-normalize_data2(data23)		#BUG source: get NaN value after "normalize_data2" because some zero rows
##############

#data_ccc<-data_CORS_cancer
data_CORS_cancer <- data23
data_ccc <- data23
list_c1<-MRHCA_IM_compute_MR(data_CORS_cancer=data_ccc,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)	#MR_1_rank_test_function_v1.1
MR_IM_result_new_c<-MRHCA_IM_compute_full_pub_new(data_CORS_cancer=data_ccc,list_c=list_c1,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
tg_key="nonono"
#list_new_c2<-Process_MR_IM_result_new(MR_IM_result_new_c,tg_key_c=tg_key,cell_type_enrich_cut=0.4,cor_cut0=0.8,num_cut=10,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
#R1_filter_step1_results_new<-R1_list_filtering_step1_new(list_new_c2,data_CORS_cancer=data_ccc,max_cut=20,cutn0=10,cut10=0.8,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)

list_new_c2<-Process_MR_IM_result_new(MR_IM_result_new_c,tg_key_c=tg_key,cell_type_enrich_cut=0.4,cor_cut0=0.7,num_cut=7,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
R1_filter_step1_results_new<-R1_list_filtering_step1_new(list_new_c2,data_CORS_cancer=data_ccc,max_cut=7,cutn0=7,cut10=0.7,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)

tg_R1_lists<-R1_filter_step1_results_new[[4]]
print(length(tg_R1_lists))
#R1_selectedCM_step2_results_new<-Step2plus_Celltype_marker_inference(tg_R1_lists,data_CORS_cancer,tg_R1_cut=6,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],resolution_level=resolution_level0)
#R1_selectedCM_step2_results_new <- Step2plus_Celltype_marker_inference_new(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=20)
#comb <- c(R1_selectedCM_step2_results_new[[1]],R1_selectedCM_step2_results_new[[2]])

# R1_selectedCM_step2_results_Hclust<-Step2plus_Celltype_marker_inference_Hclust(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=40)
# R1_selectedCM_step2_results_fixedCT<-Step2plus_Celltype_marker_inference_fixedCT(tg_R1_lists,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]])
# R1_selectedCM_step2_results_HP<-Step2plus_Celltype_marker_inference_Hclust_plus(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.5,resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.85,IM_reso_level=IM_reso_level)
# R1_selectedCM_step2_results_HCTES<-Step2plus_Celltype_marker_inference_Hclust_CTES(tg_R1_lists,data_CORS_cancer,tg_R1_cut=5,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.95,IM_reso_level=IM_reso_level,extra_ctn=3)
# CTES<-c(R1_selectedCM_step2_results_HCTES[[1]],R1_selectedCM_step2_results_HCTES[[2]])

#R1_selectedCM_step2_results_Hierarchical_CTES<-Step2plus_Celltype_marker_Hierarchical_CTES(tg_R1_lists=tg_R1_lists,data_CORS_cancer=data_CORS_cancer,data.matrix=data.matrix,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.45,tg_R1_cut=6)
R1_selectedCM_step2_results_Hierarchical_CTES<- Step2plus_Celltype_marker_Hierarchical_CTES3(tg_R1_lists,data_CORS_cancer,data.matrix,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.45,tg_R1_cut=6,ext_cn=10,CT_balance=1)  #72056_ext=6, ctBALANCE=1
CTES3<-R1_selectedCM_step2_results_Hierarchical_CTES[[1]]
c9_Hclust_CTES<-Compute_Rbase_SVD(data.matrix,CTES3)
d9<-cor(t(c9_Hclust_CTES),t(tProp))
apply(d9,2,max)


#--------------

MAT <- scData
ictd_genelist2 <- extract_list(CTES3)
#ictd_genelist <- extract_list(tg_R1_lists)
library(Seurat)
MAT_seurat_mtct2<-CreateSeuratObject(counts = MAT, project = "MAT", min.cells = 3, min.features = 200)
MAT_seurat_mtct2[["percent.mt"]] <- PercentageFeatureSet(MAT_seurat_mtct2, pattern = "^mt-")
#VlnPlot(MAT_seurat_mtct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MAT_seurat_mtct2 <- NormalizeData(MAT_seurat_mtct2)
MAT_seurat_mtct2<-FindVariableFeatures(MAT_seurat_mtct2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MAT_seurat_mtct2)
MAT_seurat_mtct2 <- ScaleData(MAT_seurat_mtct2, features = all.genes)
#MAT_seurat <- RunPCA(MAT_seurat, features = VariableFeatures(object = MAT_seurat))
MAT_seurat_mtct2 <- RunPCA(MAT_seurat_mtct2, features = ictd_genelist2)
MAT_seurat_mtct2 <- JackStraw(MAT_seurat_mtct2, num.replicate = 100)
MAT_seurat_mtct2 <- ScoreJackStraw(MAT_seurat_mtct2, dims = 1:20)
#JackStrawPlot(MAT_seurat_mtct, dims = 1:20)
#ElbowPlot(MAT_seurat_mtct)
MAT_seurat_mtct2 <- FindNeighbors(MAT_seurat_mtct2, dims = 1:20)
MAT_seurat_mtct2 <- FindClusters(MAT_seurat_mtct2, resolution = 0.5)
MAT_seurat_mtct2 <- RunTSNE(MAT_seurat_mtct2, dims = 1:20)
DimPlot(MAT_seurat_mtct2, reduction = "tsne", label=T)

MAT_ident_mtct2<-as.vector(MAT_seurat_mtct2@meta.data[,6])
names(MAT_ident_mtct2) <- rownames(MAT_seurat_mtct2@meta.data)

# #calculate silhouette
# sc_data_ori <- t(scData)
# library(clusterCrit)
# label_CTES3 <- strtoi(MAT_ident_CTES3, base = 0L)
# intCriteria(sc_data_ori,label_CTES3,"Silhouette")
# label_orig <- as.numeric(MAT_ident_orig)
# intCriteria(sc_data_ori,label_orig,"Silhouette")

library(cluster)
tsne_position <- MAT_seurat_mtct2@reductions$tsne@cell.embeddings
test_dist<-dist(tsne_position,method = "canberra")
cluster_label <- strtoi(MAT_ident_mtct2, base = 0L)
sil<-silhouette(cluster_label,test_dist)
sum(sil[,3])

#---------------------------

# automatically
for(i in 1:4)
{
 	i <- i + 1
	file_name <- paste(i,'_',names(CTES3)[i],'_mtct.On.mtct--seu.On.mtct--mtct.on.seu--seu.On.seu','.pdf',sep='')
	pdf(file_name)
	#i = 2
	mtct_gene_tmp <- CTES3[[i]]
	cell_type_tmp <- names(CTES3)[i]
	#return cluster id
	print(cell_type_tmp)
	print(mtct_gene_tmp)
	print(MAT.rna.markers[match(intersect(mtct_gene_tmp,MAT.rna.markers[, 7]), MAT.rna.markers[, 7] ),])
	cluster_id_tmp <- names(table(as.numeric(MAT.rna.markers[match(intersect(mtct_gene_tmp, MAT.rna.markers[, 7]), MAT.rna.markers[, 7] ),][, 6])-1))[1]	
	cluster_id_tmp <-  strtoi(cluster_id_tmp, base = 0L)
	seurat_gene_tmp <- MAT.rna.markers[which(MAT.rna.markers[, 6] == cluster_id_tmp), ][,7][1:6]


	FeaturePlot(MAT_seurat_mtct, features = mtct_gene_tmp, 
	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
	FeaturePlot(MAT_seurat_mtct, features = seurat_gene_tmp, 
	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)

	FeaturePlot(MAT_seurat, features = mtct_gene_tmp, 
	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
	FeaturePlot(MAT_seurat, features = seurat_gene_tmp, 
	            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
	dev.off()

}




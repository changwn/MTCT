
MAT <- data.matrix

library(Seurat)

MAT_seurat<-CreateSeuratObject(counts = MAT, project = "MAT", min.cells = 3, min.features = 200)
MAT_seurat[["percent.mt"]] <- PercentageFeatureSet(MAT_seurat, pattern = "^mt-")

VlnPlot(MAT_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MAT_seurat <- NormalizeData(MAT_seurat)

MAT_seurat<-FindVariableFeatures(MAT_seurat, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MAT_seurat)
MAT_seurat <- ScaleData(MAT_seurat, features = all.genes)

MAT_seurat <- RunPCA(MAT_seurat, features = VariableFeatures(object = MAT_seurat))

MAT_seurat <- JackStraw(MAT_seurat, num.replicate = 100)
MAT_seurat <- ScoreJackStraw(MAT_seurat, dims = 1:20)

JackStrawPlot(MAT_seurat, dims = 1:20)
ElbowPlot(MAT_seurat)

MAT_seurat <- FindNeighbors(MAT_seurat, dims = 1:20)
MAT_seurat <- FindClusters(MAT_seurat, resolution = 0.5)

MAT_seurat <- RunTSNE(MAT_seurat, dims = 1:20)

DimPlot(MAT_seurat, reduction = "tsne")


MAT<-as.matrix(MAT_seurat@assays$RNA@counts)
MAT<-sweep(MAT,2,colSums(MAT),FUN = "/")*1000000
MAT_ident<-as.vector(MAT_seurat@meta.data[,6])
Cell_anno<-Cell_anno[colnames(MAT)]
save(MAT_ident,MAT_seurat,Cell_anno,file = "BC01_use.Rdata")

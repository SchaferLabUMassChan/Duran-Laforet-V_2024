#####Analysis run with R v.4.2.2######

#load packages
library(Seurat)
library(dplyr)
library(Matrix)
library(reticulate)
library(scCustomize)
library(ggplot2)
library(sp)

### Setup object ###

## Read in data. Transpose so rows are genes and columns are cells. Remove "blank" genes. 
## Calculate total transcripts per cell.  Remove cells with ≤ 40 transcripts or volume < 50 µm3. 
## Add metadata for volume, x position, y position, condition, sample and sex. Normalize to the total transcripts/cell for each coverslip

H3777 <- read.table(file='merged_metadata_and_partitions_H3777.csv', sep=",", header=TRUE,row.names=1)
H3777 <- t(H3777)
H3777.Genes <- H3777[1:400,]
H3777[467,] <- colSums(H3777.Genes)
Real.Cells <- which(H3777[467,] > 10)
H3777 <- H3777[,Real.Cells]
Big.Cells <- which(H3777[464,]>50)
H3777 <- H3777[,Big.Cells]
volume <- H3777[464,]
head(volume)
center.x <- H3777[465,]
head(center.x)
FOV <- H3777[463,]
head(FOV)
center.y <- H3777[466,]
head(center.y)
mean.RNA <-mean(H3777[467,])
head(mean.RNA)
H3777 <- H3777[1:400,]
H3777 <- H3777/mean.RNA
H3777 <- H3777/volume
H3777_s <-CreateSeuratObject(counts=H3777)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
H3777_s@meta.data <- cbind(H3777_s@meta.data,Volumes)
H3777_s@meta.data <- cbind(H3777_s@meta.data,Center.x)
H3777_s@meta.data <- cbind(H3777_s@meta.data,Center.y)

#create metadata for conditions
H3777_s$condition <- "Aged"
H3777_s$sample <- "3777"
H3777_s$sex <- "M"


H3480 <- read.table(file='merged_metadata_and_partitions_H3480.csv', sep=",", header=TRUE,row.names=1)
H3480 <- t(H3480)
H3480.Genes <- H3480[1:400,]
H3480[467,] <- colSums(H3480.Genes)
Real.Cells <- which(H3480[467,] > 10)
H3480 <- H3480[,Real.Cells]
Big.Cells <- which(H3480[464,]>50)
H3480 <- H3480[,Big.Cells]
volume <- H3480[464,]
head(volume)
center.x <- H3480[465,]
head(center.x)
FOV <- H3480[463,]
head(FOV)
center.y <- H3480[466,]
head(center.y)
mean.RNA <-mean(H3480[467,])
head(mean.RNA)
H3480 <- H3480[1:400,]
H3480 <- H3480/mean.RNA
H3480 <- H3480/volume
H3480_s <-CreateSeuratObject(counts=H3480)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
H3480_s@meta.data <- cbind(H3480_s@meta.data,Volumes)
H3480_s@meta.data <- cbind(H3480_s@meta.data,Center.x)
H3480_s@meta.data <- cbind(H3480_s@meta.data,Center.y)

#create metadata for conditions
H3480_s$condition <- "AD"
H3480_s$sample <- "3480"
H3480_s$sex <- "M"


H879 <- read.table(file='merged_metadata_and_partitions_H879.csv', sep=",", header=TRUE,row.names=1)
H879 <- t(H879)
H879.Genes <- H879[1:400,]
H879[467,] <- colSums(H879.Genes)
Real.Cells <- which(H879[467,] > 10)
H879 <- H879[,Real.Cells]
Big.Cells <- which(H879[464,]>50)
H879 <- H879[,Big.Cells]
volume <- H879[464,]
head(volume)
center.x <- H879[465,]
head(center.x)
FOV <- H879[463,]
head(FOV)
center.y <- H879[466,]
head(center.y)
mean.RNA <-mean(H879[467,])
head(mean.RNA)
H879 <- H879[1:400,]
H879 <- H879/mean.RNA
H879 <- H879/volume
H879_s <-CreateSeuratObject(counts=H879)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
H879_s@meta.data <- cbind(H879_s@meta.data,Volumes)
H879_s@meta.data <- cbind(H879_s@meta.data,Center.x)
H879_s@meta.data <- cbind(H879_s@meta.data,Center.y)

#create metadata for conditions
H879_s$condition <- "Young"
H879_s$sample <- "879"
H879_s$sex <- "M"


H2526 <- read.table(file='merged_metadata_and_partitions_H2526.csv', sep=",", header=TRUE,row.names=1)
H2526 <- t(H2526)
H2526.Genes <- H2526[1:400,]
H2526[467,] <- colSums(H2526.Genes)
Real.Cells <- which(H2526[467,] > 10)
H2526 <- H2526[,Real.Cells]
Big.Cells <- which(H2526[464,]>50)
H2526 <- H2526[,Big.Cells]
volume <- H2526[464,]
head(volume)
center.x <- H2526[465,]
head(center.x)
FOV <- H2526[463,]
head(FOV)
center.y <- H2526[466,]
head(center.y)
mean.RNA <-mean(H2526[467,])
head(mean.RNA)
H2526 <- H2526[1:400,]
H2526 <- H2526/mean.RNA
H2526 <- H2526/volume
H2526_s <-CreateSeuratObject(counts=H2526)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
H2526_s@meta.data <- cbind(H2526_s@meta.data,Volumes)
H2526_s@meta.data <- cbind(H2526_s@meta.data,Center.x)
H2526_s@meta.data <- cbind(H2526_s@meta.data,Center.y)

#create metadata for conditions
H2526_s$condition <- "Aged"
H2526_s$sample <- "2526"
H2526_s$sex <- "F"


H5125 <- read.table(file='merged_metadata_and_partitions_H5125.csv', sep=",", header=TRUE,row.names=1)
H5125 <- t(H5125)
H5125.Genes <- H5125[1:400,]
H5125[467,] <- colSums(H5125.Genes)
Real.Cells <- which(H5125[467,] > 10)
H5125 <- H5125[,Real.Cells]
Big.Cells <- which(H5125[464,]>50)
H5125 <- H5125[,Big.Cells]
volume <- H5125[464,]
head(volume)
center.x <- H5125[465,]
head(center.x)
FOV <- H5125[463,]
head(FOV)
center.y <- H5125[466,]
head(center.y)
mean.RNA <-mean(H5125[467,])
head(mean.RNA)
H5125 <- H5125[1:400,]
H5125 <- H5125/mean.RNA
H5125 <- H5125/volume
H5125_s <-CreateSeuratObject(counts=H5125)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
H5125_s@meta.data <- cbind(H5125_s@meta.data,Volumes)
H5125_s@meta.data <- cbind(H5125_s@meta.data,Center.x)
H5125_s@meta.data <- cbind(H5125_s@meta.data,Center.y)

#create metadata for conditions
H5125_s$condition <- "Young"
H5125_s$sample <- "5125"
H5125_s$sex <- "F"


H2498 <- read.table(file='merged_metadata_and_partitions_H2498.csv', sep=",", header=TRUE,row.names=1)
H2498 <- t(H2498)
H2498.Genes <- H2498[1:400,]
H2498[467,] <- colSums(H2498.Genes)
Real.Cells <- which(H2498[467,] > 10)
H2498 <- H2498[,Real.Cells]
Big.Cells <- which(H2498[464,]>50)
H2498 <- H2498[,Big.Cells]
volume <- H2498[464,]
head(volume)
center.x <- H2498[465,]
head(center.x)
FOV <- H2498[463,]
head(FOV)
center.y <- H2498[466,]
head(center.y)
mean.RNA <-mean(H2498[467,])
head(mean.RNA)
H2498 <- H2498[1:400,]
H2498 <- H2498/mean.RNA
H2498 <- H2498/volume
H2498_s <-CreateSeuratObject(counts=H2498)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
H2498_s@meta.data <- cbind(H2498_s@meta.data,Volumes)
H2498_s@meta.data <- cbind(H2498_s@meta.data,Center.x)
H2498_s@meta.data <- cbind(H2498_s@meta.data,Center.y)

#create metadata for conditions
H2498_s$condition <- "AD"
H2498_s$sample <- "2498"
H2498_s$sex <- "F"


H3494 <- read.table(file='merged_metadata_and_partitions_H3494.csv', sep=",", header=TRUE,row.names=1)
H3494 <- t(H3494)
H3494.Genes <- H3494[1:400,]
H3494[467,] <- colSums(H3494.Genes)
Real.Cells <- which(H3494[467,] > 10)
H3494 <- H3494[,Real.Cells]
Big.Cells <- which(H3494[464,]>50)
H3494 <- H3494[,Big.Cells]
volume <- H3494[464,]
head(volume)
center.x <- H3494[465,]
head(center.x)
FOV <- H3494[463,]
head(FOV)
center.y <- H3494[466,]
head(center.y)
mean.RNA <-mean(H3494[467,])
head(mean.RNA)
H3494 <- H3494[1:400,]
H3494 <- H3494/mean.RNA
H3494 <- H3494/volume
H3494_s <-CreateSeuratObject(counts=H3494)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
H3494_s@meta.data <- cbind(H3494_s@meta.data,Volumes)
H3494_s@meta.data <- cbind(H3494_s@meta.data,Center.x)
H3494_s@meta.data <- cbind(H3494_s@meta.data,Center.y)

#create metadata for conditions
H3494_s$condition <- "AD"
H3494_s$sample <- "3494"
H3494_s$sex <- "F"


H255838 <- read.table(file='merged_metadata_and_partitions_H255838.csv', sep=",", header=TRUE,row.names=1)
H255838 <- t(H255838)
H255838.Genes <- H255838[1:400,]
H255838[467,] <- colSums(H255838.Genes)
Real.Cells <- which(H255838[467,] > 10)
H255838 <- H255838[,Real.Cells]
Big.Cells <- which(H255838[464,]>50)
H255838 <- H255838[,Big.Cells]
volume <- H255838[464,]
head(volume)
center.x <- H255838[465,]
head(center.x)
FOV <- H255838[463,]
head(FOV)
center.y <- H255838[466,]
head(center.y)
mean.RNA <-mean(H255838[467,])
head(mean.RNA)
H255838 <- H255838[1:400,]
H255838 <- H255838/mean.RNA
H255838 <- H255838/volume
H255838_s <-CreateSeuratObject(counts=H255838)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
H255838_s@meta.data <- cbind(H255838_s@meta.data,Volumes)
H255838_s@meta.data <- cbind(H255838_s@meta.data,Center.x)
H255838_s@meta.data <- cbind(H255838_s@meta.data,Center.y)

#create metadata for conditions
H255838_s$condition <- "Aged"
H255838_s$sample <- "255838"
H255838_s$sex <- "F"


#######main clustering######


#QC checks 
VlnPlot(all.combined.human.all, features=c('nFeature_RNA', 'nCount_RNA'), ncol=3, pt.size=0)

#remove non-cells
all.combined.human.all <- subset(all.combined.human.all,nFeature_RNA > 10)


VlnPlot(all.combined.human.all,features="nCount_RNA",group.by="sample",pt.size=0)


#normalize data
all.combined.human.all <- NormalizeData(all.combined.human.all)

#find variable features (reduce from default 2000 to 400 since 400 genes in the panel)

all.combined.human.all <- FindVariableFeatures(all.combined.human.all, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(all.combined.human.all)
all.combined.human.all <- ScaleData(all.combined.human.all, features = all.genes)
all.combined.human.all <- RunPCA(all.combined.human.all, features = VariableFeatures(object = all.combined.human.all),npcs = 50)
all.combined.human.all <- FindNeighbors(all.combined.human.all, dims = 1:21)
all.combined.human.all <- FindClusters(all.combined.human.all, resolution = 0.9)
all.combined.human.all <- RunUMAP(all.combined.human.all, dims = 1:21)

# cell-type annotation
Idents(all.combined.human.all) <- all.combined.human.all@meta.data$seurat_clusters
levels(all.combined.human.all)
new.cluster.ids <- c('OLs',
                     'OLs', 
                     'OLs',
                     'Endothelial Cells',
                     'Neurons',
                     'Astrocytes',
                     'Microglia',
                     'OPCs',
                     'Microglia',
                     'Astrocytes',
                     'OLs','Neurons', 
                     'Pericytes', 
                     'Neurons',
                     'Astrocytes',
                     'Neutrophils',
                     'OLs',
                     'Endothelial Cells',
                     'SMC', 
                     'T Cells',
                     'Neurons/OLs',
                     'T Cells',
                     'Neurons',
                     'OLs',
                     'Neurons',
                     'OLs',
                     'Mg/Neurons')

names(new.cluster.ids) <- levels(all.combined.human.all)
all.combined.human.all <- RenameIdents(all.combined.human.all, new.cluster.ids)
all.combined.human.all$CellType <- Idents(all.combined.human.all)
table(Idents(all.combined.human.all))


# spatial plots for cell type
for (i in "OLs") {
  cluster <- i
  sampleA <- "2526"
  x0 <- all.combined.human.all$center.x[WhichCells(object = subset(all.combined.human.all,subset = sample ==sampleA), ident = cluster)]
  y0 <- all.combined.human.all$center.y[WhichCells(object = subset(all.combined.human.all,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+ geom_point(color="#F8766D", size = 0.000001)+ggtitle(paste(sampleA,cluster,sep=" "))+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
}

for (i in "Microglia") {
  cluster <- i
  sampleA <- "2526"
  x0 <- all.combined.human.all$center.x[WhichCells(object = subset(all.combined.human.all,subset = sample ==sampleA), ident = cluster)]
  y0 <- all.combined.human.all$center.y[WhichCells(object = subset(all.combined.human.all,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+ geom_point(color="#0CB702", size = 0.000001)+ggtitle(paste(sampleA,cluster,sep=" "))+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
}


for (i in "Neurons") {
  cluster <- i
  sampleA <- "2526"
  x0 <- all.combined.human.all$center.x[WhichCells(object = subset(all.combined.human.all,subset = sample ==sampleA), ident = cluster)]
  y0 <- all.combined.human.all$center.y[WhichCells(object = subset(all.combined.human.all,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+ geom_point(color="#CD9600", size = 0.000001)+ggtitle(paste(sampleA,cluster,sep=" "))+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
}

for (i in "Astrocytes") {
  cluster <- i
  sampleA <- "2526"
  x0 <- all.combined.human.all$center.x[WhichCells(object = subset(all.combined.human.all,subset = sample ==sampleA), ident = cluster)]
  y0 <- all.combined.human.all$center.y[WhichCells(object = subset(all.combined.human.all,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+ geom_point(color="#ABA300", size = 0.000001)+ggtitle(paste(sampleA,cluster,sep=" "))+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
}

#spatial plots all cell types (repeat for each sample)
dfall<- data.frame(all.combined.human.all$center.x[WhichCells(object = subset(all.combined.human.all,subset = sample =="3777"))],all.combined.human.all$center.y[WhichCells(object = subset(all.combined.human.all,subset = sample =="3777"))],all.combined.human.all@active.ident[WhichCells(object = subset(all.combined.human.all,subset = sample =="3777"))])
x <- all.combined.human.all$center.x[WhichCells(object = subset(all.combined.human.all,subset = sample =="3777"))]
y <- all.combined.human.all$center.y[WhichCells(object = subset(all.combined.human.all,subset = sample =="3777"))]
group <- all.combined.human.all@active.ident[WhichCells(object = subset(all.combined.human.all,subset = sample =="3777"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()

#######only microglia#######


onlymicroglia <- subset(all.combined.human.all,idents = c("Microglia"))
onlymicroglia <-NormalizeData(onlymicroglia)
onlymicroglia <- FindVariableFeatures(onlymicroglia, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(onlymicroglia)
onlymicroglia <- ScaleData(onlymicroglia, features = all.genes)
onlymicroglia <- RunPCA(onlymicroglia, features = VariableFeatures(object = onlymicroglia),npcs = 100)
onlymicroglia <- FindNeighbors(onlymicroglia, dims = 1:20)
onlymicroglia <- FindClusters(onlymicroglia, resolution = 0.8)
onlymicroglia <- RunUMAP(onlymicroglia, dims = 1:20)

#remove hybrids
onlymicrogliaclean <- subset(onlymicroglia,idents = c("0", "1", "2", "3", "4", "5", "6", "11", "12"))
#re-cluster

onlymicrogliaclean <-NormalizeData(onlymicrogliaclean)
onlymicrogliaclean <- FindVariableFeatures(onlymicrogliaclean, selection.method = "vst", nfeatures = 400)

all.genes <- rownames(onlymicrogliaclean)
onlymicrogliaclean <- ScaleData(onlymicrogliaclean, features = all.genes)
onlymicrogliaclean <- RunPCA(onlymicrogliaclean, features = VariableFeatures(object = onlymicrogliaclean),npcs = 100)
onlymicrogliaclean <- FindNeighbors(onlymicrogliaclean, dims = 1:20)
onlymicrogliaclean <- FindClusters(onlymicrogliaclean, resolution = 0.8)
onlymicrogliaclean <- RunUMAP(onlymicrogliaclean, dims = 1:20)

#remove hybrids
onlymicrogliaclean <- subset(onlymicroglia,idents = c("0", "1", "2", "3", "4", "5", "9", "10", "11", "12"))

#re-cluster
onlymicrogliaclean <-NormalizeData(onlymicrogliaclean)
onlymicrogliaclean <- FindVariableFeatures(onlymicrogliaclean, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(onlymicrogliaclean)
onlymicrogliaclean <- ScaleData(onlymicrogliaclean, features = all.genes)
onlymicrogliaclean <- RunPCA(onlymicrogliaclean, features = VariableFeatures(object = onlymicrogliaclean),npcs = 100)
onlymicrogliaclean <- FindNeighbors(onlymicrogliaclean, dims = 1:20)
onlymicrogliaclean <- FindClusters(onlymicrogliaclean, resolution = 0.3)
onlymicrogliaclean <- RunUMAP(onlymicrogliaclean, dims = 1:20)
DimPlot(onlymicrogliaclean, reduction = "umap",label=TRUE)


all.markers.onlymicroglia.clean <- FindAllMarkers(onlymicrogliaclean, min.pct = 0.25, logfc.threshold = 0.25)


all.markers.onlymicroglia.clean %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
heatmap <- DoHeatmap(onlymicrogliaclean, features = top5$gene, cells = 1000) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))



senescence <- list(c(
                  'AXL',
                  'C3',
                  'CCL20',
                  'CCL3',
                  'CCL8',
                  'CD9',
                  'CSF1',
                  'CSF2',
                  'CTSB',
                  'CXCL1',
                  'CXCL10',
                  'CXCL3',
                  'CXCL8',
                  'DKK1',
                  'EDN1',
                  'ESM1',
                  'FGF1',
                  'FGF2',
                  'GDF15',
                  'HGF',
                  'HMGB1',
                  'IGFBP2',
                  'IGFBP3',
                  'IGFBP5',
                  'IL10',
                  'IL13',
                  'IL15',
                  'IL1A',
                  'IL1B',
                  'IL6',
                  'IL7',
                  'LCP1',
                  'MMP12',
                  'MMP13',
                  'MMP2',
                  'MMP3',
                  'MMP9',
                  'PECAM1',
                  'SEMA3F',
                  'SERPINE1',
                  'SPP1',
                  'TIMP2',
                  'TNF',
                  'VEGFC',
                  'VGF',
                  'WNT16'))

#calculate senescence module score
onlymicrogliaclean <- AddModuleScore(
  object = onlymicrogliaclean, 
  features = senescence, 
  ctrl = 5, 
  name = 'senescenceMg'
)


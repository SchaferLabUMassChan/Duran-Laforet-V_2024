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
## Calculate total transcripts per cell.  Remove cells with ≤ 40 transcripts or volume < 100 µm3. 
## Add metadata for volume, x position, y position, condition, sample and sex. Normalize to the total transcripts/cell for each coverslip

A1767 <- read.table(file='merged_metadata_and_partitions_Young_1.csv', sep=",", header=TRUE,row.names=1)
A1767 <- t(A1767)
A1767.Genes <- A1767[1:400,]
A1767[467,] <- colSums(A1767.Genes)
Real.Cells <- which(A1767[467,] > 40)
A1767 <- A1767[,Real.Cells]
Big.Cells <- which(A1767[464,]>100)
A1767 <- A1767[,Big.Cells]
volume <- A1767[464,]
head(volume)
center.x <- A1767[465,]
head(center.x)
FOV <- A1767[463,]
head(FOV)
center.y <- A1767[466,]
head(center.y)
mean.RNA <-mean(A1767[467,])
head(mean.RNA)
A1767 <- A1767[1:400,]
A1767 <- A1767/mean.RNA
A1767 <- A1767/volume
A1767_s <-CreateSeuratObject(counts=A1767)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
A1767_s@meta.data <- cbind(A1767_s@meta.data,Volumes)
A1767_s@meta.data <- cbind(A1767_s@meta.data,Center.x)
A1767_s@meta.data <- cbind(A1767_s@meta.data,Center.y)


A1766 <- read.table(file='merged_metadata_and_partitions_Young_2.csv', sep=",", header=TRUE,row.names=1)
A1766 <- t(A1766)
A1766.Genes <- A1766[1:400,]
A1766[467,] <- colSums(A1766.Genes)
Real.Cells <- which(A1766[467,] > 40)
A1766 <- A1766[,Real.Cells]
Big.Cells <- which(A1766[464,]>100)
A1766 <- A1766[,Big.Cells]
volume <- A1766[464,]
head(volume)
center.x <- A1766[465,]
head(center.x)
FOV <- A1766[463,]
head(FOV)
center.y <- A1766[466,]
head(center.y)
mean.RNA <-mean(A1766[467,])
head(mean.RNA)
A1766 <- A1766[1:400,]
A1766 <- A1766/mean.RNA
A1766 <- A1766/volume
A1766_s <-CreateSeuratObject(counts=A1766)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
A1766_s@meta.data <- cbind(A1766_s@meta.data,Volumes)
A1766_s@meta.data <- cbind(A1766_s@meta.data,Center.x)
A1766_s@meta.data <- cbind(A1766_s@meta.data,Center.y)


A1597 <- read.table(file='merged_metadata_and_partitions_Young_3.csv', sep=",", header=TRUE,row.names=1)
A1597 <- t(A1597)
A1597.Genes <- A1597[1:400,]
A1597[467,] <- colSums(A1597.Genes)
Real.Cells <- which(A1597[467,] > 40)
A1597 <- A1597[,Real.Cells]
Big.Cells <- which(A1597[464,]>100)
A1597 <- A1597[,Big.Cells]
volume <- A1597[464,]
head(volume)
center.x <- A1597[465,]
head(center.x)
FOV <- A1597[463,]
head(FOV)
center.y <- A1597[466,]
head(center.y)
mean.RNA <-mean(A1597[467,])
head(mean.RNA)
A1597 <- A1597[1:400,]
A1597 <- A1597/mean.RNA
A1597 <- A1597/volume
A1597_s <-CreateSeuratObject(counts=A1597)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
A1597_s@meta.data <- cbind(A1597_s@meta.data,Volumes)
A1597_s@meta.data <- cbind(A1597_s@meta.data,Center.x)
A1597_s@meta.data <- cbind(A1597_s@meta.data,Center.y)



A1645 <- read.table(file='merged_metadata_and_partitions_Aged_1.csv', sep=',', header=TRUE,row.names=1)
A1645 <- t(A1645)
A1645.Genes <- A1645[1:400,]
A1645[467,] <- colSums(A1645.Genes)
Real.Cells <- which(A1645[467,] > 40)
A1645 <- A1645[,Real.Cells]
Big.Cells <- which(A1645[464,]>100)
A1645 <- A1645[,Big.Cells]
volume <- A1645[464,]
head(volume)
center.x <- A1645[465,]
head(center.x)
FOV <- A1645[463,]
head(FOV)
center.y <- A1645[466,]
head(center.y)
mean.RNA <-mean(A1645[467,])
head(mean.RNA)
A1645 <- A1645[1:400,]
A1645 <- A1645/mean.RNA
A1645 <- A1645/volume
A1645_s <-CreateSeuratObject(counts=A1645)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
A1645_s@meta.data <- cbind(A1645_s@meta.data,Volumes)
A1645_s@meta.data <- cbind(A1645_s@meta.data,Center.x)
A1645_s@meta.data <- cbind(A1645_s@meta.data,Center.y)


A1419 <- read.table(file='merged_metadata_and_partitions_Aged_2.csv', sep=',', header=TRUE,row.names=1)
A1419 <- t(A1419)
A1419.Genes <- A1419[1:400,]
A1419[467,] <- colSums(A1419.Genes)
Real.Cells <- which(A1419[467,] > 40)
A1419 <- A1419[,Real.Cells]
Big.Cells <- which(A1419[464,]>100)
A1419 <- A1419[,Big.Cells]
volume <- A1419[464,]
head(volume)
center.x <- A1419[465,]
head(center.x)
FOV <- A1419[463,]
head(FOV)
center.y <- A1419[466,]
head(center.y)
mean.RNA <-mean(A1419[467,])
head(mean.RNA)
A1419 <- A1419[1:400,]
A1419 <- A1419/mean.RNA
A1419 <- A1419/volume
A1419_s <-CreateSeuratObject(counts=A1419)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
A1419_s@meta.data <- cbind(A1419_s@meta.data,Volumes)
A1419_s@meta.data <- cbind(A1419_s@meta.data,Center.x)
A1419_s@meta.data <- cbind(A1419_s@meta.data,Center.y)

A1643 <- read.table(file='merged_metadata_and_partitions_Aged_3.csv', sep=',', header=TRUE,row.names=1)
A1643 <- t(A1643)
A1643.Genes <- A1643[1:400,]
A1643[467,] <- colSums(A1643.Genes)
Real.Cells <- which(A1643[467,] > 40)
A1643 <- A1643[,Real.Cells]
Big.Cells <- which(A1643[464,]>100)
A1643 <- A1643[,Big.Cells]
volume <- A1643[464,]
head(volume)
center.x <- A1643[465,]
head(center.x)
FOV <- A1643[463,]
head(FOV)
center.y <- A1643[466,]
head(center.y)
mean.RNA <-mean(A1643[467,])
head(mean.RNA)
A1643 <- A1643[1:400,]
A1643 <- A1643/mean.RNA
A1643 <- A1643/volume
A1643_s <-CreateSeuratObject(counts=A1643)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
A1643_s@meta.data <- cbind(A1643_s@meta.data,Volumes)
A1643_s@meta.data <- cbind(A1643_s@meta.data,Center.x)
A1643_s@meta.data <- cbind(A1643_s@meta.data,Center.y)


AD5x5 <- read.table(file='merged_metadata_and_partitions_5xFAD_1.csv', sep=',', header=TRUE,row.names=1)
AD5x5 <- t(AD5x5)
AD5x5.Genes <- AD5x5[1:400,]
AD5x5[467,] <- colSums(AD5x5.Genes)
Real.Cells <- which(AD5x5[467,] > 40)
AD5x5 <- AD5x5[,Real.Cells]
Big.Cells <- which(AD5x5[464,]>100)
AD5x5 <- AD5x5[,Big.Cells]
volume <- AD5x5[464,]
head(volume)
center.x <- AD5x5[465,]
head(center.x)
FOV <- AD5x5[463,]
head(FOV)
center.y <- AD5x5[466,]
head(center.y)
mean.RNA <-mean(AD5x5[467,])
head(mean.RNA)
AD5x5 <- AD5x5[1:400,]
AD5x5 <- AD5x5/mean.RNA
AD5x5 <- AD5x5/volume
AD5x5_s <-CreateSeuratObject(counts=AD5x5)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
AD5x5_s@meta.data <- cbind(AD5x5_s@meta.data,Volumes)
AD5x5_s@meta.data <- cbind(AD5x5_s@meta.data,Center.x)
AD5x5_s@meta.data <- cbind(AD5x5_s@meta.data,Center.y)

AD5x4 <- read.table(file='merged_metadata_and_partitions_5xFAD_2.csv', sep=',', header=TRUE,row.names=1)
AD5x4 <- t(AD5x4)
AD5x4.Genes <- AD5x4[1:400,]
AD5x4[467,] <- colSums(AD5x4.Genes)
Real.Cells <- which(AD5x4[467,] > 40)
AD5x4 <- AD5x4[,Real.Cells]
Big.Cells <- which(AD5x4[464,]>100)
AD5x4 <- AD5x4[,Big.Cells]
volume <- AD5x4[464,]
head(volume)
center.x <- AD5x4[465,]
head(center.x)
FOV <- AD5x4[463,]
head(FOV)
center.y <- AD5x4[466,]
head(center.y)
mean.RNA <-mean(AD5x4[467,])
head(mean.RNA)
AD5x4 <- AD5x4[1:400,]
AD5x4 <- AD5x4/mean.RNA
AD5x4 <- AD5x4/volume
AD5x4_s <-CreateSeuratObject(counts=AD5x4)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
AD5x4_s@meta.data <- cbind(AD5x4_s@meta.data,Volumes)
AD5x4_s@meta.data <- cbind(AD5x4_s@meta.data,Center.x)
AD5x4_s@meta.data <- cbind(AD5x4_s@meta.data,Center.y)



AD5x2 <- read.table(file='merged_metadata_and_partitions_5xFAD_3.csv', sep=',', header=TRUE,row.names=1)
AD5x2 <- t(AD5x2)
AD5x2.Genes <- AD5x2[1:400,]
AD5x2[467,] <- colSums(AD5x2.Genes)
Real.Cells <- which(AD5x2[467,] > 40)
AD5x2 <- AD5x2[,Real.Cells]
Big.Cells <- which(AD5x2[464,]>100)
AD5x2 <- AD5x2[,Big.Cells]
volume <- AD5x2[464,]
head(volume)
center.x <- AD5x2[465,]
head(center.x)
FOV <- AD5x2[463,]
head(FOV)
center.y <- AD5x2[466,]
head(center.y)
mean.RNA <-mean(AD5x2[467,])
head(mean.RNA)
AD5x2 <- AD5x2[1:400,]
AD5x2 <- AD5x2/mean.RNA
AD5x2 <- AD5x2/volume
AD5x2_s <-CreateSeuratObject(counts=AD5x2)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
AD5x2_s@meta.data <- cbind(AD5x2_s@meta.data,Volumes)
AD5x2_s@meta.data <- cbind(AD5x2_s@meta.data,Center.x)
AD5x2_s@meta.data <- cbind(AD5x2_s@meta.data,Center.y)


WT5x12 <- read.table(file='merged_metadata_and_partitions_WT_5xFAD_1.csv', sep=',', header=TRUE,row.names=1)
WT5x12 <- t(WT5x12)
WT5x12.Genes <- WT5x12[1:400,]
WT5x12[467,] <- colSums(WT5x12.Genes)
Real.Cells <- which(WT5x12[467,] > 40)
WT5x12 <- WT5x12[,Real.Cells]
Big.Cells <- which(WT5x12[464,]>100)
WT5x12 <- WT5x12[,Big.Cells]
volume <- WT5x12[464,]
head(volume)
center.x <- WT5x12[465,]
head(center.x)
FOV <- WT5x12[463,]
head(FOV)
center.y <- WT5x12[466,]
head(center.y)
mean.RNA <-mean(WT5x12[467,])
head(mean.RNA)
WT5x12 <- WT5x12[1:400,]
WT5x12 <- WT5x12/mean.RNA
WT5x12 <- WT5x12/volume
WT5x12_s <-CreateSeuratObject(counts=WT5x12)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
WT5x12_s@meta.data <- cbind(WT5x12_s@meta.data,Volumes)
WT5x12_s@meta.data <- cbind(WT5x12_s@meta.data,Center.x)
WT5x12_s@meta.data <- cbind(WT5x12_s@meta.data,Center.y)


WT5x13 <- read.table(file='merged_metadata_and_partitions_WT_5xFAD_2.csv', sep=',', header=TRUE,row.names=1)
WT5x13 <- t(WT5x13)
WT5x13.Genes <- WT5x13[1:400,]
WT5x13[467,] <- colSums(WT5x13.Genes)
Real.Cells <- which(WT5x13[467,] > 40)
WT5x13 <- WT5x13[,Real.Cells]
Big.Cells <- which(WT5x13[464,]>100)
WT5x13 <- WT5x13[,Big.Cells]
volume <- WT5x13[464,]
head(volume)
center.x <- WT5x13[465,]
head(center.x)
FOV <- WT5x13[463,]
head(FOV)
center.y <- WT5x13[466,]
head(center.y)
mean.RNA <-mean(WT5x13[467,])
head(mean.RNA)
WT5x13 <- WT5x13[1:400,]
WT5x13 <- WT5x13/mean.RNA
WT5x13 <- WT5x13/volume
WT5x13_s <-CreateSeuratObject(counts=WT5x13)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
WT5x13_s@meta.data <- cbind(WT5x13_s@meta.data,Volumes)
WT5x13_s@meta.data <- cbind(WT5x13_s@meta.data,Center.x)
WT5x13_s@meta.data <- cbind(WT5x13_s@meta.data,Center.y)

WT5x35 <- read.table(file='merged_metadata_and_partitions_WT_5xFAD_3.csv', sep=',', header=TRUE,row.names=1)
WT5x35 <- t(WT5x35)
WT5x35.Genes <- WT5x35[1:400,]
WT5x35[467,] <- colSums(WT5x35.Genes)
Real.Cells <- which(WT5x35[467,] > 40)
WT5x35 <- WT5x35[,Real.Cells]
Big.Cells <- which(WT5x35[464,]>100)
WT5x35 <- WT5x35[,Big.Cells]
volume <- WT5x35[464,]
head(volume)
center.x <- WT5x35[465,]
head(center.x)
FOV <- WT5x35[463,]
head(FOV)
center.y <- WT5x35[466,]
head(center.y)
mean.RNA <-mean(WT5x35[467,])
head(mean.RNA)
WT5x35 <- WT5x35[1:400,]
WT5x35 <- WT5x35/mean.RNA
WT5x35 <- WT5x35/volume
WT5x35_s <-CreateSeuratObject(counts=WT5x35)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
WT5x35_s@meta.data <- cbind(WT5x35_s@meta.data,Volumes)
WT5x35_s@meta.data <- cbind(WT5x35_s@meta.data,Center.x)
WT5x35_s@meta.data <- cbind(WT5x35_s@meta.data,Center.y)

CFA_A2176 <- read.table(file='merged_metadata_and_partitions_CFA_1.csv', sep=',', header=TRUE,row.names=1)
CFA_A2176 <- t(CFA_A2176)
CFA_A2176.Genes <- CFA_A2176[1:400,]
CFA_A2176[467,] <- colSums(CFA_A2176.Genes)
Real.Cells <- which(CFA_A2176[467,] > 40)
CFA_A2176 <- CFA_A2176[,Real.Cells]
Big.Cells <- which(CFA_A2176[464,]>100)
CFA_A2176 <- CFA_A2176[,Big.Cells]
volume <- CFA_A2176[464,]
head(volume)
center.x <- CFA_A2176[465,]
head(center.x)
FOV <- CFA_A2176[463,]
head(FOV)
center.y <- CFA_A2176[466,]
head(center.y)
mean.RNA <-mean(CFA_A2176[467,])
head(mean.RNA)
CFA_A2176 <- CFA_A2176[1:400,]
CFA_A2176 <- CFA_A2176/mean.RNA
CFA_A2176 <- CFA_A2176/volume
CFA_A2176_s <-CreateSeuratObject(counts=CFA_A2176)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
CFA_A2176_s@meta.data <- cbind(CFA_A2176_s@meta.data,Volumes)
CFA_A2176_s@meta.data <- cbind(CFA_A2176_s@meta.data,Center.x)
CFA_A2176_s@meta.data <- cbind(CFA_A2176_s@meta.data,Center.y)

CFA_A2166 <- read.table(file='merged_metadata_and_partitions_CFA_2.csv', sep=',', header=TRUE,row.names=1)
CFA_A2166 <- t(CFA_A2166)
CFA_A2166.Genes <- CFA_A2166[1:400,]
CFA_A2166[467,] <- colSums(CFA_A2166.Genes)
Real.Cells <- which(CFA_A2166[467,] > 40)
CFA_A2166 <- CFA_A2166[,Real.Cells]
Big.Cells <- which(CFA_A2166[464,]>100)
CFA_A2166 <- CFA_A2166[,Big.Cells]
volume <- CFA_A2166[464,]
head(volume)
center.x <- CFA_A2166[465,]
head(center.x)
FOV <- CFA_A2166[463,]
head(FOV)
center.y <- CFA_A2166[466,]
head(center.y)
mean.RNA <-mean(CFA_A2166[467,])
head(mean.RNA)
CFA_A2166 <- CFA_A2166[1:400,]
CFA_A2166 <- CFA_A2166/mean.RNA
CFA_A2166 <- CFA_A2166/volume
CFA_A2166_s <-CreateSeuratObject(counts=CFA_A2166)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
CFA_A2166_s@meta.data <- cbind(CFA_A2166_s@meta.data,Volumes)
CFA_A2166_s@meta.data <- cbind(CFA_A2166_s@meta.data,Center.x)
CFA_A2166_s@meta.data <- cbind(CFA_A2166_s@meta.data,Center.y)

CFA_A2154 <- read.table(file='merged_metadata_and_partitions_CFA_3.csv', sep=',', header=TRUE,row.names=1)
CFA_A2154 <- t(CFA_A2154)
CFA_A2154.Genes <- CFA_A2154[1:400,]
CFA_A2154[467,] <- colSums(CFA_A2154.Genes)
Real.Cells <- which(CFA_A2154[467,] > 40)
CFA_A2154 <- CFA_A2154[,Real.Cells]
Big.Cells <- which(CFA_A2154[464,]>100)
CFA_A2154 <- CFA_A2154[,Big.Cells]
volume <- CFA_A2154[464,]
head(volume)
center.x <- CFA_A2154[465,]
head(center.x)
FOV <- CFA_A2154[463,]
head(FOV)
center.y <- CFA_A2154[466,]
head(center.y)
mean.RNA <-mean(CFA_A2154[467,])
head(mean.RNA)
CFA_A2154 <- CFA_A2154[1:400,]
CFA_A2154 <- CFA_A2154/mean.RNA
CFA_A2154 <- CFA_A2154/volume
CFA_A2154_s <-CreateSeuratObject(counts=CFA_A2154)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
CFA_A2154_s@meta.data <- cbind(CFA_A2154_s@meta.data,Volumes)
CFA_A2154_s@meta.data <- cbind(CFA_A2154_s@meta.data,Center.x)
CFA_A2154_s@meta.data <- cbind(CFA_A2154_s@meta.data,Center.y)

EAE_2163 <- read.table(file='merged_metadata_and_partitions_EAE_1.csv', sep=',', header=TRUE,row.names=1)
EAE_2163 <- t(EAE_2163)
EAE_2163.Genes <- EAE_2163[1:400,]
EAE_2163[467,] <- colSums(EAE_2163.Genes)
Real.Cells <- which(EAE_2163[467,] > 40)
EAE_2163 <- EAE_2163[,Real.Cells]
Big.Cells <- which(EAE_2163[464,]>100)
EAE_2163 <- EAE_2163[,Big.Cells]
volume <- EAE_2163[464,]
head(volume)
center.x <- EAE_2163[465,]
head(center.x)
FOV <- EAE_2163[463,]
head(FOV)
center.y <- EAE_2163[466,]
head(center.y)
mean.RNA <-mean(EAE_2163[467,])
head(mean.RNA)
EAE_2163 <- EAE_2163[1:400,]
EAE_2163 <- EAE_2163/mean.RNA
EAE_2163 <- EAE_2163/volume
EAE_2163_s <-CreateSeuratObject(counts=EAE_2163)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
EAE_2163_s@meta.data <- cbind(EAE_2163_s@meta.data,Volumes)
EAE_2163_s@meta.data <- cbind(EAE_2163_s@meta.data,Center.x)
EAE_2163_s@meta.data <- cbind(EAE_2163_s@meta.data,Center.y)

EAE_2161 <- read.table(file='merged_metadata_and_partitions_EAE_2.csv', sep=',', header=TRUE,row.names=1)
EAE_2161 <- t(EAE_2161)
EAE_2161.Genes <- EAE_2161[1:400,]
EAE_2161[467,] <- colSums(EAE_2161.Genes)
Real.Cells <- which(EAE_2161[467,] > 40)
EAE_2161 <- EAE_2161[,Real.Cells]
Big.Cells <- which(EAE_2161[464,]>100)
EAE_2161 <- EAE_2161[,Big.Cells]
volume <- EAE_2161[464,]
head(volume)
center.x <- EAE_2161[465,]
head(center.x)
FOV <- EAE_2161[463,]
head(FOV)
center.y <- EAE_2161[466,]
head(center.y)
mean.RNA <-mean(EAE_2161[467,])
head(mean.RNA)
EAE_2161 <- EAE_2161[1:400,]
EAE_2161 <- EAE_2161/mean.RNA
EAE_2161 <- EAE_2161/volume
EAE_2161_s <-CreateSeuratObject(counts=EAE_2161)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
EAE_2161_s@meta.data <- cbind(EAE_2161_s@meta.data,Volumes)
EAE_2161_s@meta.data <- cbind(EAE_2161_s@meta.data,Center.x)
EAE_2161_s@meta.data <- cbind(EAE_2161_s@meta.data,Center.y)

EAE_2565 <- read.table(file='merged_metadata_and_partitions_EAE_3.csv', sep=',', header=TRUE,row.names=1)
EAE_2565 <- t(EAE_2565)
EAE_2565.Genes <- EAE_2565[1:400,]
EAE_2565[467,] <- colSums(EAE_2565.Genes)
Real.Cells <- which(EAE_2565[467,] > 40)
EAE_2565 <- EAE_2565[,Real.Cells]
Big.Cells <- which(EAE_2565[464,]>100)
EAE_2565 <- EAE_2565[,Big.Cells]
volume <- EAE_2565[464,]
head(volume)
center.x <- EAE_2565[465,]
head(center.x)
FOV <- EAE_2565[463,]
head(FOV)
center.y <- EAE_2565[466,]
head(center.y)
mean.RNA <-mean(EAE_2565[467,])
head(mean.RNA)
EAE_2565 <- EAE_2565[1:400,]
EAE_2565 <- EAE_2565/mean.RNA
EAE_2565 <- EAE_2565/volume
EAE_2565_s <-CreateSeuratObject(counts=EAE_2565)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
EAE_2565_s@meta.data <- cbind(EAE_2565_s@meta.data,Volumes)
EAE_2565_s@meta.data <- cbind(EAE_2565_s@meta.data,Center.x)
EAE_2565_s@meta.data <- cbind(EAE_2565_s@meta.data,Center.y)

#create metadata for conditions
A1767_s$condition <- "Young"
A1767_s$sample <- "Young_1"
A1767_s$sex <- "M"


A1766_s$condition <- "Young"
A1766_s$sample <- "Young_2"
A1766_s$sex <- "M"


A1597_s$condition <- "Young"
A1597_s$sample <- "Young_3"
A1597_s$sex <- "F"


A1645_s$condition <- "Aged"
A1645_s$sample <- "Aged_1"
A1645_s$sex <- "F"


A1419_s$condition <- "Aged"
A1419_s$sample <- "Aged_2"
A1419_s$sex <- "M"


A1643_s$condition <- "Old"
A1643_s$sample <- "Aged_3"
A1643_s$sex <- "F"


AD5x5_s$condition <- "5xFAD"
AD5x5_s$sample <- "5xFAD_1"
AD5x5_s$sex <- "M"


AD5x4_s$condition <- "5xFAD"
AD5x4_s$sample <- "5xFAD_2"
AD5x4_s$sex <- "M"


AD5x2_s$condition <- "5xFAD"
AD5x2_s$sample <- "5xFAD_3"
AD5x2_s$sex <- "F"
AD5x2_s$area <- "54.865"

WT5x12_s$condition <- "WT-5xFAD"
WT5x12_s$sample <- "WT-5xFAD_1"
WT5x12_s$sex <- "F"


WT5x13_s$condition <- "WT-5xFAD"
WT5x13_s$sample <- "WT-5xFAD_2"
WT5x13_s$sex <- "F"


WT5x35_s$condition <- "WT-5xFAD"
WT5x35_s$sample <- "WT-5xFAD_3"
WT5x35_s$sex <- "M"


CFA_A2166_s$condition <- "CFA"
CFA_A2166_s$sample <- "CFA_1"
CFA_A2166_s$sex <- "F"


CFA_A2154_s$condition <- "CFA"
CFA_A2154_s$sample <- "CFA_2"
CFA_A2154_s$sex <- "F"


CFA_A2176_s$condition <- "CFA"
CFA_A2176_s$sample <- "CFA_3"
CFA_A2176_s$sex <- "F"


EAE_2163_s$condition <- "EAE"
EAE_2163_s$sample <- "EAE_1"
EAE_2163_s$sex <- "F"


EAE_2161_s$condition <- "EAE"
EAE_2161_s$sample <- "EAE_2"
EAE_2161_s$sex <- "F"


EAE_2565_s$condition <- "EAE"
EAE_2565_s$sample <- "EAE_3"
EAE_2565_s$sex <- "F"


#merge into single object 
all.combined <- merge(A1767_s, y = c(A1766_s,A1597_s,A1645_s,A1419_s,A1643_s, AD5x5_s, AD5x4_s, AD5x2_s, WT5x12_s, WT5x13_s, WT5x35_s,CFA_A2166_s,CFA_A2154_s,CFA_A2176_s,EAE_2163_s, EAE_2161_s, EAE_2565_s), add.cell.ids = c("2mo_1", "2mo_2","2mo_3","24mo_1","24mo_2","24mo_3", "5xFAD_1", "CFA_1", "5xFAD_3", "9mo_1", "9mo_2", "9mo_3", "CFA_1", "CFA_2", "CFA_3", "EAE_1", "EAE_2", "EAE_3"))

#QC checks 
VlnPlot(all.combined, features=c('nFeature_RNA', 'nCount_RNA'), ncol=3, pt.size=0)
VlnPlot(all.combined, features=c('nCount_RNA'), ncol=3, pt.size=0, split.by = 'sample')
#remove non-cells
all.combined <- subset(all.combined,nFeature_RNA > 10)

#all.combined.cropped <- subset(all.combined,nCount_RNA > 0.00075)

VlnPlot(all.combined,features="nCount_RNA",group.by="sample",pt.size=0)

########### select areas for analysis (remove non-wanted areas) ###########

#split object into individual samples
Young1 <- subset(x = all.combined, subset = (sample == "Young_1"))
Young2 <- subset(x = all.combined, subset = (sample == "Young_2"))
Young3 <- subset(x = all.combined, subset = (sample == "Young_3"))
Aged1 <- subset(x = all.combined, subset = (sample == "Aged_1"))
Aged2 <- subset(x = all.combined, subset = (sample == "Aged_2"))
Aged3 <- subset(x = all.combined, subset = (sample == "Aged_3"))
WT1 <- subset(x = all.combined, subset = (sample == "WT-5xFAD_1"))
WT2 <- subset(x = all.combined, subset = (sample == "WT-5xFAD_2"))
WT3 <- subset(x = all.combined, subset = (sample == "WT-5xFAD_3"))
AD1 <- subset(x = all.combined, subset = (sample == "5xFAD_1"))
AD2 <- subset(x = all.combined, subset = (sample == "5xFAD_2"))
AD3 <- subset(x = all.combined, subset = (sample == "5xFAD_3"))
CFA1 <- subset(x = all.combined, subset = (sample == "CFA_1"))
CFA2 <- subset(x = all.combined, subset = (sample == "CFA_2"))
CFA3 <- subset(x = all.combined, subset = (sample == "CFA_3"))
EAE1 <- subset(x = all.combined, subset = (sample == "EAE_1"))
EAE2 <- subset(x = all.combined, subset = (sample == "EAE_2"))
EAE3 <- subset(x = all.combined, subset = (sample == "EAE_3"))

#store x and y positions 
Aged1.x <- Aged1@meta.data[,5]
head(Aged1.x)
Aged1.y <- Aged1@meta.data[,6]
head(Aged1.y)

#regional polygon coordinates for selected brain region
Aged1_br_x <- c(3617.77,	7381.66,	10790.8,	10468.1,	8284.85,	6585.75,	5064.76,	4219.62,	3546.22,	2824.26,	1100.73,	1.77012,	20.9142,	823.077,	2353.97)
Aged1_br_y <- c(268.486,	852.694,	3620.82,	7115.93,	8298.64,	6954.63,	6894.63,	6564.76,	5983.67,	5758.33,	6614.24,	4787.76,	3186.64,	1790.13,	430.528)

#find region
Aged1_pol <- point.in.polygon(Aged1.x, Aged1.y, Aged1_br_x, Aged1_br_y)
Aged1.info <- cbind.data.frame(Aged1_pol)
head(Aged1.info, 50)

#assign cells to region
Aged1.info <- data.frame(Aged1.info)
Aged1@meta.data <- cbind(Aged1@meta.data,Aged1.info)
Aged1_br <- subset(Aged1,subset = Aged1_pol == 1)
Aged1_nobr <- subset(Aged1,subset = Aged1_pol == 0)
Aged1_br$Br <- "Yes"
Aged1_nobr$Br <- "No"
#check selection
dfall<- data.frame(Aged1_br$center.x[WhichCells(object = subset(Aged1_br,subset = sample =="Aged_1"))],Aged1_br$center.y[WhichCells(object = subset(Aged1_br,subset = sample =="Aged_1"))],Aged1_br@active.ident[WhichCells(object = subset(Aged1_br,subset = sample =="Aged_1"))])
x <- Aged1_br$center.x[WhichCells(object = subset(Aged1_br,subset = sample =="Aged_1"))]
y <- Aged1_br$center.y[WhichCells(object = subset(Aged1_br,subset = sample =="Aged_1"))]
group <- Aged1_br@active.ident[WhichCells(object = subset(Aged1_br,subset = sample =="Aged_1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()


#store x and y positions
Aged3.x <- Aged3@meta.data[,5]
head(Aged3.x)
Aged3.y <- Aged3@meta.data[,6]
head(Aged3.y)

#regional polygon coordinates for selected brain region
Aged3_br_x <- c(7149.28,	5507.05,	3027.68,	2081.15,	1597.77,	1996.51,	1844.37,	980.19,	1409.96,	2965.07,	6237.71,	8424.86,	8867.96,	7279.81)
Aged3_br_y <- c(9795.1,	10812.7,	10712,	9971.39,	8876.21,	6758.69,	4501.89,	2563.35,	1100.38,	33.48,	96.1328,	2287.27,	7135.5,	9894.69)

#find region
Aged3_pol <- point.in.polygon(Aged3.x, Aged3.y, Aged3_br_x, Aged3_br_y)
Aged3.info <- cbind.data.frame(Aged3_pol)
head(Aged3.info, 50)


#assign cells to region
Aged3.info <- data.frame(Aged3.info)
Aged3@meta.data <- cbind(Aged3@meta.data,Aged3.info)
Aged3_br <- subset(Aged3,subset = Aged3_pol == 1)
Aged3_nobr <- subset(Aged3,subset = Aged3_pol == 0)
Aged3_br$Br <- "Yes"
Aged3_nobr$Br <- "No"

#check selection
dfall<- data.frame(Aged3_br$center.x[WhichCells(object = subset(Aged3_br,subset = sample =="Aged_3"))],Aged3_br$center.y[WhichCells(object = subset(Aged3_br,subset = sample =="Aged_3"))],Aged3_br@active.ident[WhichCells(object = subset(Aged3_br,subset = sample =="Aged_3"))])
x <- Aged3_br$center.x[WhichCells(object = subset(Aged3_br,subset = sample =="Aged_3"))]
y <- Aged3_br$center.y[WhichCells(object = subset(Aged3_br,subset = sample =="Aged_3"))]
group <- Aged3_br@active.ident[WhichCells(object = subset(Aged3_br,subset = sample =="Aged_3"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#create final object
all.combined.cropped <-merge(Young1, y = c(Young2, Young3, WT1, WT2, WT3, Aged1_br, Aged2, Aged3_br, AD1, AD2, AD3, CFA1, CFA2, CFA3, EAE1, EAE2, EAE3))


#### clustering ######

#normalize data
all.combined.cropped <- NormalizeData(all.combined.cropped)

#find variable features (reduce from default 2000 to 400 since 400 genes in the panel)

all.combined.cropped <- FindVariableFeatures(all.combined.cropped, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(all.combined.cropped)
all.combined.cropped <- ScaleData(all.combined.cropped, features = all.genes)
all.combined.cropped <- RunPCA(all.combined.cropped, features = VariableFeatures(object = all.combined.cropped),npcs = 50)

# test of number of significant Principal Components
all.combined.cropped <- JackStraw(all.combined.cropped, num.replicate = 100,dims=40)
all.combined.cropped <- ScoreJackStraw(all.combined.cropped, dims = 1:40)
JackStrawPlot(all.combined.cropped, dims = 1:40)
#36 dims for further analysis 

all.combined.cropped <- FindNeighbors(all.combined.cropped, dims = 1:36)
all.combined.cropped <- FindClusters(all.combined.cropped, resolution = 0.8)
all.combined.cropped <- RunUMAP(all.combined.cropped, dims = 1:36)

#remove cluster 19 & 36 (imaging artifacts)
all.combined.cropped.clean <- subset(all.combined.cropped, idents = c("0",  "1" , "2",  "3" , "4" , "5" , "6",  "7" , "8",  "9" , "10", "11", "12", "13", "14", "15" ,"16" ,"17", "18",  "20", "21", "22", "23" ,"24", "25" ,"26", "27", "28",
                                                                      "29", "30", "31" ,"32", "33", "34", "35", "37" ,"38", "39", "40" ,"41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54" ,"55", "56", "57",
                                                                      "58", "59", "60", "61", "62", "63", "64", "65", "66" ,"67" ,"68" ,"69" ,"70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82" ,"83", "84", "85", "86",
                                                                      "87", "88"))

# cell-type annotation
Idents(all.combined.cropped.clean) <- all.combined.cropped.clean@meta.data$seurat_clusters
levels(all.combined.cropped.clean)
new.cluster.ids <- c('Neurons Ctx L1-4',
                     'Astrocytes',
                     'Oligodendrocytes',
                     'Neurons Ctx L5-6',
                     'Endothelial Cells',
                     'Inh. Neurons',
                     'Oligodendrocytes',
                     'Exc. Neurons Hypoth & Amyg', 
                     'Exc. Neurons Thalamus', 
                     'Microglia',
                     'Neurons Dentate Gyrus',
                     'OPCs',
                     'Inh. Neurons Caudoputamen',
                     'Neurons CA1-CA3 Hippo',
                     'Pericytes',
                     'Meninges', 
                     'Endothelial Cells',
                     'Sst+ INs',
                     'Oligodendrocytes',
                     'Neurons',
                     'Choroid Plexus',
                     'Endothelial Cells',
                     'Ependymal Cells',
                     'Habenula Neurons',
                     'Oligodendrocytes',
                     'Oligodendrocytes',
                     'Astrocytes',
                     'Microglia',
                     'Inh.Neurons L1-4',
                     'Astrocytes',
                     'Oligodendrocytes',
                     'Ependymal Cells',
                     'OPCs',
                     'Exc. Neurons Hypoth & Amyg',
                     'Ependymal Cells',
                     'Astrocytes',
                     'Neurons Ctx L1-4',
                     'Oligodendrocytes',
                     'Pericytes',
                     'Astrocytes',
                     'Neurons Ctx L1-4',
                     'Neurons Ctx L1-4',
                     'Astrocytes',
                     'Astrocytes',
                     'Oligodendrocytes',
                     'Endothelial Cells',
                     'Neurons Ctx L5-6', 
                     'Oligodendrocytes', 
                     'Neurons Ctx L5-6', 
                     'Inh. Neurons', 
                     'Oligodendrocytes', 
                     'Astrocytes', 
                     'Oligodendrocytes', 
                     'Astrocytes', 
                     'Neurons Ctx L1-4', 
                     'Neurons Ctx L1-4', 
                     'Neurons Ctx L5-6', 
                     'Habenula Neurons', 
                     'Sst+ INs', 
                     'Oligodendrocytes', 
                     'Sst+ INs', 
                     'Neurons Dentate Gyrus', 
                     'Neurons Ctx L1-4',
                     'Neurons Ctx L1-4',
                     'Endothelial Cells', 
                     'Neurons Ctx L5-6', 
                     'Neurons Ctx L1-4', 
                     'Neurons Ctx L5-6', 
                     'Neurons Ctx L1-4',
                     'Inh. Neurons', 
                     'Neurons', 
                     'Neurons',
                     'Exc. Neurons Hypoth & Amyg', 
                     'Oligodendrocytes', 
                     'Inh. Neurons Caudoputamen', 
                     'Neurons Ctx L1-4',
                     'Neurons Ctx L1-4',
                     'Neurons Ctx L1-4',
                     'Neurons Ctx L1-4',
                     'Astrocytes',
                     'Oligodendrocytes',
                     'Neurons Ctx L1-4',
                     'Exc. Neurons Hypoth & Amyg', 
                     'Exc. Neurons Hypoth & Amyg', 
                     'Neurons Ctx L1-4',
                     'Exc. Neurons Hypoth & Amyg', 
                     'Exc. Neurons Hypoth & Amyg'
)

names(new.cluster.ids) <- levels(all.combined.cropped.clean)
all.combined.cropped.clean <- RenameIdents(all.combined.cropped.clean, new.cluster.ids)
all.combined.cropped.clean$CellType <- Idents(all.combined.cropped.clean)

#spatial plots (repeat for each sample)
dfall<- data.frame(all.combined.cropped.clean$center.x[WhichCells(object = subset(all.combined.cropped.clean,subset = sample =="CFA_1"))],all.combined.cropped.clean$center.y[WhichCells(object = subset(all.combined.cropped.clean,subset = sample =="CFA_1"))],all.combined.cropped.clean@active.ident[WhichCells(object = subset(all.combined.cropped.clean,subset = sample =="CFA_1"))])
x <- all.combined.cropped.clean$center.x[WhichCells(object = subset(all.combined.cropped.clean,subset = sample =="CFA_1"))]
y <- all.combined.cropped.clean$center.y[WhichCells(object = subset(all.combined.cropped.clean,subset = sample =="CFA_1"))]
group <- all.combined.cropped.clean@active.ident[WhichCells(object = subset(all.combined.cropped.clean,subset = sample =="CFA_1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=0.000001, stroke=0)+ scale_color_manual(values=palette) +theme_classic()+coord_fixed()+NoLegend()

#calculate senescence module score#

senescence <- list(c("Axl","C3","Ccl20","Ccl3","Ccl8","Cd9","Csf1","Csf2",
                     "Ctsb","Cxcl1","Cxcl10","Cxcl2","Cxcl3","Dkk1","Edn1",
                     "Esm1","Fgf1","Fgf2","Fgf7","Gdf15","Hgf","Hmgb1","Igfbp2",
                     "Igfbp3","Igfbp5","Il10","Il13","Il15","Il1a","Il1b","Il6",
                     "Il7","Lcp1","Mmp12","Mmp13","Mmp2","Mmp3","Mmp9","Pecam1",
                     "Sema3f","Serpine1","Spp1","Spx","Timp2","Tnf","Vegfc","Vgf",
                     "Wnt16"
))

all.combined.cropped.clean <- AddModuleScore(
  object = all.combined.cropped.clean, 
  features = senescence, 
  ctrl = 5, 
  name = 'senescence'
)

# calculate threshold
young <- subset(x = all.combined.cropped.clean, subset = (condition == "Young"))
mean(young$senescence1)
sd(young$senescence1)
threesd <- (sd(young$senescence1) * 3)
threesd
cutoff <- mean(young$senescence1)+threesd
cutoff


colour_breaks <- c( 0.4138731, 25)

colours <- c( "yellow","darkblue" )

#create new metadata column
all.combined.cropped.clean$senescence <- "non senescent"
all.combined.cropped.clean$senescence[WhichCells(all.combined.cropped.clean, expression =  senescence1 > 0.4138731)] <- "senescent"
Idents(all.combined.cropped.clean) <- "senescence"
levels(all.combined.cropped.clean)

# Sescent module score feature plots (repeat for each sample)
FeaturePlot_scCustom(seurat_object = subset(all.combined.cropped.clean,subset = condition =='5xFAD') ,features = "senescence1", label = FALSE, split.by = "condition", order = T)+
  scale_colour_gradientn(
    limits  = range(all.combined.cropped.clean$senescence1),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = range(all.combined.cropped.clean$senescence1)), 1),
  )

#prepare new object for spatial plot
all.combined.cropped.clean.nonSenescent <- subset(all.combined.cropped.clean,subset = senescence1 < 0.4138731)
table(Idents(all.combined.cropped.clean.nonSenescent))
all.combined.cropped.clean.Senescent <- subset(all.combined.cropped.clean,subset = senescence1 >0.4138731)
table(Idents(all.combined.cropped.clean.Senescent))
all.combined.cropped.clean.Senescent$Senescent <- "Yes"
all.combined.cropped.clean.nonSenescent$Senescent <- "No"
all.combined.cropped.clean.Senescent_nonSenescent <- merge(all.combined.cropped.clean.nonSenescent,y = c(all.combined.cropped.clean.Senescent))

Idents(all.combined.cropped.clean.Senescent_nonSenescent)<-all.combined.cropped.clean.Senescent_nonSenescent$Senescent
levels(all.combined.cropped.clean.Senescent_nonSenescent)

#spatial plot
dfall<- data.frame(all.combined.cropped.clean.Senescent_nonSenescent$center.x[WhichCells(object = subset(all.combined.cropped.clean.Senescent_nonSenescent, sample =="Young_2"))],all.combined.cropped.clean.Senescent_nonSenescent$center.y[WhichCells(object = subset(all.combined.cropped.clean.Senescent_nonSenescent,subset = sample =="Young_2"))],all.combined.cropped.clean.Senescent_nonSenescent@active.ident[WhichCells(object = subset(all.combined.cropped.clean.Senescent_nonSenescent,subset = sample =="Young_2"))])
x <- all.combined.cropped.clean.Senescent_nonSenescent$center.x[WhichCells(object = subset(all.combined.cropped.clean.Senescent_nonSenescent,subset = sample =="Young_2"))]
y <- all.combined.cropped.clean.Senescent_nonSenescent$center.y[WhichCells(object = subset(all.combined.cropped.clean.Senescent_nonSenescent,subset = sample =="Young_2"))]
group <- all.combined.cropped.clean.Senescent_nonSenescent@active.ident[WhichCells(object = subset(all.combined.cropped.clean.Senescent_nonSenescent,subset = sample =="Young_2"))]



p <- ggplot(dfall ,aes(x,y,colour=group))+geom_point(size=0.000000001)+  theme_classic()
p + scale_color_manual(values=c( "lightgray","darkblue"))





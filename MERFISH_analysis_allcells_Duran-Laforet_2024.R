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

# Senescent module score feature plots (repeat for each sample)
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


##### only microglia ######

Idents(all.combined.cropped.clean) <- all.combined.cropped.clean@meta.data$CellType
levels(all.combined.cropped.clean)

# subset only microglia from complete dataset
onlymicroglia <- subset(all.combined.cropped.clean,idents = c("Microglia"))

Idents(onlymicroglia) <- onlymicroglia@meta.data$seurat_clusters
levels(onlymicroglia)

#subcluster microglia
onlymicroglia <-NormalizeData(onlymicroglia)
onlymicroglia <- FindVariableFeatures(onlymicroglia, selection.method = "vst", nfeatures = 400)
#scale data and run PCA
all.genes <- rownames(onlymicroglia)
onlymicroglia <- ScaleData(onlymicroglia, features = all.genes)
onlymicroglia <- RunPCA(onlymicroglia, features = VariableFeatures(object = onlymicroglia),npcs = 100)
onlymicroglia <- FindNeighbors(onlymicroglia, dims = 1:21)
onlymicroglia <- FindClusters(onlymicroglia, resolution = 2.4)
onlymicroglia <- RunUMAP(onlymicroglia, dims = 1:21)

#plot cell type markers to remove hybrids (cells expressing more than one cell-type marker due to segmentation erros)
markers <- c('Aqp4', 'Slc1a3', 'Aldh1l1', 'Cldn5', 'Pecam1', 'Slc17a6','Sst','Pvalb', 'Gad2', 'Clec7a','Cxcl3', 'Axl', 'P2ry12', 'Tmem119', 'Col4a1', 'Grm8', 'Pdgfrb', 'Sox2', 'Meg3', 'Nkg7', 'Plp1', 'Sox10', 'Pdgfra', 'Cspg4', 'Cd3e', 'Cd8a') 
DotPlot_scCustom(seurat_object = onlymicrogliaclean, features = markers, colors_use = viridis_plasma_dark_high)
#remove hybrids
onlymicrogliaclean <- subset(onlymicroglia,idents = c("0", "1", "2", "4", "5", "6","7",  "10",
                                                      "11",  "13", "14", "15",  "17","18", "19", '20',
                                                      "22", "23", "24","25", "26","27",  "28","29", "30", "31",  "33"))

Idents(onlymicrogliaclean) <- onlymicrogliaclean@meta.data$condition
levels(onlymicrogliaclean)
#re-cluster
onlymicrogliaclean <-NormalizeData(onlymicrogliaclean)
onlymicrogliaclean <- FindVariableFeatures(onlymicrogliaclean, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(onlymicrogliaclean)
onlymicrogliaclean <- ScaleData(onlymicrogliaclean, features = all.genes)
onlymicrogliaclean <- RunPCA(onlymicrogliaclean, features = VariableFeatures(object = onlymicrogliaclean),npcs = 100)
onlymicrogliaclean <- FindNeighbors(onlymicrogliaclean, dims = 1:20)
onlymicrogliaclean <- FindClusters(onlymicrogliaclean, resolution = 2.4)
onlymicrogliaclean <- RunUMAP(onlymicrogliaclean, dims = 1:20)
#Check celltype markers
DotPlot_scCustom(seurat_object = onlymicrogliaclean, features = markers, colors_use = viridis_plasma_dark_high) 
#remove hybrids
onlymicrogliaclean <- subset(onlymicrogliaclean,idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8",
                                                           "9", "10",  "11",  "13", "14","15",  
                                                           "17", "18", "19", "20", "21",  "23", "24", "25", "26", "29", "30", "31"))
#re-cluster
onlymicrogliaclean <-NormalizeData(onlymicrogliaclean)
onlymicrogliaclean <- FindVariableFeatures(onlymicrogliaclean, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(onlymicrogliaclean)
onlymicrogliaclean <- ScaleData(onlymicrogliaclean, features = all.genes)
onlymicrogliaclean <- RunPCA(onlymicrogliaclean, features = VariableFeatures(object = onlymicrogliaclean),npcs = 100)
onlymicrogliaclean <- FindNeighbors(onlymicrogliaclean, dims = 1:20)
onlymicrogliaclean <- FindClusters(onlymicrogliaclean, resolution = 0.4)
onlymicrogliaclean <- RunUMAP(onlymicrogliaclean, dims = 1:20)

#calculate plaques module score

plaques <- list(c('Cxcl3',
                  'Rbl1',
                  'Ccl20',
                  'Igfbp3',
                  'Mki67',
                  'Bmpr1a',
                  'Cst7'))


onlymicrogliaclean <- AddModuleScore(
  object = onlymicrogliaclean, 
  features = plaques, 
  ctrl = 5, 
  name = 'plaquesMg'
)

FeaturePlot_scCustom(onlymicrogliaclean, features = "plaquesMg1")

VlnPlot(onlymicrogliaclean, features = "plaquesMg1", pt.size = 0)

#remove plaques
onlymicrogliaclean.noplaques <- subset(onlymicrogliaclean, idents = c("0", "1", "2", "3", "4", "5", "7", "8"))
#re-cluster
onlymicrogliaclean.noplaques <-NormalizeData(onlymicrogliaclean.noplaques)
onlymicrogliaclean.noplaques <- FindVariableFeatures(onlymicrogliaclean.noplaques, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(onlymicrogliaclean.noplaques)
onlymicrogliaclean.noplaques <- ScaleData(onlymicrogliaclean.noplaques, features = all.genes)
onlymicrogliaclean.noplaques <- RunPCA(onlymicrogliaclean.noplaques, features = VariableFeatures(object = onlymicrogliaclean.noplaques),npcs = 100)
onlymicrogliaclean.noplaques <- FindNeighbors(onlymicrogliaclean.noplaques, dims = 1:20)
onlymicrogliaclean.noplaques <- FindClusters(onlymicrogliaclean.noplaques, resolution = 0.45)
onlymicrogliaclean.noplaques <- RunUMAP(onlymicrogliaclean.noplaques, dims = 1:20)

DimPlot(onlymicrogliaclean.noplaques, reduction = "umap",label=TRUE)
DimPlot(onlymicrogliaclean.noplaques, reduction = "umap",label=TRUE, cols = palette)
DimPlot(onlymicrogliaclean.noplaques, reduction = "umap",label=TRUE, split.by = "condition", pt.size = .3)


#apply a nFeature cutoff
onlymicrogliaclean.noplaquesfeaturecutoff <- subset(onlymicrogliaclean.noplaques, nFeature_RNA > 60)
onlymicrogliaclean.noplaquesfeaturecutoff <-NormalizeData(onlymicrogliaclean.noplaquesfeaturecutoff)
#re-cluster
onlymicrogliaclean.noplaquesfeaturecutoff <- FindVariableFeatures(onlymicrogliaclean.noplaquesfeaturecutoff, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(onlymicrogliaclean.noplaquesfeaturecutoff)
onlymicrogliaclean.noplaquesfeaturecutoff <- ScaleData(onlymicrogliaclean.noplaquesfeaturecutoff, features = all.genes)
onlymicrogliaclean.noplaquesfeaturecutoff <- RunPCA(onlymicrogliaclean.noplaquesfeaturecutoff, features = VariableFeatures(object = onlymicrogliaclean.noplaquesfeaturecutoff),npcs = 100)
onlymicrogliaclean.noplaquesfeaturecutoff <- FindNeighbors(onlymicrogliaclean.noplaquesfeaturecutoff, dims = 1:20)
onlymicrogliaclean.noplaquesfeaturecutoff <- FindClusters(onlymicrogliaclean.noplaquesfeaturecutoff, resolution = 0.4)
onlymicrogliaclean.noplaquesfeaturecutoff <- RunUMAP(onlymicrogliaclean.noplaquesfeaturecutoff, dims = 1:20)

#calculate senescence and DAM module scores

senescence <- list(c("Axl","C3","Ccl20","Ccl3","Ccl8","Cd9","Csf1","Csf2",
                     "Ctsb","Cxcl1","Cxcl10","Cxcl2","Cxcl3","Dkk1","Edn1",
                     "Esm1","Fgf1","Fgf2","Fgf7","Gdf15","Hgf","Hmgb1","Igfbp2",
                     "Igfbp3","Igfbp5","Il10","Il13","Il15","Il1a","Il1b","Il6",
                     "Il7","Lcp1","Mmp12","Mmp13","Mmp2","Mmp3","Mmp9","Pecam1",
                     "Sema3f","Serpine1","Spp1","Spx","Timp2","Tnf","Vegfc","Vgf",
                     "Wnt16"
))



DAM_features<- list(c('Abca1','Ap1b1','Apoe','Asah1','Axl','B2m','Ccl6','Cd63','Cd74',
                      'Cd9','Cebpa','Clec7a','Csf1','Cst7','Ctsb','Ctsd','Ctsl','Cxcl10',
                      'H2.Ea.ps','Hgf','Itgav','Itgax','Lipa','Lpl','Ly86','Lyz2','Mertk',
                      'Nupr1','Pla2g7','Pld3','Pparg','Ppt1','Spp1','Timp2','Trem2','Tyrobp'
)) 


onlymicrogliaclean.noplaquesfeaturecutoff <- AddModuleScore(
  object = onlymicrogliaclean.noplaquesfeaturecutoff, 
  features = senescence, 
  ctrl = 5, 
  name = 'senescenceMg_pure'
)



onlymicrogliaclean.noplaquesfeaturecutoff <- AddModuleScore(
  object = onlymicrogliaclean.noplaquesfeaturecutoff, 
  features = DAM_features, 
  ctrl = 5, 
  name = 'DAM_featuresMg_pure'
) 

#set a threshold

young <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (condition == "Young"))
mean(young$senescenceMg_pure1)
sd(young$senescenceMg_pure1)
threesd <- (sd(young$senescenceMg_pure1) * 3)
threesd
cutoff <- mean(young$senescenceMg_pure1)+threesd
cutoff

colour_breaks <- c( 0.2580792, 3)
colours <- c("yellow", "darkblue")

onlymicrogliaclean.noplaquesfeaturecutoff$senescencemg <- "non senescent"
onlymicrogliaclean.noplaquesfeaturecutoff$senescencemg[WhichCells(onlymicrogliaclean.noplaquesfeaturecutoff, expression =  senescenceMg_pure1 > 0.2580792)] <- "senescent"
Idents(onlymicrogliaclean.noplaquesfeaturecutoff) <- "senescence"
levels(onlymicrogliaclean.noplaquesfeaturecutoff)

# Senescent module score feature plots (repeat for each sample)
FeaturePlot_scCustom(seurat_object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = condition =='5xFAD') ,features = "senescenceMg_pure1", label = FALSE, pt.size = 1, split.by = "condition")+
  scale_colour_gradientn(
    limits  = range(onlymicrogliaclean.noplaquesfeaturecutoff$senescenceMg_pure1),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = range(onlymicrogliaclean.noplaquesfeaturecutoff$senescenceMg_pure1)), 1),
  )


onlymicrogliaclean.noplaquesfeaturecutoff$senescencemg <- "non senescent"
onlymicrogliaclean.noplaquesfeaturecutoff$senescencemg[WhichCells(onlymicrogliaclean.noplaquesfeaturecutoff, expression =  senescenceMg_pure1 > 0.2580792)] <- "senescent"
Idents(onlymicrogliaclean.noplaquesfeaturecutoff) <- "senescence"
levels(onlymicrogliaclean.noplaquesfeaturecutoff)

no <- subset(x=onlymicrogliaclean.noplaquesfeaturecutoff, subset = senescenceMg_pure1 < 0.2580792)
yes <- subset(x=onlymicrogliaclean.noplaquesfeaturecutoff, subset = senescenceMg_pure1 > 0.2580792)
onlymicrogliaclean.noplaquesfeaturecutoff_plot <- merge(no, y = c(yes))
levels(onlymicrogliaclean.noplaquesfeaturecutoff_plot)


Idents(onlymicrogliaclean.noplaquesfeaturecutoff_plot)<-onlymicrogliaclean.noplaquesfeaturecutoff_plot$senescencemg
levels(onlymicrogliaclean.noplaquesfeaturecutoff_plot)
table(Idents(onlymicrogliaclean.noplaquesfeaturecutoff_plot))

#spatial plots (repeat for each sample)
dfall<- data.frame(onlymicrogliaclean.noplaquesfeaturecutoff_plot$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot, sample =="Young_2"))],onlymicrogliaclean.noplaquesfeaturecutoff_plot$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))],onlymicrogliaclean.noplaquesfeaturecutoff_plot@active.ident[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))])
x <- onlymicrogliaclean.noplaquesfeaturecutoff_plot$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))]
y <- onlymicrogliaclean.noplaquesfeaturecutoff_plot$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))]
group <- onlymicrogliaclean.noplaquesfeaturecutoff_plot@active.ident[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))]



p <- ggplot(dfall ,aes(x,y,colour=group))+geom_point(size=1)+  theme_classic()
p + scale_color_manual(values=c( "lightgray","darkblue")) & NoLegend()



# Set a threshold for DAM
young <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (condition == "Young"))
mean(young$DAM_featuresMg_pure1)
sd(young$DAM_featuresMg_pure1)
threesd <- (sd(young$DAM_featuresMg_pure1) * 3)
threesd
cutoff <- mean(young$DAM_featuresMg_pure1)+threesd
cutoff


colour_breaks <- c(0.7794895, 3)
colours <- c("yellow", "red")

# DAM module score feature plots (repeat for each sample)
FeaturePlot_scCustom(seurat_object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = condition =='Young') ,features = "DAM_featuresMg_pure1", label = FALSE, pt.size = 1, split.by = "condition")+
  scale_colour_gradientn(
    limits  = range(onlymicrogliaclean.noplaquesfeaturecutoff$DAM_featuresMg_pure1),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = range(onlymicrogliaclean.noplaquesfeaturecutoff$DAM_featuresMg_pure1)), 1),
  )
onlymicrogliaclean.noplaquesfeaturecutoff$DAMmg <- "non DAM"
onlymicrogliaclean.noplaquesfeaturecutoff$DAMmg[WhichCells(onlymicrogliaclean.noplaquesfeaturecutoff, expression =  DAM_featuresMg_pure1 > 0.7794895)] <- "DAM"
Idents(onlymicrogliaclean.noplaquesfeaturecutoff) <- "DAMmg"
levels(onlymicrogliaclean.noplaquesfeaturecutoff)

no <- subset(x=onlymicrogliaclean.noplaquesfeaturecutoff, subset = DAM_featuresMg_pure1 < 0.7794895)
yes <- subset(x=onlymicrogliaclean.noplaquesfeaturecutoff, subset = DAM_featuresMg_pure1 > 0.7794895)
onlymicrogliaclean.noplaquesfeaturecutoff_plot <- merge(no, y = c(yes))
levels(onlymicrogliaclean.noplaquesfeaturecutoff_plot)


Idents(onlymicrogliaclean.noplaquesfeaturecutoff_plot)<-onlymicrogliaclean.noplaquesfeaturecutoff_plot$DAMmg
levels(onlymicrogliaclean.noplaquesfeaturecutoff_plot)
table(Idents(onlymicrogliaclean.noplaquesfeaturecutoff_plot))

#spatial plots (repeat for each sample)
dfall<- data.frame(onlymicrogliaclean.noplaquesfeaturecutoff_plot$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot, sample =="Young_2"))],onlymicrogliaclean.noplaquesfeaturecutoff_plot$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))],onlymicrogliaclean.noplaquesfeaturecutoff_plot@active.ident[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))])
x <- onlymicrogliaclean.noplaquesfeaturecutoff_plot$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))]
y <- onlymicrogliaclean.noplaquesfeaturecutoff_plot$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))]
group <- onlymicrogliaclean.noplaquesfeaturecutoff_plot@active.ident[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff_plot,subset = sample =="Young_2"))]



p <- ggplot(dfall ,aes(x,y,colour=group))+geom_point(size=1)+  theme_classic()
p + scale_color_manual(values=c( "lightgray","red")) & NoLegend()


#Find markers for each Mg cluster
Idents(onlymicrogliaclean.noplaquesfeaturecutoff) <- onlymicrogliaclean.noplaquesfeaturecutoff$seurat_clusters

markers.microglia_noplaques <- FindAllMarkers(onlymicrogliaclean.noplaquesfeaturecutoff, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.microglia_noplaques %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(onlymicrogliaclean.noplaquesfeaturecutoff, features = top5$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

#spatial plot for each Mg cluster
palette <- c("#D55E00","#E69F00" ,"#56B4E9", "#009E73", "#F0E442", "#0072B2",   "#666666","#CC79A7", "#AD7700" ,"#1C91D4" ,"#007756",
             "#06A5FF", "#D5C711", "#005685", "#A04700" ,"#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF" ,"#00F6B3", "#F4EB71" )


for (i in 0) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 10000)+ylim(0,8000)+  
          geom_point(color="#D55E00")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

for (i in 1) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 10000)+ylim(0,8000)+ 
          geom_point(color="#E69F00")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

for (i in 2) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 10000)+ylim(0,8000)+  
          geom_point(color="#56B4E9")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

for (i in 3) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 10000)+ylim(0,8000)+ 
          geom_point(color="#009E73")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

for (i in 4) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 10000)+ylim(0,8000)+
          geom_point(color="#F0E442")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}


for (i in 5) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 10000)+ylim(0,8000)+  
          geom_point(color="#0072B2")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}


for (i in 6) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 10000)+ylim(0,8000)+
          geom_point(color="#666666")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}


for (i in 7) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 10000)+ylim(0,8000)+
          geom_point(color="#CC79A7")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}


for (i in 8) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 8000)+ylim(0,11000)+  
          geom_point(color="#AD7700")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

for (i in 9) {
  cluster <- i
  sampleA <- "EAE_1"
  x0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.x[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  y0 <- onlymicrogliaclean.noplaquesfeaturecutoff$center.y[WhichCells(object = subset(onlymicrogliaclean.noplaquesfeaturecutoff,subset = sample ==sampleA), ident = cluster)]
  df<-data.frame(x0,y0)
  print(ggplot(df,aes(x0,y0))+geom_point(size=1)+xlim(0, 8000)+ylim(0,11000)+   
          geom_point(color="#1C91D4")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

#plot DAM and senescent module score
Split_FeatureScatter(seurat_object = onlymicrogliaclean.noplaquesfeaturecutoff, feature1 = "senescenceMg_pure1", feature2 = "DAM_featuresMg_pure1",split.by="condition", group.by = "ident", num_columns = 2,
                     pt.size = 1, colors_use = palette)+geom_hline(yintercept=0.7790245, size=1)+geom_vline(xintercept = 0.3596297, size=1)

#coordinates for microglia spatial analysis in different brain regions

Young1 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "Young_1"))
Young2 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "Young_2"))
Young3 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "Young_3"))
Aged1 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "Aged_1"))
Aged2 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "Aged_2"))
Aged3 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "Aged_3"))
WT1 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "WT-5xFAD_1"))
WT2 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "WT-5xFAD_2"))
WT3 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "WT-5xFAD_3"))
AD1 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "5xFAD_1"))
AD2 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "5xFAD_2"))
AD3 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "5xFAD_3"))
CFA1 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "CFA_1"))
CFA2 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "CFA_2"))
CFA3 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "CFA_3"))
EAE1 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "EAE_1"))
EAE2 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "EAE_2"))
EAE3 <- subset(x = onlymicrogliaclean.noplaquesfeaturecutoff, subset = (sample == "EAE_3"))

#Young_1

Young1.x <- Young1@meta.data[,5]
head(Young1.x)
Young1.y <- Young1@meta.data[,6]
head(Young1.y)

Young1_cx_x <- c(1477.01,	924.773,	767.014,	682.024,	688.718,	2148.21,	3331.74,	5626.46,	7568.85,	7757.95,	7402.14,	6357.34,	5507.55,	4612.61,	4762.29,	5404.77,	5990.72,	5372.95,	4937.37,	4108.96,	3225.65,	2519.82,	1680.93,	1917.17,	1852,	1288.09,	1196.99,	1245.33,	1517.51,	2053.2,	2937.18,	4142.18,	3673.16,	2882.53,	3198.68,	4043.78,	5203.37,	6115.1,	6033.79,	5475.23,	2753.69,	1389.63,	416.063,	311.687,	472.63)
Young1_cx_y <- c(5958.26,	6121.62,	6276.85,	6519.98,	6679.54,	9091.87,	9959.81,	10150,	8673.16,	7175.9,	6501.29,	6297.89,	7860.87,	8470,	8763.5,	8678.37,	8442.21,	8790.47,	8934.41,	8833.17,	8501.39,	8039.56,	6863.92,	5871.24,	5615.19,	5099.48,	4027.05,	3539.11,	2827.09,	2018.76,	1364.37,	1145.1,	1296.46,	1549.13,	1922.63,	2168.92,	2821.87,	2300.05,	729.053,	268.878,	309.189,	746.072,	2179.56,	3821.91,	6081.34)
Young1_hp_x <- c(2579.2,	2781.03,	2804.07,	2387.7,	2246.68,	2205.82,	2096.05,	1867.37,	1448.96,	1368.59,	1442.69,	1847.79,	2311.85,	2351.31,	2162.43,	2053.62,	2042.11,	2318.61,	2842.59,	3451.22,	3700.64,	3461.7,	3243.22,	3140.97,	3107.13,	2766.94)
Young1_hp_y <- c(5329.81,	5036.12,	4849.13,	3950.72,	3430.56,	2817.86,	2700.26,	2761.73,	3624.1,	4321.8,	4659.88,	5138.15,	5549.99,	5838.21,	6422.58,	6748.32,	7009.35,	7506.47,	8005.79,	8300.39,	8113.92,	7676.94,	7109.2,	6654.7,	6066.71,	5782.68)
Young1_wm_x <- c(1930.57,	1686.44,	2365.28,	2872.62,	3907.21,	4523.03,	4961.44,	5981.85,	5103.7,	5079.74,	5544.56,	5671.9,	5132.56,	4176.3,	3713.94,	3720.18,	3451.22,	2779.1,	2052.88,	2053.62,	2294.75,	2277.01,	1778.98,	1482.53,	1338.81,	1541.13,	1950.06,	2193.2,	2285.43,	2966.51,	4024.28,	3675.29,	3293.46,	2606.42,	2282.63,	2793.83,	3698.05,	2984.63,	2354.12,	1711.25,	1414.2,	1216.76,	1226.03,	1424.63,	1888.96)
Young1_wm_y <- c(5863.11,	6683.25,	7845.45,	8266.66,	8736.93,	8869.82,	8894.6,	8480.5,	8769.42,	8268.7,	7694.88, 7387.1,	7719.91,	7815.14,	7839.44,	8211.1,	8300.39,	7978.76,	7116.99,	6748.32,	6053.52,	5540.45,	5094.51,	4786.9,	4159.8,	3297.17,	2700.79,	2744.91,	3125.88,	2616.14,	2265.1,	2083.44,	1979.77,	2015.27,	2006.39,	1553.1,	1275.94,	1350.35,	1725.24,	2461.08,	3101.21,	3689.81,	4438.73,	5274.18,	5729.44)
Young1_th_x <- c(3275.84,	3133.83,	3306.05,	3662.63,	4640.67,	5422.69,	5893.89,	7152.98,	6799.51,	5221.21,	5058.06,	4603.73,	3411.91,	2403.86,	2287.82,	2617.1,	3184.86)
Young1_th_y <- c(5410.22,	6411.4,	7177.91,	7769.55,	7718.15,	6720.75,	6000.78,	5744.52,	3825.02,	3258.56,	3455.03,	3155.36,	2532.47,	3005.77,	3362.35,	4498.28,	5344.85)

Young1_pol_cx <- point.in.polygon(Young1.x, Young1.y, Young1_cx_x, Young1_cx_y)
Young1_pol_hp <- point.in.polygon(Young1.x, Young1.y, Young1_hp_x, Young1_hp_y)
Young1_pol_wm <- point.in.polygon(Young1.x, Young1.y, Young1_wm_x, Young1_wm_y)
Young1_pol_th <- point.in.polygon(Young1.x, Young1.y, Young1_th_x, Young1_th_y)

Young1.info <- cbind.data.frame(Young1_pol_cx,Young1_pol_hp,Young1_pol_wm,Young1_pol_th)
head(Young1.info, 50)


Young1.info <- data.frame(Young1.info)
Young1@meta.data <- cbind(Young1@meta.data,Young1.info)


Young1_cx <- subset(Young1,subset = Young1_pol_cx == 1)
Young1_hp <- subset(Young1,subset = Young1_pol_hp == 1)
Young1_wm <- subset(Young1,subset = Young1_pol_wm == 1)
Young1_th <- subset(Young1,subset = Young1_pol_th == 1)

Young1_cx$region <- "Cx"
Young1_hp$region <- "Hp"
Young1_wm$region <- "Wm"
Young1_th$region <- "Th"

Young1_all <- merge(Young1_cx, y = c( Young1_hp, Young1_wm, Young1_th))


Idents(Young1_all)<-Young1_all$region
levels(Young1_all)

dfall<- data.frame(Young1_all$center.x[WhichCells(object = subset(Young1_all,subset = sample =="Young_1"))],Young1_all$center.y[WhichCells(object = subset(Young1_all,subset = sample =="Young_1"))],Young1_all@active.ident[WhichCells(object = subset(Young1_all,subset = sample =="Young_1"))])
x <- Young1_all$center.x[WhichCells(object = subset(Young1_all,subset = sample =="Young_1"))]
y <- Young1_all$center.y[WhichCells(object = subset(Young1_all,subset = sample =="Young_1"))]
group <- Young1_all@active.ident[WhichCells(object = subset(Young1_all,subset = sample =="Young_1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()



#Young_2

Young2.x <- Young2@meta.data[,5]
head(Young2.x)
Young2.y <- Young2@meta.data[,6]
head(Young2.y)

Young2_cx_x <- c(6070.73,	9994.72,	9972.51,	7035.38,	5536.75,	6403.13,	7125.53,	7571.51,	8165.36,	7851.51,	8446.6,	8511.72,	8348.66,	8000.29,	7362.57,	6818.28,	6346,	5646.36,	5215.02,	4775.32,	4137.3,	3405.19,	2418.04,	1761.69,	1288.22,	1143.74,	1195.58,	1267.16,	1902.66,	2026.19,	2808.59,	2704.35,	852.021,	-24.9523,	441.378,	2813.4,	5595.17)
Young2_cx_y <- c(7584.35,	4969.54,	2604.65,	-89.8634,	950.259,	2115.81,	2358.34,	2714.63,	2690.48,	2168.94,	2821.69,	3517.25,	4136.63,	4804.57,	5545.42,	5901.67,	6089.75,	5890.98,	6009.88,	6485.5,	6587.21,	6575.71,	6265.72,	5774.94,	5136.09,	4415.49,	3938.93,	4625.06,	4598.14,	3795.81,	2384.24,	1721.86,	1943.74,	4241.78,	6589.94,	8167.88,	7800.35)
Young2_hp_x <- c(5382.08,	5739.61,	6287.18,	6723.84,	7340.03,	7432.24,	7334.06,	7014.31,	6681.52,	6409.86,	5784.6,	5620.63,	5386.27,	4998.9,	4652.19,	4169.43,	3978.14,	3353.15,	3284.41,	3820.09,	4374.39,	4782.59,	4946.59,	5299.5)
Young2_hp_y <- c(5584.2,	5643.44,	5808.76,	5746.82,	5283.26,	5015.5,	4905.95,	4842.89,	4882.18,	4897.34,	4841.02,	4889.72,	5227.14,	5306.3,	5179.67,	5460.68,	5652.68,	5939.56,	6115.15,	6355.13,	6334.04,	6060.62,	5862.54,	5616.19)
Young2_wm_x <- c(5601.92,	6235.86,	7096.91,	7844.51,	8378.9,	8479.27,	8501.47,	7907.44,	7996.47,	7707.23,	7804.49,	8219.36,	8428.77,	8336.78,	8018.64,	7713.84,	7499.85,	7193.85,	6905.44,	6590.68,	7140.86,	7443.04,	7228.22,	6657.41,	6117.81,	5585.4,	5359.39,	4906.31,	4494.89,	4055.28,	3661.23,	3360.7,	3304.66,	3601.24,	4044.64,	3682.88,	3388.13,	2906.8,	2190.18,	2383.73,	2164.42,	1814.79,	1387.8,	1198.01,	1177.07,	1397.95,	1852.88,	2514.65,	3170,	4747.67,	5238.83)
Young2_wm_y <- c(5854.69,	6072.91,	5755.06,	5044.8,	4049.1,	3528.7,	3053.96,	2673.39,	2245.76,	1656.58,	2089.03,	2699.33,	3246.83,	3918.39,	4446.02,	3617.44,	2916.62,	4225.9,	4690.04,	4874.3,	4847.69,	4996.82,	5377.11,	5793.14,	5777.81,	5628.36,	5600.36,	5906.01,	6247.38,	6358.86,	6352.7,	6210.45,	6060.07,	5822.22,	5564.84,	5663.62,	5618.19,	5283.04,	5106.51,	5760.28,	5899.4,	5664.51,	5131.39,	4453.23,	4829.51,	5374.78,	5859.1,	6298.83,	6520.1,	6500.62,	5959.63)
Young2_th_x <- c(5163.95,	5907.7,	6458.65,	6875.58,	7213.01,	6935.91,	5351.81,	4939.62,	3398.64,	2507.73,	2551.5,	3070.9,	3529.69,	3895.43,	4080.98,	4985.51)
Young2_th_y <- c(4807.59,	4838.48,	4852.26,	4684.91,	3834.23,	3298.37,	1576.21,	1152.22,	1579.64,	3794.24,	4774.22,	5346.39,	5584.34,	5586.52,	5509.22,	4855.76)


Young2_pol_cx <- point.in.polygon(Young2.x, Young2.y, Young2_cx_x, Young2_cx_y)
Young2_pol_hp <- point.in.polygon(Young2.x, Young2.y, Young2_hp_x, Young2_hp_y)
Young2_pol_wm <- point.in.polygon(Young2.x, Young2.y, Young2_wm_x, Young2_wm_y)
Young2_pol_th <- point.in.polygon(Young2.x, Young2.y, Young2_th_x, Young2_th_y)

Young2.info <- cbind.data.frame(Young2_pol_cx,Young2_pol_hp,Young2_pol_wm,Young2_pol_th)
head(Young2.info, 50)



Young2.info <- data.frame(Young2.info)
Young2@meta.data <- cbind(Young2@meta.data,Young2.info)


Young2_cx <- subset(Young2,subset = Young2_pol_cx == 1)
Young2_hp <- subset(Young2,subset = Young2_pol_hp == 1)
Young2_wm <- subset(Young2,subset = Young2_pol_wm == 1)
Young2_th <- subset(Young2,subset = Young2_pol_th == 1)

Young2_cx$region <- "Cx"
Young2_hp$region <- "Hp"
Young2_wm$region <- "Wm"
Young2_th$region <- "Th"

Young2_all <- merge(Young2_cx, y = c( Young2_hp, Young2_wm, Young2_th))


Idents(Young2_all)<-Young2_all$region
levels(Young2_all)

dfall<- data.frame(Young2_all$center.x[WhichCells(object = subset(Young2_all,subset = sample =="Young_2"))],Young2_all$center.y[WhichCells(object = subset(Young2_all,subset = sample =="Young_2"))],Young2_all@active.ident[WhichCells(object = subset(Young2_all,subset = sample =="Young_2"))])
x <- Young2_all$center.x[WhichCells(object = subset(Young2_all,subset = sample =="Young_2"))]
y <- Young2_all$center.y[WhichCells(object = subset(Young2_all,subset = sample =="Young_2"))]
group <- Young2_all@active.ident[WhichCells(object = subset(Young2_all,subset = sample =="Young_2"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()


Young3.x <- Young3@meta.data[,5]
head(Young3.x)
Young3.y <- Young3@meta.data[,6]
head(Young3.y)

Young3_cx_x <- c(1659.66,	1268.81,	1402.89,	1639.06,	2240.82,	2669.63,	3726.07,	5093.51,	5685.02,	4619.87,	6357.57,	7271.59,	8047.93,	6283.66,	4100.32,	2520.18,	445.14,	287.129,	289.058,	411.459,	1821.22,	3706.59,	6296.23,	7513.27,	7050.3,	6090.74,	5709.86,	5291.79,	4739.22,	3946.52,	2923.65,	1959.25,	1596.21,	1350.89,	1182.29,	1728.44)
Young3_cx_y <- c(54953.78,	4759.81,	4071.82,	3498.16,	2639.3,	2202.91,	1709.21,	1732.38,	1969.35,	1953.64,	4336.17,	3756.98,	2338.27,	396.733,	178.153,	845.807,	4039.66,	4622.6,	6347.47,	8401.16,	10592,	10985.7,	10834.9,	9092.35,	7359.04,	8130.97,	9264.23,	9664.75,	9949.77,	10008.5,	9633.12,	8721.8,	8136.83,	7677.17,	6514.31,	5774.3)
Young3_hp_x <- c(2261.13,	1969.9,	1558.52,	1523.62,	1823.38,	2287.87,	2897.27,	3737.91,	4195.25,	3422.88,	2908.1,	2812.61,	2662.33,	2634.5,	2598.18,	2523.55,	2791.74,	3035.53,	3466.66,	3984.48,	4779.1,	5556.38,	6008.14,	5826.02,	5610.25,	5123.97,	4749.13,	4090.05,	3554.34,	2933.76,	2420.64,	2035.6,	1708.66,	1505.04,	1553.77,	2058.52,	2209.29,	2231.4,	1925.25,	1572.53,	1535.66,	1705.16,	2106.58,	2588.99,	3257.67,	3990.98,	4323.79,	3771.93,	3299.06,	3049.37,	2857.79,	2662.33,	2668.47,	2639.22,	2681.18,	2857.75,	2996.39,	3662.36,	4323.79,	3680.33,	3051.52,	2506.89,	1819.35,	1557.75,	1520,	1727.41,	1938.05,	2239.55)
Young3_hp_y <- c(5579.09,	5070.86,	4580.03,	4190.24,	3512.49,	2889.33,	2446.93,	2191.05,	2328.97,	3008.87,	3617.17,	4305.72,	4613.37,	5398.35,	6071.81,	6463.85,	7656.15,	8177.82,	8494.11,	8642.17,	8407.91,	7959.45,	7805.76,	8591.49,	9126.82,	9586.18,	9763.02,	9808.95,	9638.34,	9278.69,	8899.68,	8442.01,	7981.97,	7335.23,	6883.14,	6434.04,	6196.38,	5549.05,	4997.65,	4669.47,	4356.93,	3659.86,	3073.03,	2625.85,	2287.15,	2208.89,	2301.02,	2778.1,	3133.56,	3412.28,	3863.8,	4613.37,	5448.1,	5468.72,	4736.45,	3901.34,	3420.25,	2869.07,	2301.02,	2255,	2350.6,	2703.47,	3459.47,	4153.82,	4385.68,	4957.87,	4985.67,	5593.01)
Young3_wm_x <- c(1697.29,	1188.27,	1410.6,	1594.35,	1871.16,	2316.83,	2878.19,	3854.59,	4689.01,	5405.65,	4455.56,	5294.82,	5324.94,	4764.76,	4227.91,	3745.37,	4440.49,	4323,	3372.55,	2581.29,	2081.26,	1666.53,	1503.89,	1755.97,	1941.56,	2192.28,	2209.72,	2194,	1961.19,	1589.56,	1480.33,	1576.28,	2112.87,	2599.16,	3211.61,	4020.94,	4558.63,	5142.63,	5518.81,	5804.69,	5803.28,	5643.69,	5317.21,	4678.29,	3978.32,	2975.51,	2060.05,	1427.68,	1318.42,	1189.81,	1632.7,	1778.41)
Young3_wm_y <- c(5523.87,	4675.53,	4060.05,	3640.09,	3162.06,	2595.07,	2114.74,	1725.28,	1703.81,	1863.34,	1918.67,	2746.4,	3282.66,	3047.93,	2915.45,	2883.84,	2358.47,	2163.69,	2220.96,	2606.32,	3097.47,	3705.18,	4348.82,	5060.25,	5073.49,	5472.96,	5858.08,	6252.79,	6464.68,	6769.87,	7140.23,	7722.89,	8599.69,	9098.69,	9492.89,	9831.51,	9873.57,	9608.22,	9324.96,	8697.94,	8917.48,	9300.88,	9595.11,	9908.33,	9971.9,	9548.99,	8852.4,	7686.96,	7277.18,	6401.93,	6001.22,	5745.46)
Young3_th_x <- c(2706.18,	2708.68,	2887.88,	3317.51,	4000.62,	4624.38,	4975.35,	5906.27,	7123.41,	6958.8,	5954.73,	4984.04,	4552.2,	3995.36,	3278.47,	2988,	2707.75,	2599.99)
Young3_th_y <- c(5505.69,	4589.35,	3864.42,	3223.19,	2959.53,	3118.04,	3713.28,	4670.21,	4752.15,	7016.03,	6933.14,	7755.96,	8368.58,	8553.01,	8252.73,	7925.17,	6940.1,	6192.99)


Young3_pol_cx <- point.in.polygon(Young3.x, Young3.y, Young3_cx_x, Young3_cx_y)
Young3_pol_hp <- point.in.polygon(Young3.x, Young3.y, Young3_hp_x, Young3_hp_y)
Young3_pol_wm <- point.in.polygon(Young3.x, Young3.y, Young3_wm_x, Young3_wm_y)
Young3_pol_th <- point.in.polygon(Young3.x, Young3.y, Young3_th_x, Young3_th_y)

Young3.info <- cbind.data.frame(Young3_pol_cx,Young3_pol_hp,Young3_pol_wm,Young3_pol_th)
head(Young3.info, 50)



Young3.info <- data.frame(Young3.info)
Young3@meta.data <- cbind(Young3@meta.data,Young3.info)


Young3_cx <- subset(Young3,subset = Young3_pol_cx == 1)
Young3_hp <- subset(Young3,subset = Young3_pol_hp == 1)
Young3_wm <- subset(Young3,subset = Young3_pol_wm == 1)
Young3_th <- subset(Young3,subset = Young3_pol_th == 1)

Young3_cx$region <- "Cx"
Young3_hp$region <- "Hp"
Young3_wm$region <- "Wm"
Young3_th$region <- "Th"

Young3_all <- merge(Young3_cx, y = c( Young3_hp, Young3_wm, Young3_th))


Idents(Young3_all)<-Young3_all$region
levels(Young3_all)

dfall<- data.frame(Young3_all$center.x[WhichCells(object = subset(Young3_all,subset = sample =="Young_3"))],Young3_all$center.y[WhichCells(object = subset(Young3_all,subset = sample =="Young_3"))],Young3_all@active.ident[WhichCells(object = subset(Young3_all,subset = sample =="Young_3"))])
x <- Young3_all$center.x[WhichCells(object = subset(Young3_all,subset = sample =="Young_3"))]
y <- Young3_all$center.y[WhichCells(object = subset(Young3_all,subset = sample =="Young_3"))]
group <- Young3_all@active.ident[WhichCells(object = subset(Young3_all,subset = sample =="Young_3"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#Aged_1

Aged1.x <- Aged1@meta.data[,5]
head(Aged1.x)
Aged1.y <- Aged1@meta.data[,6]
head(Aged1.y)

Aged1_cx_x <- c(5932.61,	6516.76,	7324.24,	8021.71,	8576.96,	9001.85,	9304.4,	8476.33,	9443.66,	9377.21,	8946.04,	8243.92,	7768.8,	10692.9,	10671.4,	7272.77,	2879.78,	851.196,	-22.4738,	205.035,	1462.93,	1419.77,	976.2,	898.366,	917.743,	1021.72,	1388.4,	2285.88,	3384.53,	4188.04,	5081.4,	5587.85)
Aged1_cx_y <- c(2263.19,	1984.3,	2332.85,	2729.92,	3336.13,	3863.31,	4522.97,	5420.88,	5761.83,	6111.13,	6795.5,	7411.77,	8385.46,	7193.08,	3330,	721.004,	152.991,	1849.12,	3714,	5774.25,	6500.6,	5796.75,	5040.45,	4621.96,	4076.2,	3463.29,	2888.98,	2143.27,	1706.87,	1578.88,	1585.31,	2156.53)
Aged1_hp_x <- c(5637.24,	5117.56,	4608.15,	3363.04,	2704.27,	2338.71,	1256.97,	1073.51,	1124.55,	1982.75,	2747.58,	2696.57,	2451.81,	2492.52,	2799.24,	3211.18,	3595.91,	4161.71,	4977.67,	5379.98,	5878.57,	6123.11,	7152.86,	7704.78,	8022.2,	8175.38,	8048.64,	7601.39,	6982.99,	7593,	8748.11,	9038.83,	9253.48,	9257.13,	9148.6,	8746.16,	8315.21,	7671.01,	7092.97,	6850.08,	6322.66,	5840.8)
Aged1_hp_y <- c(2661.65,	2137.04,	1829.26,	2031.62,	2529.92,	2875.9,	3734.22,	4202.49,	4888.66,	6265.79,	5643.65,	4715.42,	4060.5,	3636.66,	3027.17,	2823.95,	2793.23,	2989.95,	3350.79,	3213.16,	3398.4,	3618.03,	3785.09,	3966.03,	4422.9,	4920.17,	5346.51,	5845.41,	6791.82,	7568,	6837.75,	6454.95,	6004.17,	5583.35,	5128.76,	4401.11,	3473.01,	2811.99,	2505.62,	2418.23,	2538.56,	2713.48)
Aged1_wm_x <- c(5818.82,	6514.51,	7073.94,	7743.72,	3857.33,	8943.2,	9310,	9465.24,	9410.78,	9105.25,	8803.07,	8198.33,	8584.71,	9038.83,	9253.59,	9297.9,	9163.38,	8681.05,	8262.41,	7752.18,	7171.29,	6606.78,	5855.27,	5514.4,	4768.04,	4188.59,	3300.92,	2556.11,	1843.06,	1209.94,	1030.95,	1114.09,	1267.51,	1009.81,	929.74,	968.922,	1102.53,	1281.78,	2129.89,	2723.34,	3608.07,	4567.59,	5201.99,	5726.31)
Aged1_wm_y <- c(2198.32,	1957.34,	2240.25,	2595.37,	3072.5,	3778.97,	4669.21,	5442.55,	6006.4,	6557.21,	6983.96,	7401.13,	7048.78,	6454.95,	6042.01,	5722.8,	5136.31,	4110.79,	3353.25,	2856.34,	2517.15,	2425.44,	2691.72,	2573.35,	1876.75,	1801.36,	2052.29,	2657.49,	3154.47,	3719.45,	4196.36,	4956.54,	5347.76,	4911.91,	4383.03,	3968.63,	3476.67,	3051.83,	2356.96,	1934.14,	1666.28,	1525.78,	1722.62,	2260.46)
Aged1_th_x <- c(5665.59,	7093.65,	7702.76,	8002.54,	8027.96,	7754.66,	6876.43,	6401.91,	2898.12,	2729.7,	2497.25,	2703.92,	3065.13,	3393.54,	3923.16,	4657.29,	5413.05)
Aged1_th_y <- c(3597.82,	3788.81,	4021.2,	4578.39,	5066.57,	5578.82,	6698.67,	6847.5,	5925.95,	4467.55,	3932.75,	3472.75,	3029.24,	2883.26,	2928.83,	3252.74,	3572.93)


Aged1_pol_cx <- point.in.polygon(Aged1.x, Aged1.y, Aged1_cx_x, Aged1_cx_y)
Aged1_pol_hp <- point.in.polygon(Aged1.x, Aged1.y, Aged1_hp_x, Aged1_hp_y)
Aged1_pol_wm <- point.in.polygon(Aged1.x, Aged1.y, Aged1_wm_x, Aged1_wm_y)
Aged1_pol_th <- point.in.polygon(Aged1.x, Aged1.y, Aged1_th_x, Aged1_th_y)

Aged1.info <- cbind.data.frame(Aged1_pol_cx,Aged1_pol_hp,Aged1_pol_wm,Aged1_pol_th)
head(Aged1.info, 50)



Aged1.info <- data.frame(Aged1.info)
Aged1@meta.data <- cbind(Aged1@meta.data,Aged1.info)


Aged1_cx <- subset(Aged1,subset = Aged1_pol_cx == 1)
Aged1_hp <- subset(Aged1,subset = Aged1_pol_hp == 1)
Aged1_wm <- subset(Aged1,subset = Aged1_pol_wm == 1)
Aged1_th <- subset(Aged1,subset = Aged1_pol_th == 1)

Aged1_cx$region <- "Cx"
Aged1_hp$region <- "Hp"
Aged1_wm$region <- "Wm"
Aged1_th$region <- "Th"

Aged1_all <- merge(Aged1_cx, y = c( Aged1_hp, Aged1_wm, Aged1_th))


Idents(Aged1_all)<-Aged1_all$region
levels(Aged1_all)

dfall<- data.frame(Aged1_all$center.x[WhichCells(object = subset(Aged1_all,subset = sample =="Aged_1"))],Aged1_all$center.y[WhichCells(object = subset(Aged1_all,subset = sample =="Aged_1"))],Aged1_all@active.ident[WhichCells(object = subset(Aged1_all,subset = sample =="Aged_1"))])
x <- Aged1_all$center.x[WhichCells(object = subset(Aged1_all,subset = sample =="Aged_1"))]
y <- Aged1_all$center.y[WhichCells(object = subset(Aged1_all,subset = sample =="Aged_1"))]
group <- Aged1_all@active.ident[WhichCells(object = subset(Aged1_all,subset = sample =="Aged_1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#Aged_2

Aged2.x <- Aged2@meta.data[,5]
head(Aged2.x)
Aged2.y <- Aged2@meta.data[,6]
head(Aged2.y)

Aged2_cx_x <- c(-41.456,	48.4751,	1435.09,	3108.67,	6128349,	6825.78,	6340.85,	5373.31,	40963,	3864.09,	3081.9,	3405.89,	4903.45,	4168.43,	3167.18,	2749.23,	1928.76,	1541.83,	1262.36,	1080.05,	1070.28,	1084.53,	1471.5,	1526.91,	1340.35,	1552.94,	2226.41,	3195.64,	3860.9,	4529.39,	5470.54,	4762.78,	4694.56,	5548.35,	5891.4,	6549.95,	7365.03,	7243.73,	4026.8,	2878.16,	2204.4,	1475.75,	1054.59,	530.375)
Aged2_cx_y <- c(5368.52,	2561.28,	721.767,	89.1081,	763.077,	1814.15,	3278.58,	3032.04,	2028.56,	1908.44,	1533.57,	1411.15,	1297.91,	1236.59,	1288.09,	1462.13,	2034.55,	2548.82,	3068.25,	3602.56,	4180.05,	4920.12,	5503.43,	5773.26,	6674.66,	7405.6,	8340.47,	8893.68,	9091.17,	9093.93,	8743.4,	8924.47,	8416,	7419.46,	6983.34,	6655.88,	7813.72,	9103.47,	10306.9,	10201.6,	9792.5,	9041.83,	8186.58,	6759.85)
Aged2_hp_x <- c(2042.07,	1742.93,	1602.24,	1674.62,	2266.86,	2721.12,	2882.52,	2871.72,	2772,	2870.42,	2826.2,	2523.41,	2464.15,	2670.43,	2369.45,	2164.5,	1913.85,	1509.71,	1244.79,	1177.52,	1285.41,	1924.14)
Aged2_hp_y <- c(5493.69,	6222.66,	6808.93,	7231.28,	8040.7,	8161.46,	7817.69,	7363.66,	6985.44,	6132.35,	5959.33,	5664.48,	5187.95,	4744.87,	3982.26,	2886.05,	2760.86,	3074.95,	3638.06,	4071.52,	4551.61,	5366.96)
Aged2_wm_x <- c(1520.73,	1065.33,	1040.08,	1111.54,	1414.43,	1814.84,	2495.04,	3274.32,	3729.1,	3191.53,	2745.77,	2303.29,	2941.61,	3838.56,	3184.93,	2343.08,	2197.1,	1837.01,	1466.61,	1201.17,	1208.9,	1647.37,	1902.23,	1936.8,	1655.56,	1593.49,	1697.49,	2219.56,	2689.21,	2914.55,	2929.08,	2481.27,	4420.67,	3850.27,	3223.55,	3617.13,	4153.89,	4338.72,	5450.15,	4883.34,	4445.27,	3697.26,	2958.15,	1961.6,	1590.44,	1291.96,	1506.36)
Aged2_wm_y <- c(5522.22,	4956.2,	4034.16,	3455.08,	2791.42,	2234.06,	1610.48,	1326.35,	1323.95,	1450.42,	1671.34,	1995.96,	2048.69,	2417.07,	2540.19,	3318.01,	2927.76,	2735.95,	3114.74,	3643.48,	4403.11,	5114.2,	5412.9,	5645.71,	6360.56,	6816.76,	7369.99,	8034.48,	8212.97,	7994.4,	7477.45,	7774.94,	7946.03,	8557.02,	8712.98,	8990.08,	9038.11,	8908.63,	8717.81,	8963.05,	9142.47,	9100.81,	8955.37,	8133.34,	7493.04,	6510.68,	5713.6)
Aged2_th_x <- c(3063.02,	2875.82,	2824.62,	2806.42,	3022.7,	3916.48,	5311.05,	5927.43,	6741.46,	6741.65,	6350.89,	5274.03,	3469.09,	2945.8,	2452.44,	2377.14,	2398.67,	2701.08,	3085.16)
Aged2_th_y <- c(5447.11,	6024.84,	6464.6,	6896.3,	7461.06,	7781.2,	6953.33,	6772.36,	5657.68,	4531.31,	3940.12,	3049.61,	2607.62,	2960.01,	3300.38,	3649.97,	3970.18,	4717.58,	5208.26)


Aged2_pol_cx <- point.in.polygon(Aged2.x, Aged2.y, Aged2_cx_x, Aged2_cx_y)
Aged2_pol_hp <- point.in.polygon(Aged2.x, Aged2.y, Aged2_hp_x, Aged2_hp_y)
Aged2_pol_wm <- point.in.polygon(Aged2.x, Aged2.y, Aged2_wm_x, Aged2_wm_y)
Aged2_pol_th <- point.in.polygon(Aged2.x, Aged2.y, Aged2_th_x, Aged2_th_y)

Aged2.info <- cbind.data.frame(Aged2_pol_cx,Aged2_pol_hp,Aged2_pol_wm,Aged2_pol_th)
head(Aged2.info, 50)



Aged2.info <- data.frame(Aged2.info)
Aged2@meta.data <- cbind(Aged2@meta.data,Aged2.info)


Aged2_cx <- subset(Aged2,subset = Aged2_pol_cx == 1)
Aged2_hp <- subset(Aged2,subset = Aged2_pol_hp == 1)
Aged2_wm <- subset(Aged2,subset = Aged2_pol_wm == 1)
Aged2_th <- subset(Aged2,subset = Aged2_pol_th == 1)

Aged2_cx$region <- "Cx"
Aged2_hp$region <- "Hp"
Aged2_wm$region <- "Wm"
Aged2_th$region <- "Th"

Aged2_all <- merge(Aged2_cx, y = c( Aged2_hp, Aged2_wm, Aged2_th))


Idents(Aged2_all)<-Aged2_all$region
levels(Aged2_all)

dfall<- data.frame(Aged2_all$center.x[WhichCells(object = subset(Aged2_all,subset = sample =="Aged_2"))],Aged2_all$center.y[WhichCells(object = subset(Aged2_all,subset = sample =="Aged_2"))],Aged2_all@active.ident[WhichCells(object = subset(Aged2_all,subset = sample =="Aged_2"))])
x <- Aged2_all$center.x[WhichCells(object = subset(Aged2_all,subset = sample =="Aged_2"))]
y <- Aged2_all$center.y[WhichCells(object = subset(Aged2_all,subset = sample =="Aged_2"))]
group <- Aged2_all@active.ident[WhichCells(object = subset(Aged2_all,subset = sample =="Aged_2"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#Aged_3

Aged3.x <- Aged3@meta.data[,5]
head(Aged3.x)
Aged3.y <- Aged3@meta.data[,6]
head(Aged3.y)

Aged3_cx_x <- c(1176.76,	1300.86,	1388.06,	2037.64,	3828.18,	4351.11,	6028.18,	7572.91,	8251.1,	8523.2,	8606.21,	7222.43,	4245.82,	2235.6,	1680.23,	1490.36,	1919.57,	2581.07,	3004.15,	3864.04,	4366.55,	4953.21,	5721.23,	6335.98,	6684.14,	7050.27,	7273.31,	7363.66,	6933.02,	6926.76,	7349.59,	7163.56,	6865.39,	6459.32,	5843.93,	5252.36,	4540.73,	2739.77,	2134.91,	1440.52,	1196.86)
Aged3_cx_y <- c(2368.16,	1730.47,	1091.35,	324.94,	109.92,	-58.7144,	169.085,	1142.55,	2735.62,	4395.32,	7953.17,	9920.08,	10935.2,	10198.2,	8949.29,	7868.46,	8221.47,	9045.47,	9371.36,	9572.76,	9701.32,	9548.33,	9178.42,	8714.03,	8206.52,	7528.72,	6650.3,	6107.32,	5269.06,	4962.05,	4362.91,	3607.93,	2900.42,	2300.19,	1742.15,	1371.48,	1160.02,	1482.43,	2031.64,	2960.29,	2464.43)
Aged3_hp_x <- c(6447.53,	6762.56,	6995.16,	6931.49,	6640.48,	6168.18,	5757.85,	5607.96,	5711.02,	5591.52,	5680.36,	6225.98,	6433.95,	5967.09,	5700.31,	5649.44,	5815.73,	5877.33,	5783.93,	5943.65,	6259.79,	6470.83,	6970.3,	7037.08,	6817.99,	6460.26)
Aged3_hp_y <- c(5060.86,	4564.49,	4003.62,	3514.37,	3014.69,	2477.35,	2369.08,	2551.89,	3498.27,	4537.59,	4816.02,	5071.71,	5193.83,	5429.87,	5679.71,	5958.65,	6837.51,	7499.02,	8050.83,	8090.38,	7990.13,	7808.57,	6648.58,	6371.4,	5865.13,	5393.53)
Aged3_wm_x <- c(6944.44,	7239.16,	7200.01,	7041.17,	6804.34,	6319.49,	5687.3,	4783.13,	4117.57,	2834.55,	2205.77,	1733.67,	2099.21,	2527.54,	3424.19,	3982.43,	3912.17,	4988.18,	5491.17,	5671.71,	5933.23,	6511.34,	6890.37,	7037.99,	6447.53,	6484.03,	7099.02,	7106.72,	6575.42,	5582.09,	4803.78,	4347.49,	3505.65,	2466.21,	2892.19,	3905.12,	4405.29,	5334.31,	6256.31,	7020.52,	7159.15,	7346.88,	7295.73,	6909.64)
Aged3_wm_y <- c(4862.9,	4278.31,	3848.37,	3271.66,	2819.8,	2175.57,	1631.57,	1221.27,	1179,	1484.93,	1999.53,	2631.81,	2298.88,	1871.81,	1587.41,	1719.84,	2322.78,	2507.74,	2752.71,	2443.51,	2345.74,	2780.16,	3369.26,	3823.95,	5060.86,	5389.47,	6213.81,	6604.94,	7898.7,	8612.5,	9272.55,	9447.78,	9342.45,	8747.84,	9274.63,	9536.93,	9662.12,	9366.48,	8772.19,	7765.09,	7151.66,	6209.63,	5854.57,	5291.78)
Aged3_th_x <- c(5606.23,	5708.96,	5563.49,	5031.71,	3412.99,	2304.84,	2079.81,	2572.86,	2417.35,	3790.3,	4611.16,	5071.36,	5760.89,	5811.42,	5949.44)
Aged3_th_y <- c(4755.48,	3528.7,	2842.82,	2525.53,	3097.5,	4108.76,	4709.94,	5319.22,	6265.89,	7685.58,	8193.13,	8131.19,	7525.87,	7333.18,	5958.65)


# find layer yes/no using point in polygon
Aged3_pol_cx <- point.in.polygon(Aged3.x, Aged3.y, Aged3_cx_x, Aged3_cx_y)
Aged3_pol_hp <- point.in.polygon(Aged3.x, Aged3.y, Aged3_hp_x, Aged3_hp_y)
Aged3_pol_wm <- point.in.polygon(Aged3.x, Aged3.y, Aged3_wm_x, Aged3_wm_y)
Aged3_pol_th <- point.in.polygon(Aged3.x, Aged3.y, Aged3_th_x, Aged3_th_y)

Aged3.info <- cbind.data.frame(Aged3_pol_cx,Aged3_pol_hp,Aged3_pol_wm,Aged3_pol_th)
head(Aged3.info, 50)



Aged3.info <- data.frame(Aged3.info)
Aged3@meta.data <- cbind(Aged3@meta.data,Aged3.info)


Aged3_cx <- subset(Aged3,subset = Aged3_pol_cx == 1)
Aged3_hp <- subset(Aged3,subset = Aged3_pol_hp == 1)
Aged3_wm <- subset(Aged3,subset = Aged3_pol_wm == 1)
Aged3_th <- subset(Aged3,subset = Aged3_pol_th == 1)

Aged3_cx$region <- "Cx"
Aged3_hp$region <- "Hp"
Aged3_wm$region <- "Wm"
Aged3_th$region <- "Th"

Aged3_all <- merge(Aged3_cx, y = c( Aged3_hp, Aged3_wm, Aged3_th))


Idents(Aged3_all)<-Aged3_all$region
levels(Aged3_all)

dfall<- data.frame(Aged3_all$center.x[WhichCells(object = subset(Aged3_all,subset = sample =="Aged_3"))],Aged3_all$center.y[WhichCells(object = subset(Aged3_all,subset = sample =="Aged_3"))],Aged3_all@active.ident[WhichCells(object = subset(Aged3_all,subset = sample =="Aged_3"))])
x <- Aged3_all$center.x[WhichCells(object = subset(Aged3_all,subset = sample =="Aged_3"))]
y <- Aged3_all$center.y[WhichCells(object = subset(Aged3_all,subset = sample =="Aged_3"))]
group <- Aged3_all@active.ident[WhichCells(object = subset(Aged3_all,subset = sample =="Aged_3"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()


#5xFAD_1

AD1.x <- AD1@meta.data[,5]
head(AD1.x)
AD1.y <- AD1@meta.data[,6]
head(AD1.y)

AD1_cx_x <- c(7374.62,	9117.71,	9594.49,	8176.72,	5557.02,	4978.53,	4700.28,	4268.89,	4900.38,	6307.33,	7053.57,	7539.83,	7133.12,	6603.14,	5906.64,	6492.8,	7284.1,	7604.7,	7835.23,	7925.25,	7900.77,	7892.3,	7522.29,	6984.77,	6249.41,	5957.36,	5803,	4921.89,	4147.35,	3337.73,	2664.55,	2073.73,	372.291,	963.667,	3042.24,	5342.74,	7013.89)
AD1_cx_y <- c(1701.89,	3713.88,	6313.6,	8711.9,	9086.89,	9003.17,	8749.28,	7700.76,	7158.84,	7141.74,	6633.63,	6796.64,	7296.52,	7522.89,	8215.94,	7812.68,	7311.94,	6961.55,	6529.36,	5966.41,	5303.65,	4680.77,	3876.32,	3170.82,	3118.54,	2880.94,	2183.32,	1707.35,	1430.27,	1395.47,	1524.68,	1750.08,	1727.27,	738.08,	64.4413,	296.103,	1699.05)
AD1_hp_x <- c(5575.58,	5514.1,	5205.87,	4027.28,	3815.04,	4230.69,	4604.14,	4835.22,	5349.07,	5549.15,	5665.3,	5802.15,	6311.49,	7020.9,	7422.39,	7600.96,	7599.49,	7353.57,	6957.6,	6325.42,	6057.32,	5921.15)
AD1_hp_y <- c(3135.29,	2341.05,	2053.05,	1703.44,	1865.1,	2535,	3035.95,	3496.98,	3627.06,	3747.73,	4179.84,	4329.52,	4483.54,	4931.31,	5118.72,	4978.14,	4670.31,	4007.76,	3560.72,	3383.81,	3348.55,	3272.9)
AD1_wm_x <- c(6025.66,	5784.91,	4856.3,	4117.33,	3002.74,	2197.56,	3615.87,	2954.79,	2852.49,	2989.47,	4084.85,	3783.48,	4051.63,	4822.14,	5442.89,	5775.58,	6015.31,	6861.36,	7246.5,	7628.45,	7575.81,	7313.35,	6701.99,	6845.7,	7200.62,	7713.23,	7732.35,	7565.9,	7121.53,	7591.49,	7760.77,	7938.95,	7968.45,	7792.97,	7347.31,	6881.14,	6219.09,	6014.41)
AD1_wm_y <- c(2981.52,	2204.09,	1691.69,	1431.3,	1445.15,	1775.42,	1650.86,	2093.25,	2467.96,	2596.14,	2470.32,	1869.54,	1684.31,	1863.28,	2195.32,	3135.29,	3319.98,	3483.41,	3833.49,	4784.31,	5061.72,	5100.17,	4746.75,	5580.7,	5660.46,	5414.15,	6275.13,	6869.43,	7343.39,	6989.77,	6567.58,	5946.33,	5409.42,	4497.13,	3683.96,	3155.59,	3125.35,	2950.8)
AD1_th_x <- c(5074.02,	4581.42,	4033.21,	3227.77,	2891.28,	2717.16,	2324.27,	2047.5,	3617.14,	6510.99,	6771.14,	6771.15,	6324.02,	5167.58)
AD1_th_y <- c(4089.83,	3000.98,	2533.09,	2671.94,	2989.54,	3348.26,	4590.68,	6126.02,	7698.59,	6204.95,	5886.06,	4959.17,	4539.12,	4086.83)


AD1_pol_cx <- point.in.polygon(AD1.x, AD1.y, AD1_cx_x, AD1_cx_y)
AD1_pol_hp <- point.in.polygon(AD1.x, AD1.y, AD1_hp_x, AD1_hp_y)
AD1_pol_wm <- point.in.polygon(AD1.x, AD1.y, AD1_wm_x, AD1_wm_y)
AD1_pol_th <- point.in.polygon(AD1.x, AD1.y, AD1_th_x, AD1_th_y)

AD1.info <- cbind.data.frame(AD1_pol_cx,AD1_pol_hp,AD1_pol_wm,AD1_pol_th)
head(AD1.info, 50)



AD1.info <- data.frame(AD1.info)
AD1@meta.data <- cbind(AD1@meta.data,AD1.info)


AD1_cx <- subset(AD1,subset = AD1_pol_cx == 1)
AD1_hp <- subset(AD1,subset = AD1_pol_hp == 1)
AD1_wm <- subset(AD1,subset = AD1_pol_wm == 1)
AD1_th <- subset(AD1,subset = AD1_pol_th == 1)

AD1_cx$region <- "Cx"
AD1_hp$region <- "Hp"
AD1_wm$region <- "Wm"
AD1_th$region <- "Th"

AD1_all <- merge(AD1_cx, y = c( AD1_hp, AD1_wm, AD1_th))


Idents(AD1_all)<-AD1_all$region
levels(AD1_all)

dfall<- data.frame(AD1_all$center.x[WhichCells(object = subset(AD1_all,subset = sample =="5xFAD_1"))],AD1_all$center.y[WhichCells(object = subset(AD1_all,subset = sample =="5xFAD_1"))],AD1_all@active.ident[WhichCells(object = subset(AD1_all,subset = sample =="5xFAD_1"))])
x <- AD1_all$center.x[WhichCells(object = subset(AD1_all,subset = sample =="5xFAD_1"))]
y <- AD1_all$center.y[WhichCells(object = subset(AD1_all,subset = sample =="5xFAD_1"))]
group <- AD1_all@active.ident[WhichCells(object = subset(AD1_all,subset = sample =="5xFAD_1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#5xFAD_2

AD2.x <- AD2@meta.data[,5]
head(AD2.x)
AD2.y <- AD2@meta.data[,6]
head(AD2.y)

AD2_cx_x <- c(1626.35,	209.656,	13.2123,	1645.82,	3998.42,	5679.83,	5231.16,	2874.36,	2345.96,	2830.15,	330.09,	2925.48,	2272.75,	1901.72,	1564.34,	1497.99,	1827.46,	2346.9,	2659.32,	2921.27,	3142.4,	3326.85,	3482.02,	4146.31,	5055.02,	6108.64,	6987.34,	7694.42,	8283.71,	7557.25,	7271.45,	7257.25,	7436.09,	8307.11,	10060.3,	8946.33,	5195.76,	2681.17,	1691.64)
AD2_cx_y <- c(7565.66,	4982.76,	2357.2,	888.802,	129.422,	913.401,	2418.94,	2306.14,	2036.65,	1826.31,	1381,	1665.05,	1909.47,	2338.8,	3047.38,	4354.78,	5365.57,	6120.72,	6193.97,	6306.94,	6588.18,	7310.15,	7412.68,	8009.85,	8388.03,	8522.87,	8260.11,	7656.15,	7019.41,	7644.43,	7860.04,	6308.4,	5240.45,	4095.36,	5747.89,	9604.58,	10213.1,	8758.93,	7569.99)
AD2_hp_x <- c(3386.73,	3499.19,	3647.31,	4038.86,	4931.5,	5147.87,	4551.46,	4356.25,	3933.64,	3746.12,	3620.85,	3539.79,	2854.4,	2038.06,	1853.4,	2016.86,	2353.31,	3336.41)
AD2_hp_y <- c(6254.51,	6942.95,	7353.8,	7708.76,	8045.84,	7873.78,	6889.82,	6271.22,	6026.33,	6008.76,	5771.76,	5272.53,	4820.52,	4225.94,	4423.76,	5318.19,	5761.69,	6162.23)
AD2_wm_x <- c(3069.41,	3287.46,	3512,	4126.03,	5016.33,	5980.34,	6683.06,	7313.81,	6599.2,	6342.36,	5049.77,	4893.59,	5148.87,	4957.26,	4115.09,	3637.87,	3394.57,	3268.47,	2459.73,	2016.86,	1849.44,	2088.74,	2580.59,	2527.73,	2010.56,	1798.46,	1749.18,	1970.53,	2366.48,	2504.62,	2097.85,	1870.43,	1650.75,	1606.67,	1554.8,	1565.78,	1850.8,	2240,	2329.04,	2991.91)
AD2_wm_y <- c(6452.81,	7248.04,	7426.74,	7981.19,	8354.62,	8524.23,	8408.61,	8008.41,	8239.82,	7243.88,	7385.335,	7369.57,	7873.78,	8088.72,	7749.07,	7299.78,	6361.86,	6146.13,	5893.68,	5318.19,	4276.48,	4178.98,	4557.17,	2676.74,	3149.65,	3827.56,	3059.02,	2428.8,	1981.36,	1858.75,	2041.52,	2400.8,	2872.08,	3223.97,	3986.71,	4545.39,	5412.26,	5956.21,	6120.07,	6320.36)
AD2_th_x <- c(4103.77,	4105.99,	4347.89,	4572.23,	5016.95,	6211.29,	7244.35,	7395.96,	8233.11,	5659.38,	3465.45,	2811.8,	2562.1,	2764.72,	3981.55,	4048.23)
AD2_th_y <- c(5683.98,	6079.31,	6248.25,	6859.75,	7358.72,	7302.47,	5005.45,	3929.07,	2949.28,	2329.89,	3122.03,	3179.39,	4298.98,	4724.29,	5466.93,	5595.86)


AD2_pol_cx <- point.in.polygon(AD2.x, AD2.y, AD2_cx_x, AD2_cx_y)
AD2_pol_hp <- point.in.polygon(AD2.x, AD2.y, AD2_hp_x, AD2_hp_y)
AD2_pol_wm <- point.in.polygon(AD2.x, AD2.y, AD2_wm_x, AD2_wm_y)
AD2_pol_th <- point.in.polygon(AD2.x, AD2.y, AD2_th_x, AD2_th_y)

AD2.info <- cbind.data.frame(AD2_pol_cx,AD2_pol_hp,AD2_pol_wm,AD2_pol_th)
head(AD2.info, 50)



AD2.info <- data.frame(AD2.info)
AD2@meta.data <- cbind(AD2@meta.data,AD2.info)


AD2_cx <- subset(AD2,subset = AD2_pol_cx == 1)
AD2_hp <- subset(AD2,subset = AD2_pol_hp == 1)
AD2_wm <- subset(AD2,subset = AD2_pol_wm == 1)
AD2_th <- subset(AD2,subset = AD2_pol_th == 1)

AD2_cx$region <- "Cx"
AD2_hp$region <- "Hp"
AD2_wm$region <- "Wm"
AD2_th$region <- "Th"

AD2_all <- merge(AD2_cx, y = c( AD2_hp, AD2_wm, AD2_th))


Idents(AD2_all)<-AD2_all$region
levels(AD2_all)

dfall<- data.frame(AD2_all$center.x[WhichCells(object = subset(AD2_all,subset = sample =="5xFAD_2"))],AD2_all$center.y[WhichCells(object = subset(AD2_all,subset = sample =="5xFAD_2"))],AD2_all@active.ident[WhichCells(object = subset(AD2_all,subset = sample =="5xFAD_2"))])
x <- AD2_all$center.x[WhichCells(object = subset(AD2_all,subset = sample =="5xFAD_2"))]
y <- AD2_all$center.y[WhichCells(object = subset(AD2_all,subset = sample =="5xFAD_2"))]
group <- AD2_all@active.ident[WhichCells(object = subset(AD2_all,subset = sample =="5xFAD_2"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()


#5xFAD_3

AD3.x <- AD3@meta.data[,5]
head(AD3.x)
AD3.y <- AD3@meta.data[,6]
head(AD3.y)

AD3_cx_x <- c(3417.45,	7337.21,	9571.91,	10396.6,	9615.22,	7715.09,	7694.82,	8287.77,	8578.29,	8904.58,	8717.72,	8219.85,	7308.53,	6670.15,	5747.58,	5384.67,	4924.09,	4141.79,	3808.46,	3156.78,	2475.85,	1699.95,	1432.03,	1446.26,	1789.35,	1523.82,	1918.95,	2351.25,	3751.45,	4897.05,	3040.98,	1785.42,	758.15,	115.922,	889.028,	3153.97)
AD3_cx_y <- c(1287.17,	143.433,	1245.01,	3767.07,	6411.18,	6257.22,	4815.56,	3547.87,	3291.8,	4059.53,	3301.76,	2454.3,	1890.97,	1684.96,	1696.64,	1752.28,	1983.35,	2675.14,	2803.84,	2715.79,	3301.78,	4299.08,	5271.18,	6260.03,	6766.55,	5802.65,	5865.23,	6512.04,	7020.14,	7961.26,	8902.05,	8624.69,	7615.17,	4519.41,	2796.54,	1308.4)
AD3_hp_x <- c(4256.87,	4992.34,	5783.74,	6444.87,	7048,	7104.38,	7027.74,	6483.73,	5897.32,	5082.3,	4832.99,	4623.39,	4338.23,	4092.44,	3794.59,	3708.68,	3120.26,	2255.39,	1923.61,	1940.97,	2254.07,	2869.47,	3274.15,	4141.01)
AD3_hp_y <- c(3108.33,	2117.93,	1864.95,	1883.66,	2206.28,	2406.68,	2530.39,	2726.07,	2946.34,	3580.79,	3583.83,	3485.4,	3603.81,	4026.66,	4146.55,	4128.21,	4314.49,	4862.52,	4764.78,	4435.21,	3798.91,	3197.46,	3062.46,	3165.86)
AD3_wm_x <- c(3918.59,	3247.83,	2388.22,	1836.79,	1516.14,	1367.03,	1445.83,	1894.74,	1518.11,	1508.49,	1976.16,	2527.6,	2520.66,	2084.89,	1902.21,	2015.2,	2517.18,	3078.42,	4104.73,	4318.56,	4992.34,	5924.18,	7031.79,	7132.75,	7005.58,	6582.67,	7464.32,	7885.16,	7660.58,	7751.94,	8509.33,	8754.22,	8322.69,	7866.3,	7220.25,	6422.36,	5767.15,	4531.43,	4198.98)
AD3_wm_y <- c(2832.13,	2775.65,	3407.73,	4126.1,	4864.09,	5542.7,	6233.15,	6944.39,	6061.13,	5770.39,	5900.76,	5885.72,	4699.94,	4873.35,	4657.7,	4104.54,	3472.27,	3099.77,	3117.46,	3041.92,	2117.93,	1836.64,	2095.89,	2392.66,	25476.01,	2717.11,	3384.54,	3161.58,	2515.99,	2298.86,	3054.72,	3555.8,	2666.58,	2245.15,	1872.7,	1699.94,	1711.02,	2197.24,	2657.33)
AD3_th_x <- c(4234.01,	3821.84,	3141.7,	2617.67,	2727.4,	4387.59,	5599.95,	7639.85,	7093.85,	7406.34,	6394.45,	5937.25,	4843.41,	4793.49)
AD3_th_y <- c(3837.7,	4162.88,	4322.84,	4706.66,	5760.11,	6609.62,	8212.66,	7144.82,	4645.14,	3557.62,	2790.12,	2952.47,	3749.63,	4153.81)


AD3_pol_cx <- point.in.polygon(AD3.x, AD3.y, AD3_cx_x, AD3_cx_y)
AD3_pol_hp <- point.in.polygon(AD3.x, AD3.y, AD3_hp_x, AD3_hp_y)
AD3_pol_wm <- point.in.polygon(AD3.x, AD3.y, AD3_wm_x, AD3_wm_y)
AD3_pol_th <- point.in.polygon(AD3.x, AD3.y, AD3_th_x, AD3_th_y)

AD3.info <- cbind.data.frame(AD3_pol_cx,AD3_pol_hp,AD3_pol_wm,AD3_pol_th)
head(AD3.info, 50)



AD3.info <- data.frame(AD3.info)
AD3@meta.data <- cbind(AD3@meta.data,AD3.info)


AD3_cx <- subset(AD3,subset = AD3_pol_cx == 1)
AD3_hp <- subset(AD3,subset = AD3_pol_hp == 1)
AD3_wm <- subset(AD3,subset = AD3_pol_wm == 1)
AD3_th <- subset(AD3,subset = AD3_pol_th == 1)

AD3_cx$region <- "Cx"
AD3_hp$region <- "Hp"
AD3_wm$region <- "Wm"
AD3_th$region <- "Th"

AD3_all <- merge(AD3_cx, y = c( AD3_hp, AD3_wm, AD3_th))


Idents(AD3_all)<-AD3_all$region
levels(AD3_all)

dfall<- data.frame(AD3_all$center.x[WhichCells(object = subset(AD3_all,subset = sample =="5xFAD_3"))],AD3_all$center.y[WhichCells(object = subset(AD3_all,subset = sample =="5xFAD_3"))],AD3_all@active.ident[WhichCells(object = subset(AD3_all,subset = sample =="5xFAD_3"))])
x <- AD3_all$center.x[WhichCells(object = subset(AD3_all,subset = sample =="5xFAD_3"))]
y <- AD3_all$center.y[WhichCells(object = subset(AD3_all,subset = sample =="5xFAD_3"))]
group <- AD3_all@active.ident[WhichCells(object = subset(AD3_all,subset = sample =="5xFAD_3"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()


#WT-5xFAD_1

WT1.x <- WT1@meta.data[,5]
head(WT1.x)
WT1.y <- WT1@meta.data[,6]
head(WT1.y)

WT1_cx_x <- c(7198.98,	5563.76,	3907.76,	2810.74,	2152.36,	479.183,	-58.0796,	-30.6971,	1033.93,	1835.56,	2612.7,	3189.63,	3352.58,	3106.69,	2684.75,	1888.19,	2679.44,	3111.74,	3374.06,	4156.91,	4884.92,	5515.68,	5645.57,	5837.82,	5550.52,	5867.42,	6162.81,	6254.58,	6236.05,	5691.9,	5122.33,	4356.45,	3781.91,	3175.82,	3976.96,	4750.06,	4054.39,	3316.53,	1617.5,	908.426,	1830.23,	4004.45,	5675.21,	7077.81,	7198.98)
WT1_cx_y <- c(5694.38,	9372.21,	9938.86,	10010.3,	9830.3,	8899.65,	8293.16,	7673.09,	6476.21,	7530.16,	8288.41,	8530.04,	8775.63,	8904.75,	8858.47,	8661.32,	8946.85,	9001.66,	8975.31,	8653.18,	8068.58,	7376.32,	7002.42,	6502.37,	5770.91,	4956.72,	4649.13,	3778.36,	3427.29,	2374.4,	1704.91,	1246.07,	1114.14,	1162.42,	1219.57,	1568.39,	1796.59,	1993.41,	3118.07,	1484.13,	311.756,	162.714,	260.734,	2454.97,	5694.38)
WT1_hp_x <- c(5203.97,	5962.23,	6097.57,	6039.5,	5382.73,	4917.14,	5019.68,	4621.81,	4814.29,	4747.71,	4353.4,	4307.45,	3990.55,	4094.91,	4322.72,	4842.67,	5318.34,	5518.56,	5515.68,	5175.58)
WT1_hp_y <- c(5239.01,	4373.72,	4096.71,	3511.67,	2373.03,	2473.77,	3541.35,	4717.02,	5090.77,	5558.34,	5927.43,	6977.89,	7935.61,	8202.62,	8181.77,	7919.58,	7357.62,	6887.5,	6473.4,	5564.76)
WT1_wm_x <- c(5720.74,	6137.52,	6225.64,	6214.8,	5691.9,	5122.33,	4125.52,	4702.88,	4742.8,	3912.33,	3150.64,	4711.78,	4934.45,	4936.49,	5122.25,	5382.73,	5916.75,	6114.18,	5256.53,	5212.75,	5539.97,	5332.55,	4385.85,	4059.43,	3981.09,	4039.53,	2312.53,	3351.29,	3106.69,	2601.66,	2273.96,	3061.63,	3634.94,	4883.21,	5457.14,	5663.38,	5799.73,	5608.8,	5640.26)
WT1_wm_y <- c(5161.73,	4780.58,	3891.76,	3378.73,	2374.4,	1704.91,	1169.06,	1563.12,	1804.24,	1809.52,	2069.84,	2677.05,	2842.62,	2438.4,	2317.98,	2379.03,	3152.38,	3942.37,	5230.73,	5614.16,	6652.92,	7316.61,	8173.24,	8183.84,	7897.52,	7663.9,	7996.46,	8717.51,	8904.75,	8862.68,	8815.68,	8997.16,	8877.55,	8051.54,	7397.6,	6754.79,	6329.08,	5851.63,	5610.35)
WT1_th_x <- c(4157.76,	4635.84,	4620.88,	4961.08,	4872.3,	3835.55,	684.38,	188.534,	2658.57,	3854.72,	4101.31,	4277.81,	4419.94,	4168.51)
WT1_th_y <- c(5222,	4853.69,	4638.96,	3600.39,	2881.19,	2411.73,	3675.23,	4765.98,	7569.54,	7748.76,	7542.55,	7000.66,	5648.39,	5366.06)


WT1_pol_cx <- point.in.polygon(WT1.x, WT1.y, WT1_cx_x, WT1_cx_y)
WT1_pol_hp <- point.in.polygon(WT1.x, WT1.y, WT1_hp_x, WT1_hp_y)
WT1_pol_wm <- point.in.polygon(WT1.x, WT1.y, WT1_wm_x, WT1_wm_y)
WT1_pol_th <- point.in.polygon(WT1.x, WT1.y, WT1_th_x, WT1_th_y)

WT1.info <- cbind.data.frame(WT1_pol_cx,WT1_pol_hp,WT1_pol_wm,WT1_pol_th)
head(WT1.info, 50)



WT1.info <- data.frame(WT1.info)
WT1@meta.data <- cbind(WT1@meta.data,WT1.info)


WT1_cx <- subset(WT1,subset = WT1_pol_cx == 1)
WT1_hp <- subset(WT1,subset = WT1_pol_hp == 1)
WT1_wm <- subset(WT1,subset = WT1_pol_wm == 1)
WT1_th <- subset(WT1,subset = WT1_pol_th == 1)

WT1_cx$region <- "Cx"
WT1_hp$region <- "Hp"
WT1_wm$region <- "Wm"
WT1_th$region <- "Th"

WT1_all <- merge(WT1_cx, y = c( WT1_hp, WT1_wm, WT1_th))


Idents(WT1_all)<-WT1_all$region
levels(WT1_all)

dfall<- data.frame(WT1_all$center.x[WhichCells(object = subset(WT1_all,subset = sample =="WT-5xFAD_1"))],WT1_all$center.y[WhichCells(object = subset(WT1_all,subset = sample =="WT-5xFAD_1"))],WT1_all@active.ident[WhichCells(object = subset(WT1_all,subset = sample =="WT-5xFAD_1"))])
x <- WT1_all$center.x[WhichCells(object = subset(WT1_all,subset = sample =="WT-5xFAD_1"))]
y <- WT1_all$center.y[WhichCells(object = subset(WT1_all,subset = sample =="WT-5xFAD_1"))]
group <- WT1_all@active.ident[WhichCells(object = subset(WT1_all,subset = sample =="WT-5xFAD_1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#WT-5xFAD_2

WT2.x <- WT2@meta.data[,5]
head(WT2.x)
WT2.y <- WT2@meta.data[,6]
head(WT2.y)

WT2_cx_x <- c(7315.5,	7873.64,	5816.47,	3192.24,	2157.93,	2727.24,	2938.6,	3775.59,	4355.04,	5105.02,	5588.38,	5744.36,	5270.39,	4886.37,	5552.18,	5896.68,	6355.35,	6619.73,	6817.08,	6742.81,	6635.4,	5909.21,	6129.49,	5917.78,	6059.92,	5853.01,	5309.5,	4571.07,	3649.66,	2696.28,	2078.86,	2763.2,	3193.35,	3060.7,	2452.75,	2019.76,	1492.75,	-165.3,1666.96,	4237.19,	7181.85)
WT2_cx_y <- c(4339.93,	8575.56,	10886.8,	10802,	9434.8,	7723.51,	7726.39,	8366,	9087.96,	8984.14,	8749.53,	8831.84,	9224.27,	9424.98,	9306.03,	8972.79,	8305.66,	7701.23,	6715.16,	6030.46,	5547.89,	4981.33,	4709.2,	4260.79,	3722.35,	3348.04,	2517.55,	2005.44,	1645.37,	1503.89,	1759.5,	1648.29,	1746.73,	2021.11,	2488.67,	3231.37,	4605.52,	3381.84,	138.417,	139.603,	3058.77)
WT2_hp_x <- c(5582.66,	6546.49,	6664.98,	6509.43,	5933.82,	5800.6,	5701.23,	5734.06,	5780.18,	5545.57,	5158.46,	5103.63,	5227.28,	5097.26,	4720.78,	4454.18,	4008.72,	4036.54,	4286.2,	4840.39,	5415.54,	5682.74,	5472.61)
WT2_hp_y <- c(4965.5,	5927.31,	6681.86,	7466.6,	8250.66,	8280.95,	8112.69,	7554.24,	7127.49,	6461.07,	5701.29,	5332.59,	5132.54,	4750.58,	4536.15,	3447.99,	2644.07,	2315.02,	2235.65,	2441.77,	2913.42,	3424.86,	4625.4)
WT2_wm_x <- c(6113.5,	6532.98,	6823.08,	6578.45,	6190.95,	5655.73,	5066.96,	5675.74,	5872.89,	5133.85,	4411.21,	5212.81,	5661.48,	5675.88,	6049.58,	6332.83,	6606.31,	6657.28,	6544.44,	5615.92,	5441.94,	5699.45,	5470.84,	4568.96,	4101.23,	3959.59,	4159.49,	2711.29,	2536.13,	3648.15,	2836,	2236,	2696.28,	3098.43,	4182.09,	5040.93,	5626.63,	6019.69,	5857.87,	6129.49)
WT2_wm_y <- c(5121.08,	5440.63,	6602.68,	7743.8,	8533.11,	9234.91,	9413.04,	8982.56,	8639.31,	8931.65,	8445.81,	8162.36,	7668.99,	8297.84,	8278.08,	7902.39,	7229.35,	6384.57,	5899.34,	5025.44,	4598.29,	3607.22,	2970.66,	2269.93,	2231.01,	2426.58,	2865.85,	2741.87,	2488.94,	1835.56,	1626.87,	1760.31,	1503.89,	1516.58,	1891.8,	2386.91,	2900.41,	3724.56,	4026.92,	4709.2)
WT2_th_x <- c(4658.82,	4754.42,	4435.55,	4112.13,	2998.86,	2048.71,	933.013,	2384.63,	2997.33,	4885.25,	5447.84,	5732.16,	5690.09,	5060.93,	4745.03)
WT2_th_y <- c(5085.1,	4606.14,	3386.86,	2864.71,	2868.58,	4698.78,	5078.82,	7632.54,	7155.94,	8236.74,	7856.49,	7436.58,	6917.15,	5406.6,	5349.22)


WT2_pol_cx <- point.in.polygon(WT2.x, WT2.y, WT2_cx_x, WT2_cx_y)
WT2_pol_hp <- point.in.polygon(WT2.x, WT2.y, WT2_hp_x, WT2_hp_y)
WT2_pol_wm <- point.in.polygon(WT2.x, WT2.y, WT2_wm_x, WT2_wm_y)
WT2_pol_th <- point.in.polygon(WT2.x, WT2.y, WT2_th_x, WT2_th_y)

WT2.info <- cbind.data.frame(WT2_pol_cx,WT2_pol_hp,WT2_pol_wm,WT2_pol_th)
head(WT2.info, 50)



WT2.info <- data.frame(WT2.info)
WT2@meta.data <- cbind(WT2@meta.data,WT2.info)


WT2_cx <- subset(WT2,subset = WT2_pol_cx == 1)
WT2_hp <- subset(WT2,subset = WT2_pol_hp == 1)
WT2_wm <- subset(WT2,subset = WT2_pol_wm == 1)
WT2_th <- subset(WT2,subset = WT2_pol_th == 1)

WT2_cx$region <- "Cx"
WT2_hp$region <- "Hp"
WT2_wm$region <- "Wm"
WT2_th$region <- "Th"

WT2_all <- merge(WT2_cx, y = c( WT2_hp, WT2_wm, WT2_th))


Idents(WT2_all)<-WT2_all$region
levels(WT2_all)

dfall<- data.frame(WT2_all$center.x[WhichCells(object = subset(WT2_all,subset = sample =="WT-5xFAD_2"))],WT2_all$center.y[WhichCells(object = subset(WT2_all,subset = sample =="WT-5xFAD_2"))],WT2_all@active.ident[WhichCells(object = subset(WT2_all,subset = sample =="WT-5xFAD_2"))])
x <- WT2_all$center.x[WhichCells(object = subset(WT2_all,subset = sample =="WT-5xFAD_2"))]
y <- WT2_all$center.y[WhichCells(object = subset(WT2_all,subset = sample =="WT-5xFAD_2"))]
group <- WT2_all@active.ident[WhichCells(object = subset(WT2_all,subset = sample =="WT-5xFAD_2"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#WT-5xFAD_3

WT3.x <- WT3@meta.data[,5]
head(WT3.x)
WT3.y <- WT3@meta.data[,6]
head(WT3.y)

WT3_cx_x <- c(7346.26,	9383.82,	8440.56,	6649.57,	4412.58,	3743.65,	4044.54,	6023.73,	6711.93,	7194.63,	7269.16,	6879.96,	7306,	7586.46,	7751.78,	7757.33,	7390.57,	6948.13,	6247.9,	6066.7,	5831.59,	5331.32,	4718.81,	3518.72,	2361.05,	1713.2,	2245.86,	2858.75,	2136.32,	1992.93,	1687.8,	1087.69,	102.084,	823.714,	3002.49,	5573.25,	7346.26)
WT3_cx_y <- c(2751,	6195.21,	8897.47,	9328.11,	9402.09,	7833.23,	7409.97,	7722.75,	7559.4,	7340.45,	7507.94,	8044.32,	7680.56,	7139.59,	6614.79,	5785.5,	4802.15,	4170.89,	3836.45,	3543.95,	2836.02,	2387.3,	1895.16,	1668.58,	1675.36,	2327.57,	1896.51,	1910.53,	2825.72,	3948.01,	4842.73,	5225.06,	3662.13,	1178.45,	-20.1032,	1014.03,	2751)
WT3_hp_x <- c(5836.8,	6471.74,	6827.47,	7238.99,	7514.95,	7550.98,	7421.55,	7144.64,	6896.56,	6724.36,	6678.77,	6201.08,	5688.95,	5576,	5639.89,	5336.83,	4968.57,	4660.08,	4350.47,	3873.94,	3574.02,	3592.88,	376145,	4491.27,	5038.07,	5470.88,	5658.55,	5785.55)
WT3_hp_y <- c(4086.46,	4354.38,	4453.83,	4871.62,	5630.28,	6020.39,	6635.91,	6790.25,	6647.94,	6157.12,	5958.37,	5358.43,	4894.5,	4668.78,	4356.23,	4130.99,	4036.26,	3550.42,	3095.44,	2624.22,	2301.76,	2108.83,	1950.36,	2058.95,	2323.2,	2761.21,	3477.05,	3973.23)
WT3_wm_x <- c(6079.07,	5934.79,	5573.96,	4992.7,	4359.6,	3834.02,	3157.87,	2728.63,	2409.39,	3152.09,	2524.94,	2662.38,	3716.44,	3565.94,	3745.12,	4230.02,	4872.75,	5298.94,	5547.7,	5768.51,	5999.93,	6331.5,	6827.47,	7171.8,	7551.51,	7477.14,	7312.37,	7125.25,	6930.96,	6789.78,	6195.29,	6456.1,	7066.42,	7417.02,	7269.16,	6952.3,	7397.14,	7600.91,	7747.84,	7718.05,	7628.83,	7345.75,	7005.98,	6936.29,	6236.27)
WT3_wm_y <- c(3699.09,	3032.81,	2698.01,	2128.26,	1806.75,	1682.15,	1642.19,	1685.95,	1830.29,	1797.96,	2532.24,	2875.93,	2592.76,	2252.02,	1956.8,	1960.31,	2218.13,	2560.81,	2880.49,	3863.64,	4135.5,	4324.1,	4453.83,	4761.2,	5651.65,	6560.46,	6760.29,	6812.86,	6710.66,	6483.28,	7185.09,	7515.71,	7453.19,	6997.11,	7507.94,	8030.01,	7554.96,	6983.29,	6434.09,	5819.07,	5443.6,	4835.08,	4407.8,	4222.62,	3862.21)
WT3_th_x <- c(5067.06,	4885.14,	4660.08,	4356.62,	3813.43,	2967.82,	2636.65,	1924.82,	1849.81,	2319.26,	3660.95,	4864.61,	6120.17,	6725.63,	6685.08,	6114.52,	5560.95,	4978.93)
WT3_th_y <- c(4085.14,	3964.46,	3550.42,	3113.35,	2612.55,	2953.82,	4272.2,	5369.62,	5736.38,	6524.89,	7280.3,	6987.48,	7175.45,	6444.56,	5981.57,	5271.82,	4713.08,	4573.51)


WT3_pol_cx <- point.in.polygon(WT3.x, WT3.y, WT3_cx_x, WT3_cx_y)
WT3_pol_hp <- point.in.polygon(WT3.x, WT3.y, WT3_hp_x, WT3_hp_y)
WT3_pol_wm <- point.in.polygon(WT3.x, WT3.y, WT3_wm_x, WT3_wm_y)
WT3_pol_th <- point.in.polygon(WT3.x, WT3.y, WT3_th_x, WT3_th_y)

WT3.info <- cbind.data.frame(WT3_pol_cx,WT3_pol_hp,WT3_pol_wm,WT3_pol_th)
head(WT3.info, 50)



WT3.info <- data.frame(WT3.info)
WT3@meta.data <- cbind(WT3@meta.data,WT3.info)


WT3_cx <- subset(WT3,subset = WT3_pol_cx == 1)
WT3_hp <- subset(WT3,subset = WT3_pol_hp == 1)
WT3_wm <- subset(WT3,subset = WT3_pol_wm == 1)
WT3_th <- subset(WT3,subset = WT3_pol_th == 1)

WT3_cx$region <- "Cx"
WT3_hp$region <- "Hp"
WT3_wm$region <- "Wm"
WT3_th$region <- "Th"

WT3_all <- merge(WT3_cx, y = c( WT3_hp, WT3_wm, WT3_th))


Idents(WT3_all)<-WT3_all$region
levels(WT3_all)

dfall<- data.frame(WT3_all$center.x[WhichCells(object = subset(WT3_all,subset = sample =="WT-5xFAD_3"))],WT3_all$center.y[WhichCells(object = subset(WT3_all,subset = sample =="WT-5xFAD_3"))],WT3_all@active.ident[WhichCells(object = subset(WT3_all,subset = sample =="WT-5xFAD_3"))])
x <- WT3_all$center.x[WhichCells(object = subset(WT3_all,subset = sample =="WT-5xFAD_3"))]
y <- WT3_all$center.y[WhichCells(object = subset(WT3_all,subset = sample =="WT-5xFAD_3"))]
group <- WT3_all@active.ident[WhichCells(object = subset(WT3_all,subset = sample =="WT-5xFAD_3"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#CFA_1

CFA1.x <- CFA1@meta.data[,5]
head(CFA1.x)
CFA1.y <- CFA1@meta.data[,6]
head(CFA1.y)

CFA1_cx_x <- c(7685.95,	6119.55,	4330.95,	2203.96,	395.24,	-33.3766,	366.337,	1278.25,	1773.26,	1657.01,	2555.92,	2548.08,	1426.76,	2077.2,	2644.19,	2986.45,	3689.26,	4369.09,	5298.02,	5752.19,	6002.42,	5891.72,	6137.12,	6747.15,	6820.07,	6986.99,	7068.71,	7006.24,	6554.21,	6133.92,	5026.51,	5987.37,	4242.76,	2559.11,	2506.18,	2514.42,	2643.05,	2803.77,	5092.96,	6560.24,	6887.8,	7219.44,	7876.06,	8361.05)
CFA1_cx_y <- c(6668.74,	8945.71,	9915.03,	9941.51,	8757.4,	8110.75,	6201.96,	5589.79,	6823.04,	7102.69,	8168.28,	8487.72,	8132.48,	8496.2,	8643.15,	8686.2,	8658.03,	8465.34,	7878.77,	7196.7,	6805.94,	6106.35,	5872.17,	5497.56,	5071.64,	4630.89,	3902.64,	3286.55,	2474.89,	1803.79,	1224.63,	1895.06,	1985.66,	2583.36,	1257.18,	995.796,	672.323,	532.336,	196.295,	773.337,	964.737,	1346.58,	2330.89,	5458.06)
CFA1_hp_x <- c(5573.93,	5489.22,	5634.71,	5244.26,	4513.118,	3662.19,	3319.73,	3201.99,	3698.61,	4095.49,	4269.7,	4703.5,	4808.87,	5119.02,	5411.33,	5228.27,	5245.68,	5791.81,	6022.31,	6099.35,	6298.99,	6442.06,	6645.88,	6814.72,	6770.77,	6524.95,	5685.56)
CFA1_hp_y <- c(5880.66,	6064.93,	6894.79,	7678.64,	8159.14,	8361.48,	8225.85,	8005.01,	7688.67,	7437.77,	7163.3,	6061.51,	5960.63,	5821.62,	5349.44,	5032.92,	4962.52,	4098.56,	3305.76,	2924.87,	2932.91,	2994.08,	3293.78,	3850.22,	4741.81,	5146.75,	5603.19)
CFA1_wm_x <- c(6091.73,	6654.6,	6841.76,	6996.46,	6999.01,	6924.89,	6714.56,	6120.71,	5446.02,	6118.83,	5385.21,	5387.13,	5974.31,	6104.01,	6325.9,	6650.34,	6856.09,	6897.9,	6756.7,	6419.06,	5659.6,	5502.13,	5645.54,	5399.1,	4887.95,	4205.84,	3541.06,	3220.34,	3201.99,	3672.81,	2568.61,	2109.17,	3201.29,	2377.5,	1666.95,	2315.98,	2829.74,	3371.14,	3876.62,	4605.22,	5181.18,	5590.24,	5922.93,	5824.55)
CFA1_wm_y <- c(5896.72,	5477.92,	4974.14,	4468.59,	3739.16,	3299.77,	2773.89,	1852.97,	1421.85,	1977.04,	1981.24,	2678.82,	3311.96,	2859.8,	2906.19,	3215.32,	3688.25,	4262.43,	4831.55,	5255.1,	5635.99,	5974.89,	6874.04,	7507.79,	8020.45,	8320.44,	8370.45,	8160.43,	8005.01,	7682.56,	7434.67,	7638.93,	8440.45,	8519.41,	8231.61,	8537.71,	8670.26,	8653.65,	8625.16,	8360.95,	7957.27,	7396.63,	6842.76,	6128.27)
CFA1_th_x <- c(4655.14,	5227.52,	5208.11,	5798.65,	5951.57,	5879.21,	4933.78,	3034.19,	2574.94,	1749.07,	1471.86,	1164.5,	963.418,	1255.45,	1859.82,	2083.14,	2442.5,	3422.14,	3747.56,	3945.2,	4162.02,	4744.3,	4600.74)
CFA1_th_y <- c(5262.12,	5093.95,	4978.19,	4069.82,	3383.78,	3124.96,	2490.07,	3115.73,	2808.07,	2913.49,	3144.34,	3530.77,	4388.3,	5386.29,	5498.31,	6205.53,	7180.76,	7634.75,	7624.94,	7525.62,	7316.62,	5935.05,	5447.5)


CFA1_pol_cx <- point.in.polygon(CFA1.x, CFA1.y, CFA1_cx_x, CFA1_cx_y)
CFA1_pol_hp <- point.in.polygon(CFA1.x, CFA1.y, CFA1_hp_x, CFA1_hp_y)
CFA1_pol_wm <- point.in.polygon(CFA1.x, CFA1.y, CFA1_wm_x, CFA1_wm_y)
CFA1_pol_th <- point.in.polygon(CFA1.x, CFA1.y, CFA1_th_x, CFA1_th_y)

CFA1.info <- cbind.data.frame(CFA1_pol_cx,CFA1_pol_hp,CFA1_pol_wm,CFA1_pol_th)
head(CFA1.info, 50)



CFA1.info <- data.frame(CFA1.info)
CFA1@meta.data <- cbind(CFA1@meta.data,CFA1.info)


CFA1_cx <- subset(CFA1,subset = CFA1_pol_cx == 1)
CFA1_hp <- subset(CFA1,subset = CFA1_pol_hp == 1)
CFA1_wm <- subset(CFA1,subset = CFA1_pol_wm == 1)
CFA1_th <- subset(CFA1,subset = CFA1_pol_th == 1)

CFA1_cx$region <- "Cx"
CFA1_hp$region <- "Hp"
CFA1_wm$region <- "Wm"
CFA1_th$region <- "Th"

CFA1_all <- merge(CFA1_cx, y = c( CFA1_hp, CFA1_wm, CFA1_th))


Idents(CFA1_all)<-CFA1_all$region
levels(CFA1_all)

dfall<- data.frame(CFA1_all$center.x[WhichCells(object = subset(CFA1_all,subset = sample =="CFA_1"))],CFA1_all$center.y[WhichCells(object = subset(CFA1_all,subset = sample =="CFA_1"))],CFA1_all@active.ident[WhichCells(object = subset(CFA1_all,subset = sample =="CFA_1"))])
x <- CFA1_all$center.x[WhichCells(object = subset(CFA1_all,subset = sample =="CFA_1"))]
y <- CFA1_all$center.y[WhichCells(object = subset(CFA1_all,subset = sample =="CFA_1"))]
group <- CFA1_all@active.ident[WhichCells(object = subset(CFA1_all,subset = sample =="CFA_1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#CFA_2

CFA2.x <- CFA2@meta.data[,5]
head(CFA2.x)
CFA2.y <- CFA2@meta.data[,6]
head(CFA2.y)

CFA2_cx_x <- c(1085.04,	76.942,	126.892,	1963.72,	5823.2,	6611.79,	5611.18,	5210.48,	3221.6,	2142.1,	2349.15,	2540.04,	2102.63,	1551.45,	1292.93,	1236.27,	1321.15,	1622.39,	1651.97,	2324.81,	2486.52,	2385.48,	2702.83,	3725.23,	4784.1,	5419.86,	6274.23,	6777.84,	5970.09,	6572.11,	7058.54,	8166.7,	8940.32,	7477.6,	5626.41,	4245.15,	1870.09)
CFA2_cx_y <- c(3537.35,	5416.33,	8061.51,	9850.86,	10403.7,	8768.11,	7533.02,	7475.57,	8196.51,	8064.1,	8485.43,	8635.31,	8470.23,	7595.85,	6630.88,	6135.29,	5493.08,	4793.86,	4619.69,	4282.65,	4018.63,	3250.11,	2873.98,	1803.71,	1536.16,	1494.91,	1714.05,	2007.19,	1730.78,	2471.61,	5050.9,	4894.8,	2449.93,	650.228,	-130.926,	-106.465,	1806.29)
CFA2_hp_x <- c(2713.12,	1788.46,	1444.52,	1449.83,	1659.93,	1998.27,	2219.95,	2458.17,	2886.44,	3091.98,	3089.82,	3036.74,	3166.76,	3600.77,	3705.75,	3961.2,	4119.82,	4262.81,	4632.46,	5130.4,	5113.22,	4712.66,	4278.75,	3512.04,	2904.53,	2700.02,	2800.14,	2807.76)
CFA2_hp_y <- c(4446.65,	4914.31,	5523.02,	6354.14,	7084.58,	7174.65,	6562.28,	6075.86,	5523.57,	5256.86,	4963.8,	4759.06,	4432.31,	4299.67,	4143.58,	3423.02,	3057.92,	2852.94,	2550.46,	2233.67,	2004.16,	1868.61,	1872.82,	2225.32,	2792.55,	3312.63,	3935.09,	4359.83)
CFA2_wm_x <- c(2320.94,	1713.64,	1341,	1245.5,	1376.81,	1905.84,	2403.27,	1939.4,	2960.01,	3093.31,	2232.37,	2188.87,	2024.74,	1691.86,	1466.8,	1416.73,	1460.67,	1707.95,	2707.24,	2809.74,	2706.14,	2890.72,	3504.33,	4408.25,	5056.8,	5171.46,	4568.11,	5360.09,	5958.55,	6368.54,	6539.27,	6470.91,	6130.12,	5065.73,	4981.32,	5581.23,	6269.16,	5906.58,	5629.72,	5157.15,	4441.73,	3617.02,	3039.34,	2598.83,	2442,	2469.92)
CFA2_wm_y <- c(4326.17,	4628.22,	5487.98,	6148.66,	7041.87,	8140.38,	8579.68,	7938.6,	8210.48,	7897.03,	6958.2,	6794.33,	7148.04,	7121.51,	6574.08,	5867.22,	5449.05,	5010.65,	4486.74,	4169.04,	3234.13,	2782.79,	2191.86,	1852.43,	1991.27,	2149.54,	2599.15,	2621.91,	3035.37,	3204.35,	2813.92,	2385.21,	2142.18,	1883.03,	1728.56,	1620.08,	1784.34,	1587.35,	1516.65,	1495.23,	1608.91,	1921.47,	2447.32,	3071.88,	3266.8,	4088.1)
CFA2_th_x <- c(3641.36,	3098.22,	3072.17,	2465.61,	2228.94,	2232.37,	2765.75,	3329.58,	3789.05,	5118.63,	6194.48,	7365.65,	6302.46,	5989.56,	5641.91,	5283.49,	4813.61,	4352.99,	3924.17,	3599.12,	3697.82)
CFA2_th_y <- c(4931.06,	5124.2,	5292.43,	6073.33,	6616.17,	6958.2,	7520.1,	7835.01,	7596.48,	7142.8,	7799.58,	5502.42,	4822.65,	3162.34,	2817.06,	2636.67,	2577.5,	2793.26,	3552.43,	4356.89,	4818.07)


CFA2_pol_cx <- point.in.polygon(CFA2.x, CFA2.y, CFA2_cx_x, CFA2_cx_y)
CFA2_pol_hp <- point.in.polygon(CFA2.x, CFA2.y, CFA2_hp_x, CFA2_hp_y)
CFA2_pol_wm <- point.in.polygon(CFA2.x, CFA2.y, CFA2_wm_x, CFA2_wm_y)
CFA2_pol_th <- point.in.polygon(CFA2.x, CFA2.y, CFA2_th_x, CFA2_th_y)

CFA2.info <- cbind.data.frame(CFA2_pol_cx,CFA2_pol_hp,CFA2_pol_wm,CFA2_pol_th)
head(CFA2.info, 50)



CFA2.info <- data.frame(CFA2.info)
CFA2@meta.data <- cbind(CFA2@meta.data,CFA2.info)


CFA2_cx <- subset(CFA2,subset = CFA2_pol_cx == 1)
CFA2_hp <- subset(CFA2,subset = CFA2_pol_hp == 1)
CFA2_wm <- subset(CFA2,subset = CFA2_pol_wm == 1)
CFA2_th <- subset(CFA2,subset = CFA2_pol_th == 1)

CFA2_cx$region <- "Cx"
CFA2_hp$region <- "Hp"
CFA2_wm$region <- "Wm"
CFA2_th$region <- "Th"

CFA2_all <- merge(CFA2_cx, y = c( CFA2_hp, CFA2_wm, CFA2_th))


Idents(CFA2_all)<-CFA2_all$region
levels(CFA2_all)

dfall<- data.frame(CFA2_all$center.x[WhichCells(object = subset(CFA2_all,subset = sample =="CFA_2"))],CFA2_all$center.y[WhichCells(object = subset(CFA2_all,subset = sample =="CFA_2"))],CFA2_all@active.ident[WhichCells(object = subset(CFA2_all,subset = sample =="CFA_2"))])
x <- CFA2_all$center.x[WhichCells(object = subset(CFA2_all,subset = sample =="CFA_2"))]
y <- CFA2_all$center.y[WhichCells(object = subset(CFA2_all,subset = sample =="CFA_2"))]
group <- CFA2_all@active.ident[WhichCells(object = subset(CFA2_all,subset = sample =="CFA_2"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#CFA_3

CFA3.x <- CFA3@meta.data[,5]
head(CFA3.x)
CFA3.y <- CFA3@meta.data[,6]
head(CFA3.y)

CFA3_cx_x <- c(5955.11,	3898.2,	1431.51,	730.506,	106.451,	77.8216,	274.282,	1504.06,	2429.17,	1995.18,	1806.67,	1793.223,	1380.89,	1211.77,	1439.73,	1758.68,	2368.75,	3230.91,	4172.87,	5216.58,	5528.08,	5693.06,	6216.4,	6764.22,	7095.36,	7826.41,	8398.26,	8556.26,	8337.56,	8088.42,	8266.81,	7485.63,	5459.78,	4911.28,	6404.54,	8614.14,	9950.53,	9965.78,	9139.31,	7435.17,	6006.66)
CFA3_cx_y <- c(1076.3,	525.246,	414.32,	1654.62,	2863.55,	3998.14,	5242.51,	6264.8,	5779.27,	4180.72,	3053.17,	2393.55,	3052.74,	3261.79,	2629.94,	2243.3,	1697.71,	1385.12,	1399.81,	1781.67,	1802.53,	2297,	2601.65,	2528.43,	2828.84,	3373.96,	4221.48,	5471.15,	6304.34,	6666.35,	5980.01,	6457.79,	6887.61,	7446.87,	8650.86,	8093.01,	6284.38,	4421.52,	2872.13,	1844.96,	1116.02)
CFA3_hp_x <- c(5569.96,	5052.86,	4566.43,	3386.33,	2806.37,	2459.61,	2632.97,	3116.92,	3634.96,	4272.57,	4859.71,	5254.3,	5851.03,	5866.66,	6580.13,	7063.16,	7433.5,	7749.58,	8045.13,	8167.44,	8162.33,	7743.9,	7203.54,	6840.51,	6502.53,	5874.89)
CFA3_hp_y <- c(2830.37,	1955.43,	1663.79,	1582.72,	1819.27,	2255.03,	2475.81,	2403.76,	2384.14,	2701.31,	3089.5,	3075.1,	3388.68,	3727.14,	3987.74,	4266.91,	4578.94,	5093.73,	4997.55,	4736.72,	4442.47,	3503.21,	3047.64,	2906.54,	2963.31,	2985.4)
CFA3_wm_x <- c(5747.68,	5409.44,	5093.43,	4199.84,	3214.76,	2700.12,	2028.28,	1351.06,	1792,	1837.11,	2209.91,	2822.03,	2481.79,	2444.07,	2695.88,	3303.98,	4188.09,	4969.21,	5553.64,	5984.12,	6585.99,	7034.95,	7693.25,	8129.78,	8210.43,	7999.94,	7811.82,	7589.17,	7230.35,	7617.63,	8306.75,	8265.91,	7952.8,	8228.6,	8457.26,	8506.48,	8466.39,	8286.36,	7679.44,	6795.5,	6194.41)
CFA3_wm_y <- c(2376.47,	1787.96,	1756.23,	1433.59,	1444.23,	1508.43,	1987.68,	2946.43,	2319.72,	3524.17,	3470.87,	2559.69,	2464.92,	2246.39,	1824.16,	1543.13,	1535.05,	1829.57,	2799.92,	2984.69,	2934.04,	2931.81,	3426.56,	4241.38,	4741.66,	5130.47,	5136.45,	4870.44,	6101.19,	6271.39,	5674.11,	6171.93,	6921.47,	6363.13,	5837.59,	5170.14,	4822.76,	4159.5,	3256.04,	2602.45,	2618.93)
CFA3_th_x <- c(5196.88,	4999.94,	4242.08,	3515.11,	2860.45,	2530.68,	2342.67,	2860.06,	2635,	2784.25,	3756.47,	5120.1,	5290.18,	6956.6,	7586.12,	7427.92,	7033.86,	5808.79,	5225.12)
CFA3_th_y <- c(3856.46,	3209.53,	2707.55,	2411.86,	2620.16,	3001.78,	3723.26,	5229.75,	5676.96,	6520.29,	7845.75,	7032.31,	6490.16,	6159.17,	4961.47,	4593.04,	4254.37,	3711.64,	3865.97)


CFA3_pol_cx <- point.in.polygon(CFA3.x, CFA3.y, CFA3_cx_x, CFA3_cx_y)
CFA3_pol_hp <- point.in.polygon(CFA3.x, CFA3.y, CFA3_hp_x, CFA3_hp_y)
CFA3_pol_wm <- point.in.polygon(CFA3.x, CFA3.y, CFA3_wm_x, CFA3_wm_y)
CFA3_pol_th <- point.in.polygon(CFA3.x, CFA3.y, CFA3_th_x, CFA3_th_y)

CFA3.info <- cbind.data.frame(CFA3_pol_cx,CFA3_pol_hp,CFA3_pol_wm,CFA3_pol_th)
head(CFA3.info, 50)



CFA3.info <- data.frame(CFA3.info)
CFA3@meta.data <- cbind(CFA3@meta.data,CFA3.info)


CFA3_cx <- subset(CFA3,subset = CFA3_pol_cx == 1)
CFA3_hp <- subset(CFA3,subset = CFA3_pol_hp == 1)
CFA3_wm <- subset(CFA3,subset = CFA3_pol_wm == 1)
CFA3_th <- subset(CFA3,subset = CFA3_pol_th == 1)

CFA3_cx$region <- "Cx"
CFA3_hp$region <- "Hp"
CFA3_wm$region <- "Wm"
CFA3_th$region <- "Th"

CFA3_all <- merge(CFA3_cx, y = c( CFA3_hp, CFA3_wm, CFA3_th))


Idents(CFA3_all)<-CFA3_all$region
levels(CFA3_all)

dfall<- data.frame(CFA3_all$center.x[WhichCells(object = subset(CFA3_all,subset = sample =="CFA_3"))],CFA3_all$center.y[WhichCells(object = subset(CFA3_all,subset = sample =="CFA_3"))],CFA3_all@active.ident[WhichCells(object = subset(CFA3_all,subset = sample =="CFA_3"))])
x <- CFA3_all$center.x[WhichCells(object = subset(CFA3_all,subset = sample =="CFA_3"))]
y <- CFA3_all$center.y[WhichCells(object = subset(CFA3_all,subset = sample =="CFA_3"))]
group <- CFA3_all@active.ident[WhichCells(object = subset(CFA3_all,subset = sample =="CFA_3"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()


#EAE_1

EAE1.x <- EAE1@meta.data[,5]
head(EAE1.x)
EAE1.y <- EAE1@meta.data[,6]
head(EAE1.y)

EAE1_cx_x <- c(4518.89,	6362.06,	8031.43,	9588.26,	10237.3,	10272.9,	8947.71,	6925.37,	6951.96,	7694.15,	8487.8,	8503.24,	8953.6,	8841.56,	8303.07,	7409.55,	6401.93,	5309.52,	5141.6,	4896.37,	4400.89,	3894.6,	2857.08,	2186.5,	1636.27,	1266.74,	1421.36,	1460.51,	2053.74,	2927.1,	4225.01,	2934.79,	2312.84,	1411.49,	124.683,	179.357,	366.513,	1953.38,	2974.35,	4213.99)
EAE1_cx_y <- c(7496.53,	7520.03,	7325.79,	6321.03,	5257.39,	2881.02,	1636.72,	1665,	2269.67,	2872.86,	4081.34,	5203.86,	4397.87,	4888.89,	5623.51,	6319.65,	6600.88,	6649.34,	6297.26,	6159.81,	6144.37,	6348.48,	5943.38,	5507.4,	4726.69,	3701.4,	2676.98,	3631.75,	2822.28,	2253.57,	1567.18,	549.845,	462.413,	1437.17,	2625.49,	4838.32,	5173.32,	7025.98,	7428.33,	7429.34)
EAE1_hp_x <- c(4995.69,	5583.98,	6224.49,	6869.85,	7709.55,	7746.88,	7214.42,	6598,	5470.39,	5205.23,	4515.51,	4420.51,	3157.14,	2712.99,	2098.11,	1865.61,	1960.12,	2547.42,	3145.09,	3707,	4659.47)
EAE1_hp_y <- c(5765.73,	6306.56,	6452.18,	6394.84,	5760.66,	5387.69,	5408.06,	5416.15,	5178.29,	5369.69,	5228.28,	5006.1,	4813.97,	4659.43,	4090.19,	4207.05,	4769,	5558.05,	5891.23,	6014.08,	5692.31)
EAE1_wm_x <- c(5012.32,	5409.89,	6557.39,	7479.79,	8469.11,	8954.95,	9077.45,	8416.16,	8510.57,	8000.03,	7494.86,	7353.87,	7753.15,	7771.26,	7358.72,	6568.4,	5832.08,	5453.64,	5009.59,	4556.79,	3831.95,	3332.08,	2586.41,	2032.77,	1845.74,	1930.85,	2413.07,	2360.33,	2400.78,	1712.86,	1542.98,	1333.05,	1408.34,	1272.65,	1339.16,	1592.26,	1980.72,	2668.88,	4005.12,	4521.06)
EAE1_wm_y <- c(6202.48,	6588.34,	6591.1,	6270.78,	5444.45,	4691.49,	3893.56,	5270.37,	4009.53,	4091.38,	5230.13,	5365.64,	5314.44,	5699.09,	6152.68,	6473.14,	6425.62,	6171.53,	5792,	5727.89,	6022,	5989.21,	5613.53,	5007.4,	4374.53,	4098.74,	4312.85,	3645.08,	2894.62,	3075.76,	3929.1,	3435.2,	2940.12,	3451.73,	3937.2,	4553.89,	5297.25,	5828.06,	6274.11,	6113.3)
EAE1_th_x <- c(5114.98,	5429.98,	6605.83,	7111.81,	7567.21,	7990.84,	6786.13,	6741.79,	5852.87,	4734.57,	4050.76,	4145.11,	2564.15,	2438.29,	2361.39,	2511.16,	2707.65,	4474.92,	4922.65)
EAE1_th_y <- c(4780.82,	5171.9,	5406.29,	5405.05,	5074.33,	4017.09,	2623.01,	1563.14,	1204.43,	1164.22,	1874.17,	2127.96,	3061.84,	3486,	3708.36,	4356.95,	4622.71,	4990.57,	4728.61)


EAE1_pol_cx <- point.in.polygon(EAE1.x, EAE1.y, EAE1_cx_x, EAE1_cx_y)
EAE1_pol_hp <- point.in.polygon(EAE1.x, EAE1.y, EAE1_hp_x, EAE1_hp_y)
EAE1_pol_wm <- point.in.polygon(EAE1.x, EAE1.y, EAE1_wm_x, EAE1_wm_y)
EAE1_pol_th <- point.in.polygon(EAE1.x, EAE1.y, EAE1_th_x, EAE1_th_y)

EAE1.info <- cbind.data.frame(EAE1_pol_cx,EAE1_pol_hp,EAE1_pol_wm,EAE1_pol_th)
head(EAE1.info, 50)



EAE1.info <- data.frame(EAE1.info)
EAE1@meta.data <- cbind(EAE1@meta.data,EAE1.info)


EAE1_cx <- subset(EAE1,subset = EAE1_pol_cx == 1)
EAE1_hp <- subset(EAE1,subset = EAE1_pol_hp == 1)
EAE1_wm <- subset(EAE1,subset = EAE1_pol_wm == 1)
EAE1_th <- subset(EAE1,subset = EAE1_pol_th == 1)

EAE1_cx$region <- "Cx"
EAE1_hp$region <- "Hp"
EAE1_wm$region <- "Wm"
EAE1_th$region <- "Th"

EAE1_all <- merge(EAE1_cx, y = c( EAE1_hp, EAE1_wm, EAE1_th))


Idents(EAE1_all)<-EAE1_all$region
levels(EAE1_all)

dfall<- data.frame(EAE1_all$center.x[WhichCells(object = subset(EAE1_all,subset = sample =="EAE_1"))],EAE1_all$center.y[WhichCells(object = subset(EAE1_all,subset = sample =="EAE_1"))],EAE1_all@active.ident[WhichCells(object = subset(EAE1_all,subset = sample =="EAE_1"))])
x <- EAE1_all$center.x[WhichCells(object = subset(EAE1_all,subset = sample =="EAE_1"))]
y <- EAE1_all$center.y[WhichCells(object = subset(EAE1_all,subset = sample =="EAE_1"))]
group <- EAE1_all@active.ident[WhichCells(object = subset(EAE1_all,subset = sample =="EAE_1"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()

#EAE_2

EAE2.x <- EAE2@meta.data[,5]
head(EAE2.x)
EAE2.y <- EAE2@meta.data[,6]
head(EAE2.y)

EAE2_cx_x <- c(390.082,	965.789,	2380.28,	4138.25,	6505.82,	7818.56,	6377.77,	5542.11,	4968.91,	3879.66,	4528.37,	5430.94,	4619.13,	3695.12,	2554.26,	1718.59,	1322.28,	1639.48,	1577.48,	1029.69,	947.058,	999.563,	1373.43,	2138.7,	2889.53,	3590.02,	2852.6,	3322.15,	5708.13,	6177.69,	3583.58,	2653.16,	1309.92,	823.68,	93.237,	38.103,	83.1406,	136.512)
EAE2_cx_y <- c(6980.85,	8499,	10001.3,	10547.6,	9402.24,	8247.01,	6382.88,	7477.57,	8301.28,	8892.86,	9037.26,	8858.98,	9212.24,	9128.94,	8548.22,	7536.86,	6777.29,	6288.53,	5711.54,	5295.1,	4365.94,	3539.57,	2654.13,	1881.25,	1379.5,	1249.87,	1687.44,	1967.45,	3399.98,	852.591,	-2.69029,	191.177,	1049.82,	1350.47,	2992.66,	4005.24,	4906.23,	5451.72)
EAE2_hp_x <- c(2158.8,	1770.07,	1980,	2851.17,	3058.36,	3460.92,	3029.96,	2938.84,	2914.07,	2729.96,	2538.72,	2440.21,	2655.05,	2210.45,	2128.49,	2160.19,	2273.11,	2041.61,	1564.24,	1227.47,	1066.34,	1254.05,	2097.88)
EAE2_hp_y <- c(5872.96,	6821.45,	7640.24,	8386.8,	8460.1,	8287.6,	7320.26,	6804.67,	6156.82,	5966.89,	5887.08,	5482.76,	5013.93,	4004.56,	3492.01,	3213.52,	2688.7,	2628.23,	2882.02,	3539.98,	4330.85,	5010.37,	5685.51)
EAE2_wm_x <- c(1690.84,	1477.01,	1839.51,	2277.79,	2849.48,	3675.2,	4332.51,	4825.11,	4334.08,	3647.85,	4439.42,	3344.55,	3489.19,	3056.51,	2335.2,	1831.74,	1740.06,	2104.57,	2052.39,	1281.11,	1047.42,	1251.84,	1632.35,	2183.01,	2301.07,	2171,	3487.62,	2196.26,	2701.663,	3395.66,	2647.31,	2194.89,	1532.15,	1145.2,	970.811,	982.297,	1059.44,	1584.66)
EAE2_wm_y <- c(6080.71,	6759.29,	7706.49,	8203.79,	8624.67,	8998.48,	9138.44,	99047.39,	9069.7,	8891.47,	8182.61,	7933.95,	8340.35,	8496.82,	8051.43,	7384.34,	6726.1,	5918.64,	5654.16,	5091.54,	4369.06,	3447.08,	2773.81,	2571.29,	2768.36,	3210.37,	2451.98,	2090.47,	1615.11,	1303.29,	1549.71,	1881.71,	2543.7,	3177.73,	3867.64,	4718.15,	5278.48,	5720.7)
EAE2_th_x <- c(3128.49,	2813.85,	2931.36,	2973.39,	3045.7,	3320.06,	4818.67,	5609.06,	6008.44,	6504.91,	6582.74,	6595.04,	6041.3,	5086.9,	4402.27,	3590.34,	2255.96,	2160.45,	2222.43,	2665.83,	3039.42)
EAE2_th_y <- c(5551.49,	6034.64,	6182.26,	6921.03,	7320.64,	7876.88,	7890.44,	6279.83,	6274.73,	5658.7,	5412.96,	4894.01,	3294.45,	3654.71,	3147.9,	2548.07,	3204.89,	3562.64,	3988.86,	5032.58,	5483.5)


EAE2_pol_cx <- point.in.polygon(EAE2.x, EAE2.y, EAE2_cx_x, EAE2_cx_y)
EAE2_pol_hp <- point.in.polygon(EAE2.x, EAE2.y, EAE2_hp_x, EAE2_hp_y)
EAE2_pol_wm <- point.in.polygon(EAE2.x, EAE2.y, EAE2_wm_x, EAE2_wm_y)
EAE2_pol_th <- point.in.polygon(EAE2.x, EAE2.y, EAE2_th_x, EAE2_th_y)

EAE2.info <- cbind.data.frame(EAE2_pol_cx,EAE2_pol_hp,EAE2_pol_wm,EAE2_pol_th)
head(EAE2.info, 50)



EAE2.info <- data.frame(EAE2.info)
EAE2@meta.data <- cbind(EAE2@meta.data,EAE2.info)


EAE2_cx <- subset(EAE2,subset = EAE2_pol_cx == 1)
EAE2_hp <- subset(EAE2,subset = EAE2_pol_hp == 1)
EAE2_wm <- subset(EAE2,subset = EAE2_pol_wm == 1)
EAE2_th <- subset(EAE2,subset = EAE2_pol_th == 1)

EAE2_cx$region <- "Cx"
EAE2_hp$region <- "Hp"
EAE2_wm$region <- "Wm"
EAE2_th$region <- "Th"

EAE2_all <- merge(EAE2_cx, y = c( EAE2_wm, EAE2_hp, EAE2_th))


Idents(EAE2_all)<-EAE2_all$region
levels(EAE2_all)

dfall<- data.frame(EAE2_all$center.x[WhichCells(object = subset(EAE2_all,subset = sample =="EAE_2"))],EAE2_all$center.y[WhichCells(object = subset(EAE2_all,subset = sample =="EAE_2"))],EAE2_all@active.ident[WhichCells(object = subset(EAE2_all,subset = sample =="EAE_2"))])
x <- EAE2_all$center.x[WhichCells(object = subset(EAE2_all,subset = sample =="EAE_2"))]
y <- EAE2_all$center.y[WhichCells(object = subset(EAE2_all,subset = sample =="EAE_2"))]
group <- EAE2_all@active.ident[WhichCells(object = subset(EAE2_all,subset = sample =="EAE_2"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()


#EAE_3

EAE3.x <- EAE3@meta.data[,5]
head(EAE3.x)
EAE3.y <- EAE3@meta.data[,6]
head(EAE3.y)

EAE3_cx_x <- c(723.503,	61.2397,	120.568,	377.209,	1985.44,	2769.11,	4443.64,	4813.7,	5306.06,	3437.26,	2623.46,	1696.85,	1920.64,	2294.95,	1718.94,	1364.03,	1214.21,	1412.08,	1898.67,	2508.58,	2879.21,	2889.52,	3594.54,	4512.14,	5574.83,	6572.93,	7501.87,	6360.55,	6730.28,	7157.51,	7038.05,	7459.34,	8447.8,	9032.81,	8898.86,	7207.95,	6719.83,	5440.29,	4032.35,	1868.16)
EAE3_cx_y <- c(3270.52,	5378.61,	6763.59,	7463.98,	9423.83,	9612.95,	9392.81,	9184.32,	7322.19,	7649.81,	7334.55,	6894.48,	7521.19,	7909.9,	7319.17,	6546.04,	5393.57,	4530.04,	3659.67,	3526.25,	3157.65,	2598.56,	1903.94,	1367.28,	1250.47,	1471.44,	2190,	1664.16,	2355.42,	3827.14,	4930.68,	5188.81,	4474.15,	3195.86,	2042.31,	688.22,	404.925,	116.728,	141.19,	1932.45)
EAE3_hp_x <- c(3018.25,	1907.14,	1414.79,	1453.55,	1794.29,	2052.62,	2197.2,	2709.93,	3327.04,	3182.95,	3649.34,	3972.09,	4415.75,	4963.16,	5608.71,	5632.86,	5181.95,	4346.31,	3466.64,	3283.31,	3230.14)
EAE3_hp_y <- c(3704.74,	4077.85,	4933.26,	6077.53,	6413.93,	6277.72,	5639.84,	4912.15,	4406.89,	4051.68,	3611.12,	3647.54,	2906.3,	2341.63,	2028.34,	1770.27,	1528.24,	1656.16,	2189.26,	2539.07,	3496.2)
EAE3_wm_x <- c(2650.35,	1970.13,	1372.47,	1203.52,	1473.41,	2025.4,	1595.5,	2389.09,	2671.66,	212282,	1966.81,	1732.43,	1451.27,	1342.48,	1567.68,	1952.9,	2939.99,	3210.7,	3321.07,	4033.23,	5200.79,	5648.95,	5700.6,	5300.84,	6580.65,	6378.75,	5849.49,	6632.91,	7093.04,	6661.02,	6006.07,	5238.19,	4288.74,	3554.2,	2864.35,	2929.95,	2789.61)
EAE3_wm_y <- c(3470.47,	3707.64,	4634.44,	5578.74,	6771.16,	7706.7,	6884.77,	7270.79,	6873.05,	6024.4,	6392.54,	6432.1,	6169.54,	5240.92,	4526.69,	3996.98,	3727.82,	3391.79,	2388.06,	1729.63,	1507.97,	1740.85,	1945.61,	2188.53,	2707.41,	1792.27,	1540.87,	1635.53,	1911.44,	1549.41,	1337.6,	1322.58,	1500.26,	1944,	2566.63,	3035.59,	3308.22)
EAE3_th_x <- c(3984.24,	3328.9,	2363.79,	2153.81,	2209.44,	2769.75,	3167.36,	4841.94,	5282.25,	6392.61,	7434.04,	6645.29,	6604.48,	6170.19,	5466.96,	5110.74,	4849.54,	3964.66,	3966.05)
EAE3_th_y <- c(4353.86,	4369.99,	5359.04,	6008.43,	6190.78,	6959.31,	7135.3,	6856.66,	7341.03,	7197.13,	5989.12,	4803.62,	3071.92,	2522.67,	2193.97,	2283.25,	2470.67,	3693.67,	4228.38)


EAE3_pol_cx <- point.in.polygon(EAE3.x, EAE3.y, EAE3_cx_x, EAE3_cx_y)
EAE3_pol_hp <- point.in.polygon(EAE3.x, EAE3.y, EAE3_hp_x, EAE3_hp_y)
EAE3_pol_wm <- point.in.polygon(EAE3.x, EAE3.y, EAE3_wm_x, EAE3_wm_y)
EAE3_pol_th <- point.in.polygon(EAE3.x, EAE3.y, EAE3_th_x, EAE3_th_y)

EAE3.info <- cbind.data.frame(EAE3_pol_cx,EAE3_pol_hp,EAE3_pol_wm,EAE3_pol_th)
head(EAE3.info, 50)



EAE3.info <- data.frame(EAE3.info)
EAE3@meta.data <- cbind(EAE3@meta.data,EAE3.info)


EAE3_cx <- subset(EAE3,subset = EAE3_pol_cx == 1)
EAE3_hp <- subset(EAE3,subset = EAE3_pol_hp == 1)
EAE3_wm <- subset(EAE3,subset = EAE3_pol_wm == 1)
EAE3_th <- subset(EAE3,subset = EAE3_pol_th == 1)

EAE3_cx$region <- "Cx"
EAE3_hp$region <- "Hp"
EAE3_wm$region <- "Wm"
EAE3_th$region <- "Th"

EAE3_all <- merge(EAE3_cx, y = c( EAE3_wm, EAE3_hp, EAE3_th))


Idents(EAE3_all)<-EAE3_all$region
levels(EAE3_all)

dfall<- data.frame(EAE3_all$center.x[WhichCells(object = subset(EAE3_all,subset = sample =="EAE_3"))],EAE3_all$center.y[WhichCells(object = subset(EAE3_all,subset = sample =="EAE_3"))],EAE3_all@active.ident[WhichCells(object = subset(EAE3_all,subset = sample =="EAE_3"))])
x <- EAE3_all$center.x[WhichCells(object = subset(EAE3_all,subset = sample =="EAE_3"))]
y <- EAE3_all$center.y[WhichCells(object = subset(EAE3_all,subset = sample =="EAE_3"))]
group <- EAE3_all@active.ident[WhichCells(object = subset(EAE3_all,subset = sample =="EAE_3"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=.1)+  theme_classic()+coord_fixed()+NoLegend()




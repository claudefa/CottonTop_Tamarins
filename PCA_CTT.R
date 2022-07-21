library("ggplot2")
library("tidyr")
library("ggConvexHull")
library("scales")
library("ggrepel")
library("RColorBrewer")
library("viridis")
library("readxl")

setwd("~/Documents/OneDrive - upf.edu/Doctorat/Cotton-top/")
setwd("~/Documents/Feina/Cotton-top/")

metadata <- read_excel("Shotgun_Metadata.xlsx")

# PCA with ANGSD with all samples 1st Scaffold --------
names <- read.table("Files/Samples")
m <- as.matrix(read.table("Files/CTT_CM038391.1_covmatrix.cov"))

#get the eigen vectors
e <- eigen(m)

#convert to ggplot dataframe
eigenvectors <- data.frame(e$vectors)
eigenvalues <- data.frame(e$values)

row.names(eigenvectors) <- names$V1
eigenvectors$Full_ID <- names$V1


eigenvectors <- merge(eigenvectors, metadata[ , 1:12], by="Full_ID")                            

# Eigenvalues Percentage each PC
# PC1 and PC2
eigenvalues$e.values[1]/sum(eigenvalues$e.values)
eigenvalues$e.values[2]/sum(eigenvalues$e.values)
eigenvalues$e.values[3]/sum(eigenvalues$e.values)
eigenvalues$e.values[4]/sum(eigenvalues$e.values)

centroids <- aggregate(cbind(X1,X2)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X1","X2") 

# Plot PC1 and PC2
g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X1, X2,  color=Site, shape=Type, size=Coverage)) + 
  ylab("PC2 (7.91%)") + xlab("PC1 (14.17%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 - 787,262 SNPs)") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X1, X2, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X1, X2, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(X3,X4)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X3","X4") 

g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X3, X4,  color=Site, shape=Type), size=4) + ylab("PC4 (3.68%)") + xlab("PC3 (4.08%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 - 787,262 SNPs)") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X3, X4, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X3, X4, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_pc3pc4.pdf", height=8, width=15)
g
dev.off()



# PCA with ANGSD with all samples 1st Scaffold MAF 0.05--------
names <- read.table("Files/Samples")
m <- as.matrix(read.table("Files/CTT_CM038391.1_maf0.05_covmatrix.cov"))

#get the eigen vectors
e <- eigen(m)

#convert to ggplot dataframe
eigenvectors <- data.frame(e$vectors)
eigenvalues <- data.frame(e$values)

row.names(eigenvectors) <- names$V1
eigenvectors$Full_ID <- names$V1


eigenvectors <- merge(eigenvectors, metadata[ , 1:12], by="Full_ID")                            

# Eigenvalues Percentage each PC
# PC1 and PC2
eigenvalues$e.values[1]/sum(eigenvalues$e.values)
eigenvalues$e.values[2]/sum(eigenvalues$e.values)
eigenvalues$e.values[3]/sum(eigenvalues$e.values)
eigenvalues$e.values[4]/sum(eigenvalues$e.values)

centroids <- aggregate(cbind(X1,X2)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X1","X2") 

# Plot PC1 and PC2
g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X1, X2,  color=Site, shape=Type), size=4) + 
  ylab("PC2 (7.90%)") + xlab("PC1 (13.82%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 - 500,285 SNPs) MAF 0.05") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X1, X2, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X1, X2, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_MAF0.05.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(X3,X4)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X3","X4") 

g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X3, X4,  color=Site, shape=Type), size=4) + ylab("PC4 (3.64%)") + xlab("PC3 (4.11%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 - 500,285 SNPs) MAF 0.05") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X3, X4, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X3, X4, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_pc3pc4_MAF0.05.pdf", height=8, width=15)
g
dev.off()




# PCA with ANGSD WITHOUT transitions --------
names <- read.table("Files/Samples")
m <- as.matrix(read.table("Files/CTT_CM038391.1_rmTrans_covmatrix.cov"))

#get the eigen vectors
e <- eigen(m)

#convert to ggplot dataframe
eigenvectors <- data.frame(e$vectors)
eigenvalues <- data.frame(e$values)
row.names(eigenvectors) <- names$V1
eigenvectors$Full_ID <- names$V1
eigenvectors <- merge(eigenvectors, metadata[ , 1:12], by="Full_ID")                            


# Eigenvalues Percentage each PC
# PC1 and PC2
eigenvalues$e.values[1]/sum(eigenvalues$e.values)
eigenvalues$e.values[2]/sum(eigenvalues$e.values)
eigenvalues$e.values[3]/sum(eigenvalues$e.values)
eigenvalues$e.values[4]/sum(eigenvalues$e.values)

centroids <- aggregate(cbind(X1,X2)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X1","X2") 

# Plot PC1 and PC2
g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X1, X2,  color=Site, shape=Type), size=4) + ylab("PC2 (8.68%)") + xlab("PC1 (16.24%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 - 196,585 SNPs) - NO TRANSITIONS") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X1, X2, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X1, X2, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_rmTrans.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(X3,X4)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X3","X4") 

g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X3, X4,  color=Site, shape=Type), size=4) + ylab("PC4 (3.85%)") + xlab("PC3 (4.19%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 -  196,585 SNPs) - NO TRANSITIONS") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X3, X4, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X3, X4, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_pc3pc4_rmTrans.pdf", height=8, width=15)
g
dev.off()


# PCA with ANGSD WITHOUT TRIM 10BP --------
names <- read.table("Files/Samples")
m <- as.matrix(read.table("Files/CTT_CM038391.1_trim_covmatrix.cov"))

#get the eigen vectors
e <- eigen(m)

#convert to ggplot dataframe
eigenvectors <- data.frame(e$vectors)
eigenvalues <- data.frame(e$values)
row.names(eigenvectors) <- names$V1
eigenvectors$Full_ID <- names$V1
eigenvectors <- merge(eigenvectors, metadata[ , 1:12], by="Full_ID")                            


# Eigenvalues Percentage each PC
# PC1 and PC2
eigenvalues$e.values[1]/sum(eigenvalues$e.values)
eigenvalues$e.values[2]/sum(eigenvalues$e.values)
eigenvalues$e.values[3]/sum(eigenvalues$e.values)
eigenvalues$e.values[4]/sum(eigenvalues$e.values)

centroids <- aggregate(cbind(X1,X2)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X1","X2") 

# Plot PC1 and PC2
g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X1, X2,  color=Site, shape=Type), size=4) + ylab("PC2 (4,94%)") + xlab("PC1 (14.43%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 -  667,814 SNPs) - Trim 10bp") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X1, X2, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X1, X2, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_trim10bp.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(X3,X4)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X3","X4") 

g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X3, X4,  color=Site, shape=Type), size=4) + ylab("PC4 (3.68%)") + xlab("PC3 (4.08%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 - 667,814 SNPs) - Trim 10bp") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X3, X4, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X3, X4, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_pc3pc4_trim10bp.pdf", height=8, width=15)
g
dev.off()


# PCA from snpAD + plink PROJECTED

pca_plink <- read.table("Files/CTT_CM038391.1_projected.eigenvec")
val_plink <- read.table("Files/CTT_CM038391.1_projected.eigenval")
val_plink[1,1]/sum(val_plink$V1) *100
val_plink[2,1]/sum(val_plink$V1) *100
val_plink[3,1]/sum(val_plink$V1) *100
val_plink[4,1]/sum(val_plink$V1) *100

pca_plink$ID_vcf <- pca_plink$V1
pca_plink <- merge(pca_plink, metadata[ , 1:12], by="ID_vcf")                            
centroids <- aggregate(cbind(V3,V4)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V3","V4") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V3, V4,  color=Site, shape=Type), size=4) + ylab("PC2 (8.64%)") + xlab("PC1 (17.94%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 -  115,559 SNPs) - VCF") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V3, V4,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V3, V4,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_VCF.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(V5,V6)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V5","V6") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V5, V6,  color=Site, shape=Type), size=4) + ylab("PC4 (6.50%)") + xlab("PC3 (6.74%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 -  115,559 SNPs) - VCF") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V5, V6,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V5, V6,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g
pdf("Plots/PCA_allsamples_allsites_pc3pc4_VCF.pdf", height=8, width=15)
g
dev.off()


# PCA from snpAD + plink

pca_plink <- read.table("Files/CTT_CM038391.1.eigenvec")
val_plink <- read.table("Files/CTT_CM038391.1.eigenval")
val_plink[1,1]/sum(val_plink$V1) *100
val_plink[2,1]/sum(val_plink$V1) *100
val_plink[3,1]/sum(val_plink$V1) *100
val_plink[4,1]/sum(val_plink$V1) *100

pca_plink$ID_vcf <- pca_plink$V1
pca_plink <- merge(pca_plink, metadata[ , 1:12], by="ID_vcf")                            
centroids <- aggregate(cbind(V3,V4)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V3","V4") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V3, V4,  color=Site, shape=Type), size=4) + ylab("PC2 (8.64%)") + xlab("PC1 (17.94%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 -  115,559 SNPs) - VCF") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V3, V4,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V3, V4,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g
pdf("Plots/PCA_allsamples_allsites_VCF.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(V5,V6)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V5","V6") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V5, V6,  color=Site, shape=Type), size=4) + ylab("PC4 (6.50%)") + xlab("PC3 (6.74%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 -  115,559 SNPs) - VCF") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V5, V6,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V5, V6,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g
pdf("Plots/PCA_allsamples_allsites_pc3pc4_VCF.pdf", height=8, width=15)
g
dev.off()




# PCA from snpAD + plink PROJECTING MODERN INTO HISTORICAL

pca_plink <- read.table("Files/CTT_CM038391.1_projected.eigenvec")
val_plink <- read.table("Files/CTT_CM038391.1_projected.eigenval")
val_plink[1,1]/sum(val_plink$V1) *100
val_plink[2,1]/sum(val_plink$V1) *100
val_plink[3,1]/sum(val_plink$V1) *100
val_plink[4,1]/sum(val_plink$V1) *100


pca_plink$ID_vcf <- pca_plink$V1
pca_plink <- merge(pca_plink, metadata[ , 1:12], by="ID_vcf")                            
centroids <- aggregate(cbind(V3,V4)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V3","V4") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V3, V4,  color=Site, shape=Type), size=4) + ylab("PC2 (8.75%)") + xlab("PC1 (14.73%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 -   SNPs) - VCF") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V3, V4,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V3, V4,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g
pdf("Plots/PCA_allsamples_allsites_VCF_projected.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(V5,V6)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V5","V6") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V5, V6,  color=Site, shape=Type), size=4) + ylab("PC4 (6.63%)") + xlab("PC3 (7.49%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples (CM038391.1 -  115,559 SNPs) - VCF") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V5, V6,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V5, V6,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g
pdf("Plots/PCA_allsamples_allsites_pc3pc4_VCF_projected.pdf", height=8, width=15)
g
dev.off()



# F INBREEDING -----------

g <- ggplot(metadata, aes(Site, F, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() +
  theme(legend.position = "none")
g

pdf("Plots/Inbreeding_F.pdf", height=4, width=11)
g
dev.off()

g <- ggplot(metadata[metadata$Coverage>5,], aes(Site, F, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() +ggtitle("Only Samples > 5x")+
  theme(legend.position = "none")
g

pdf("Plots/Inbreeding_F_highCov.pdf", height=4, width=11)
g
dev.off()

g <- ggplot(metadata, aes(Coverage, F, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() 
g

pdf("Plots/Inbreeding_FvsCoverage.pdf", height=5, width=6)
g
dev.off()

g <- ggplot(metadata, aes(Missingess, F, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() + xlab("Missingness")
g

pdf("Plots/Inbreeding_FvsMissingness.pdf", height=5, width=6)
g
dev.off()


g <- ggplot(metadata, aes(Missingess, Coverage, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() + xlab("Missingness")
g

pdf("Plots/CoveragevsMissingness.pdf", height=5, width=6)
g
dev.off()


# Heterozygosity
g <- ggplot(metadata, aes(Site, Heterozygosity, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() +
  theme(legend.position = "none")
g

pdf("Plots/Heterozygosity.pdf", height=4, width=11)
g
dev.off()

g <- ggplot(metadata[metadata$Coverage>5,], aes(Site, Heterozygosity, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() + ggtitle("Only Samples > 5x")+
  theme(legend.position = "none")
g

pdf("Plots/Heterozygosity_highCov.pdf", height=4, width=11)
g
dev.off()

g <- ggplot(metadata, aes(Missingess, Heterozygosity, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() + xlab("Missingness")
g

pdf("Plots/HeterozygosityvsMissingness.pdf", height=5, width=6)
g
dev.off()

g <- ggplot(metadata, aes(Coverage, Heterozygosity, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() 
g

pdf("Plots/HeterozygosityvsCoverage.pdf", height=5, width=6)
g
dev.off()

# RUNS OF HOMOZYGOSITY by Window heterozygosity

# Per window
het <- list.files("Files/Het/", pattern = ".het", full.names = TRUE)
ldf_het <- lapply(het, read.table)

samplesid <- unlist(lapply(1:length(basename(het)), function(i) strsplit(basename(het)[i],"_")[[1]][1]))
         
for (i in 1:length(samplesid)){
  ID_vcf <- samplesid[i]
  ldf_het[[i]]<- cbind.data.frame(ldf_het[[i]], ID_vcf=as.character(ID_vcf))
}

df_het <- do.call(rbind.data.frame,ldf_het)

df_het[df_het$Callable < 60000,]$Heterozygosity <- "NA"
df_het <- df_het[complete.cases(df_het),]

df_het1 <- merge(df_het, metadata[ , 1:12], by="ID_vcf")                            


colnames(df_het1) <- c("ID_vcf","Scaffold", "Pos", "NumHet", "Callable", "Heterozygosity","Sample","ID",
                      "CollectionYear","Type","Material","Provience","Site","Latitude","Longitud","Sex",
                      "Coverage")
df_het <- df_het1
df_het$Sample <- gsub("CEI_051","test", df_het$Sample)
df_het$Sample <- gsub("SJN_001","CEI_051_c", df_het$Sample)
df_het$Sample <- gsub("test","SJN_001_c", df_het$Sample)



# Density of callable size
g <- ggplot(df_het, aes(Callable, color=Sample))+
  geom_density()+ facet_wrap(.~Site)+
  theme(panel.grid.major = element_line(colour = "grey"),
        legend.position = "none",
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

g

# Correlation of callable size heterozygosity
g <- ggplot(df_het, aes(Callable, Heterozygosity, color=Sample))+
  geom_point()+ facet_wrap(.~Sample)+
  theme(panel.grid.major = element_line(colour = "grey"),
        legend.position = "none",
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

g

median(df_het[df_het$ID_vcf=="1.variant",]$Heterozygosity)
mean(df_het[df_het$ID_vcf=="1.variant",]$Heterozygosity)

metadata[metadata$ID_vcf=="1.variant",]$Heterozygosity


# Select more than 1 standard deviation from the mean in each subspecies. 


c <- mean(df_new$Heterozygosity) - sd(df_new$Heterozygosity)

#remove windows with > 2*SD of heterozygosity 
#df_new$Heterozygosity <- as.numeric(df_new$Heterozygosity)

#df_het[df_het$Heterozygosity > 0.004,]$Heterozygosity <- "NA"
df_new <- df_het[complete.cases(df_het),]
df_new$Heterozygosity <- as.numeric(df_new$Heterozygosity)

completeScaffolds <- unique(df_new$Scaffold)


for (i in 1:length(completeScaffolds)){   
  chrom <- completeScaffolds[i]
  df <- df_new[df_new$Scaffold == chrom,]
  g <- ggplot(df, aes(Heterozygosity, color=Sample))+
    geom_density()+ facet_wrap(.~Site)+ xlim(c(0,0.003))+
    ggtitle(paste(chrom)) +
    theme(panel.grid.major = element_line(colour = "grey"),
          legend.position = "none", axis.text.x = element_text(angle=45, hjust = 1),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  name=paste("Plots/Het_windows/Het_",chrom,"_density.pdf", sep="")
  pdf(name, width = 10,height = 6)
  print(g)
  dev.off()
}


for (i in 1:length(completeScaffolds)){   
  chrom <- completeScaffolds[i]
  df <- df_new[df_new$Scaffold == chrom,]
  g <- ggplot(df, aes(Pos,Heterozygosity))
  g <- g + geom_point(size=0.1) + facet_wrap(.~Sample, ncol=1,strip.position = "right")+
    ggtitle(paste(chrom)) + ylim(0,0.0071)+
    geom_point(data=df[df$Heterozygosity<0.00001,], size=0.5, color="red", aes(Pos,Heterozygosity)) +
    scale_x_continuous(name = "Position", waiver(), labels=comma, expand = c(0,0))+
    theme(panel.grid.major = element_line(colour = "grey"),
          legend.position = "none",
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  name=paste("Plots/Het_windows/Het_",chrom,".pdf", sep="")
  pdf(name, width = 6,height = 20)
  print(g)
  dev.off()
}


#FST--------
library(corrplot)
Fst <- read.table("Files/Fst_total", header = FALSE) # without related


colnames(Fst) <- c("Pop1","Pop2","UnWeighted","Weighted")
Fst$Pop1 <- as.character(Fst$Pop1)
Fst$Pop2 <- as.character(Fst$Pop2)
notcomp <- c("Unknown","Cauca","Arjona","Caracas")
Fst6 <- Fst[-which(Fst$Pop1%in%notcomp),]
Fst6 <- Fst6[-which(Fst6$Pop2%in%notcomp),]

Fst6[Fst6$UnWeighted<0,]$UnWeighted <- 0
Fst6[Fst6$Weighted<0,]$Weighted <- 0
Fst_UnWeighted <- Fst6[,1:3,]
Fst_Weighted <- cbind.data.frame(Fst6[,1:2,], Fst6[,4])
colnames(Fst_Weighted) <- c("Pop1","Pop2","Weighted")

# spread
x <- spread(Fst_Weighted, Pop1, Weighted, fill = 0, convert = FALSE)
names <- x$Pop2
x <- x[,2:ncol(x)]
rownames(x) <- names
y <- data.matrix(x)
pdf("Plots/Fst_weighted_noRelated.pdf", height=8, width = 8 )
corrplot(y,addrect = 3, na.label = "NA",
         tl.col = "black", type="lower",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90,
         pch.cex = 1.0, method="color",is.corr = FALSE,
         pch.col = "purple")

dev.off()
pdf("Plots/Fst_weighted_noRelated_fullmatrix.pdf", height=8, width = 8 )
corrplot(y,addrect = 3, na.label = "NA",
         tl.col = "black", type="lower",#order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90,
         pch.cex = 1.0, method="color",is.corr = FALSE,
         pch.col = "purple")

dev.off()

library("ape")
tree <-nj(y)
pdf("Plots/Tree_fst.pdf")
plot(tree, type = "unr", show.tip.lab = TRUE, font=1, cex=0.5)
dev.off()



### Relatedness --------
# Relatedness
library(tidyr)
library(gplots)
library(Matrix)

rel_vcf <- read.table("Files/CTT_rel.ml", header = TRUE)
Names <- read.table("Files/Samples")$V1


Names<- gsub("CEI_051","test", Names)
Names <- gsub("SJN_001","CEI_051_c", Names)
Names <- gsub("test","SJN_001_c", Names)

# For two outbred individuals we write (k0, k1, k2) 
# for the probabilities that they have exactly 0, 1, and 2 alleles IBD. Then θ = k1/4 + k2/2.
rel_vcf$a <- rel_vcf$a +1
rel_vcf$b <- rel_vcf$b +1
mat <- sparseMatrix(i=rel_vcf$a, j=rel_vcf$b, x=rel_vcf$theta, symmetric = TRUE)
diag(mat) <- 0.5
mat <- as.matrix(mat)

colnames(mat) <- Names
row.names(mat) <- Names

x <- mat
colnames(mat)

mate2 <- as.data.frame(mat)
mate2$id1 <- row.names(mate2)
new <-gather(mate2, key, value, -id1)
colnames(new) <- c("ID1","ID2", "Kinship")

length(rel_vcf[rel_vcf$theta>0.046875 & rel_vcf$theta < 0.09276,]$a)
length(rel_vcf[rel_vcf$theta>0.09276 & rel_vcf$theta < 0.1875,]$a)
length(rel_vcf[rel_vcf$theta>0.1875 & rel_vcf$theta < 0.375,]$a)

length(rel_vcf[rel_vcf$theta>0.046875 & rel_vcf$theta < 0.5,]$a)

length(rel_vcf[rel_vcf$theta < 0.5,]$a)



g <- ggplot(new, aes(ID1, ID2, fill=Kinship))
g <- g + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0)) + 
  scale_fill_viridis_c() + ylab("") + xlab("")
g

pdf("Plots/RelatednessGlobal.pdf", height=10, width = 10)
g
dev.off()
pdf("Plots/RelatednessGlobal_corrplot.pdf", height=10, width = 10)
corrplot(x,addrect = 3, na.label = "NA",
         tl.col = "black", type="lower",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90,
         pch.cex = 1.0, method="color",is.corr = FALSE,
         pch.col = "purple")
dev.off()

library("ggplot2")
library("tidyr")
library("ggConvexHull")
library("scales")
library("ggrepel")
library("RColorBrewer")
library("viridis")
library("readxl")

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
metadata <- read_excel("Shotgun_Metadata.xlsx")

# PCA with ANGSD with all samples all scaffolds --------
names <- read.table("Files/Samples_Bams")
m <- as.matrix(read.table("Files/CTT_allsamples_covmatrix.cov"))

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
g <- g +   geom_point(data=eigenvectors, aes(X1, X2,  color=Site, shape=Type)) + 
  ylab("PC2 (5.4%)") + xlab("PC1 (18.41%)") + 
  ggtitle("PCA - Cotton-Top Tamarin N=34 + S. geoffroyi (8,702,060 SNPs)") + theme_test() + 
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
g <- g +   geom_point(data=eigenvectors, aes(X3, X4,  color=Site, shape=Type), size=4) + 
  ylab("PC4 (4.02%)") + xlab("PC3 (4.24%)") + 
  ggtitle("PCA - Cotton-Top Tamarin N=34 + S. geoffroyi (8,702,060 SNPs)") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X3, X4, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X3, X4, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_allsites_pc3pc4.pdf", height=8, width=15)
g
dev.off()

# PCA with ANGSD with only CTT samples nomaf--------
names <- read.table("Files/Samples_Bams_onlyCTT")
m <- as.matrix(read.table("Files/CTT_only_covmatrix.cov"))

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
  ylab("PC2 (5.56%)") + xlab("PC1 (13.06%)") + 
  ggtitle("PCA - Cotton-Top Tamarin N=34 (8,001,210 SNPs)") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X1, X2, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X1, X2, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_onlyCTT_allsites.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(X3,X4)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X3","X4") 

g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X3, X4,  color=Site, shape=Type), size=4) +
  ylab("PC4 (3.87%)") + xlab("PC3 (4.51%)") + 
  ggtitle("PCA - Cotton-Top Tamarin N=34 (8,001,210 SNPs)") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X3, X4, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X3, X4, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_onlyCTT_allsites_pc3pc4.pdf", height=8, width=15)
g
dev.off()

# PCA with ANGSD with only CTT samples maf 0.05 --------
names <- read.table("Files/Samples_Bams_onlyCTT")
m <- as.matrix(read.table("Files/CTT_only_maf_covmatrix.cov"))

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
  ylab("PC2 (5.56%)") + xlab("PC1 (13.04%)") + 
  ggtitle("PCA - Cotton-Top Tamarin N=34 (4,705,720 SNPs)") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X1, X2, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X1, X2, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_onlyCTT_allsites_maf0.05.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(X3,X4)~Site,eigenvectors,mean)
colnames(centroids) <- c("Site","X3","X4") 

g <- ggplot()
g <- g +   geom_point(data=eigenvectors, aes(X3, X4,  color=Site, shape=Type), size=4) +
  ylab("PC4 (3.84%)") + xlab("PC3 (4.50%)") + 
  ggtitle("PCA - Cotton-Top Tamarin N=34 (4,705,720 SNPs)") + theme_test() + 
  geom_convexhull(data=eigenvectors, aes(X3, X4, fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(X3, X4, label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_onlyCTT_allsites_maf0.05_pc3pc4.pdf", height=8, width=15)
g
dev.off()


# PCA from snpAD + plink all autosomes ----------------
pca_plink <- read.table("Files/CTT_filter_onlyCTT.eigenvec")
val_plink <- read.table("Files/CTT_filter_onlyCTT.eigenval")
val_plink[1,1]/sum(val_plink$V1) *100
val_plink[2,1]/sum(val_plink$V1) *100
val_plink[3,1]/sum(val_plink$V1) *100
val_plink[4,1]/sum(val_plink$V1) *100

pca_plink$Full_ID <- pca_plink$V1
pca_plink <- merge(pca_plink, metadata[ , 1:12], by="Full_ID")                            
centroids <- aggregate(cbind(V3,V4)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V3","V4") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V3, V4,  color=Site, shape=Type), size=4) + ylab("PC2 (7.65%)") + xlab("PC1 (12.32%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 34 samples (4,522,615 SNPs)") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V3, V4,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V3, V4,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g
pdf("Plots/PCA_allsamples_allsites_VCF_allscaffolds.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(V5,V6)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V5","V6") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V5, V6,  color=Site, shape=Type), size=4) + ylab("PC4 (6.24%)") + xlab("PC3 (6.63%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 34 samples (4,522,615 SNPs)") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V5, V6,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V5, V6,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g
pdf("Plots/PCA_allsamples_allsites_pc3pc4_VCF_allscaffolds.pdf", height=8, width=15)
g
dev.off()

# PCA from snpAD + plink PROJECTING MODERN INTO HISTORICAL-----

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


# PCA from snpAD + S. geoffryi -------------

pca_plink <- read.table("Files/CTT_allsamples_filter_singletons.eigenvec") #not really working two outlier
val_plink <- read.table("Files/CTT_allsamples_filter_singletons.eigenval") 

pca_plink <- read.table("Files/CTT_allsamples_filter.eigenvec")
val_plink <- read.table("Files/CTT_allsamples_filter.eigenval")
val_plink[1,1]/sum(val_plink$V1) *100
val_plink[2,1]/sum(val_plink$V1) *100
val_plink[3,1]/sum(val_plink$V1) *100
val_plink[4,1]/sum(val_plink$V1) *100


pca_plink$Full_ID <- pca_plink$V1
pca_plink <- merge(pca_plink, metadata[ , 1:12], by="Full_ID")                            
centroids <- aggregate(cbind(V3,V4)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V3","V4") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V3, V4,  color=Site, shape=Type), size=4) + ylab("PC2 (7.52%)") + xlab("PC1 (18.11%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples + 1 S. geoffroyi (5,252,046 SNPs)") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V3, V4,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V3, V4,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_SG_maf0.05.pdf", height=8, width=15)
g
dev.off()

centroids <- aggregate(cbind(V5,V6)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V5","V6") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V5, V6,  color=Site, shape=Type), size=4) + ylab("PC4 (5.81%)") + xlab("PC3 (6.15%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 36 samples + 1 S. geoffroyi (5,252,046 SNPs)") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V5, V6,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V5, V6,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) 
#scale_color_manual(values=c("#f5ba18", "#451cb8") ) 
g

pdf("Plots/PCA_allsamples_SG_maf0.05_pc3pc4.pdf", height=8, width=15)
g
dev.off()



# Admixture 0.05 no prunning --------
# K = 2

order=c("Saguinus_geoffroyi", "Tulenapa","Mutatá", "Turbo","Tierralta",
        "Planeta Rica","Caracas","Cauca" ,"Coloso","San Juán", "Arjona","Ceibal",
        "Unknown")

file_Q=paste0("Files/Admixture/CTT_allsamples_filter.2.Q")
file_Q=paste0("Files/Admixture/CTT_allsamples_filter_singletons_LD.2.Q") # no maf, keep singletons 
file_Q=paste0("Files/Admixture/CTT_onlyCTT_filter.2.Q") # only CTT no maf

tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata[metadata$Site!="Saguinus_geoffroyi",]$Full_ID, 
                             metadata[metadata$Site!="Saguinus_geoffroyi",]$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k2 <- ggplot(df, aes(Sample,value, fill=key))
k2 <- k2 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k2

pdf("Plots/OnlyCTT_k2.pdf", width = 10, height = 3)
k2
dev.off()

# K = 3
file_Q=paste0("Files/Admixture/CTT_onlyCTT_filter.3.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata[metadata$Site!="Saguinus_geoffroyi",]$Full_ID, 
                             metadata[metadata$Site!="Saguinus_geoffroyi",]$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Sample", "Site")
df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k3 <- ggplot(df, aes(Sample,value, fill=key))
k3 <- k3 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k3
pdf("Plots/OnlyCTT_k3.pdf", width = 10, height = 3)
k3
dev.off()

# K = 4
file_Q=paste0("Files/Admixture/CTT_onlyCTT_filter.4.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata[metadata$Site!="Saguinus_geoffroyi",]$Full_ID, 
                             metadata[metadata$Site!="Saguinus_geoffroyi",]$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k4 <- ggplot(df, aes(Sample,value, fill=key))
k4 <- k4 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k4

pdf("Plots/OnlyCTT_k4.pdf", width = 10, height = 3)
k4
dev.off()

# K = 5
file_Q=paste0("Files/Admixture/CTT_onlyCTT_filter.5.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata[metadata$Site!="Saguinus_geoffroyi",]$Full_ID, 
                             metadata[metadata$Site!="Saguinus_geoffroyi",]$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Q5","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k5 <- ggplot(df, aes(Sample,value, fill=key))
k5 <- k5 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k5
pdf("Plots/OnlyCTT_k5.pdf", width = 10, height = 3)
k5
dev.off()

# K = 6
file_Q=paste0("Files/Admixture/CTT_onlyCTT_filter.6.Q")
tbl=read.table(file_Q)

mergedAdmixtureTable = cbind(tbl, metadata[metadata$Site!="Saguinus_geoffroyi",]$Full_ID, 
                             metadata[metadata$Site!="Saguinus_geoffroyi",]$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Q5","Q6","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k6 <- ggplot(df, aes(Sample,value, fill=key))
k6 <- k6 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
k6

pdf("Plots/OnlCTT_k6.pdf", width = 10, height = 3)
k6
dev.off()

# K = 7
file_Q=paste0("Files/Admixture/CTT_onlyCTT_filter.7.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata[metadata$Site!="Saguinus_geoffroyi",]$Full_ID, 
                             metadata[metadata$Site!="Saguinus_geoffroyi",]$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k7 <- ggplot(df, aes(Sample,value, fill=key))
k7 <- k7 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
k7

pdf("Plots/OnlyCTT_k7.pdf", width = 10, height = 3)
k7
dev.off()

# K = 8
file_Q=paste0("Files/Admixture/CTT_onlyCTT_filter.8.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata[metadata$Site!="Saguinus_geoffroyi",]$Full_ID, 
                             metadata[metadata$Site!="Saguinus_geoffroyi",]$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k8 <- ggplot(df, aes(Sample,value, fill=key))
k8 <- k8 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
k8

pdf("Plots/OnlyCTT_k8.pdf", width = 10, height = 3)
k8
dev.off()



# Admixture Singletons Subset no prunning --------
# K = 2
samples_subset <- read.table("Files/Admixture/Samples_subset")
order=c("Saguinus_geoffroyi", "Tulenapa","Mutatá", "Turbo","Tierralta",
        "Planeta Rica","Caracas","Cauca" ,"Coloso","San Juán", "Arjona","Ceibal",
        "Unknown")
file_Q=paste0("Files/Admixture/CTT_subset_filter.2.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, samples_subset$V1)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Full_ID")

mergedAdmixtureTable <- merge(mergedAdmixtureTable, metadata[ , c("Full_ID","Site")], by="Full_ID")                            
colnames(mergedAdmixtureTable) <- c("Sample","Q1","Q2","Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k2 <- ggplot(df, aes(Sample,value, fill=key))
k2 <- k2 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k2

pdf("Plots/All_k2.pdf", width = 25, height = 3)
k2
dev.off()

# K = 3
file_Q=paste0("Files/Admixture/CTT_allsamples_filter.3.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k3 <- ggplot(df, aes(Sample,value, fill=key))
k3 <- k3 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k3
pdf("Plots/All_k3.pdf", width = 25, height = 3)
k3
dev.off()

# K = 4
file_Q=paste0("Files/Admixture/CTT_allsamples_filter.4.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k4 <- ggplot(df, aes(Sample,value, fill=key))
k4 <- k4 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k4

pdf("Plots/All_k4.pdf", width = 25, height = 3)
k4
dev.off()

# K = 5
file_Q=paste0("Files/Admixture/CTT_allsamples_filter.5.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Q5","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k5 <- ggplot(df, aes(Sample,value, fill=key))
k5 <- k5 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k5
pdf("Plots/All_k5.pdf", width = 25, height = 3)
k5
dev.off()

# K = 6
file_Q=paste0("Files/Admixture/CTT_allsamples_filter.6.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Q5","Q6","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k6 <- ggplot(df, aes(Sample,value, fill=key))
k6 <- k6 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
k6

pdf("Plots/All_k6.pdf", width = 25, height = 3)
k6
dev.off()

# K = 7
file_Q=paste0("Files/Admixture/CTT_allsamples_filter.7.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k7 <- ggplot(df, aes(Sample,value, fill=key))
k7 <- k7 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
k7

pdf("Plots/All_k7.pdf", width = 25, height = 3)
k7
dev.off()






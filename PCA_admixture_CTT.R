library("ggplot2")
library("tidyr")
library("ggConvexHull")
library("scales")
library("ggrepel")
library("RColorBrewer")
library("viridis")
library("readxl")
library("ggpubr")

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
metadata <- read_excel("Paper/SupplementaryTables_v1.1.xlsx", sheet = 4)
colors <- c("#645244ff", "#832232ff", "#645244ff","#70a37f","#b3dec1ff","#fbefa6ff","#fab2eaff", "#6d466bff", "#508991ff", "#0a4a33ff", "#eb9486ff","#645244ff" )

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
  theme(legend.position="right", legend.title = element_blank()) +
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)
  
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
  theme(legend.position="right", legend.title = element_blank()) +
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)
g

pdf("Plots/PCA_onlyCTT_allsites_maf0.05_pc3pc4.pdf", height=8, width=15)
g
dev.off()


# PCA from snpAD + plink all autosomes maf 0.05 ----------------
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

pc1 <- ggplot()
pc1 <- pc1 +   geom_point(data=pca_plink, aes(V3, V4, color=Site, shape=Type), size=4) + ylab("PC2 (7.65%)") + xlab("PC1 (12.32%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 34 samples (4,522,615 SNPs)") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V3, V4,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V3, V4,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  +
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  theme(legend.position="right", legend.title = element_blank()) 
pc1
pdf("Plots/PCA_allsamples_allsites_VCF_allscaffolds.pdf", height=8, width=15)
pc1
dev.off()

centroids <- aggregate(cbind(V5,V6)~Site,pca_plink,mean)
colnames(centroids) <- c("Site","V5","V6") 

g <- ggplot()
g <- g +   geom_point(data=pca_plink, aes(V5, V6,  color=Site, shape=Type), size=4) + ylab("PC4 (6.24%)") + xlab("PC3 (6.63%)") + 
  ggtitle("PCA - Cotton-Top Tamarin 34 samples (4,522,615 SNPs)") + theme_test() + 
  geom_convexhull(data=pca_plink, aes(V5, V6,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V5, V6,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1)  + 
  theme(legend.position="right", legend.title = element_blank()) +
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)
g
pdf("Plots/PCA_allsamples_allsites_pc3pc4_VCF_allscaffolds.pdf", height=8, width=15)
g
dev.off()

# ADMIXTURE with maf 0.05 --------
# Cross validation 
allfiles <- list.files("Files/Admixture/CV/", pattern=".out", full.names=TRUE) 
CV_list <- lapply(allfiles,
                 FUN = function(files) {
                   read.table(files, header = FALSE, sep = ",")
                 })
rep <- gsub(".out","",basename(allfiles))

for (i in 1:length(rep)){
  CV_list[[i]]$Rep <- rep[i]
  CV_list[[i]]$K <- c(1,10,2:9)
}

CV_df <- do.call(rbind, CV_list)

# all
g <- ggplot(CV_df, aes(K, V1, color=Rep))
g <- g  +geom_point() +theme_classic() + ylab("CV Error")
g

pdf("Plots/CV_error.pdf", height = 4, width = 5)
g
dev.off()
# chosen
CV_df[CV_df$V1==min(CV_df[CV_df$K==2,]$V1),]
CV_df[CV_df$V1==min(CV_df[CV_df$K==3,]$V1),]
CV_df[CV_df$V1==min(CV_df[CV_df$K==4,]$V1),]
CV_df[CV_df$V1==min(CV_df[CV_df$K==5,]$V1),]
CV_df[CV_df$V1==min(CV_df[CV_df$K==6,]$V1),]
CV_df[CV_df$V1==min(CV_df[CV_df$K==7,]$V1),]
CV_df[CV_df$V1==min(CV_df[CV_df$K==8,]$V1),]
CV_df[CV_df$V1==min(CV_df[CV_df$K==9,]$V1),]
CV_df[CV_df$V1==min(CV_df[CV_df$K==10,]$V1),]


order=c("Ceibal","San Juan","Colosó", "Caracas",
        "Planeta Rica","Tierralta","Turbo","Mutatá","Tulenapa",
        "Arjona_Uncertain","Cauca_Uncertain","Unknown")

# K = 2
file_Q=paste0("Files/Admixture/Rep6/CTT_onlyCTT_filter.2.Q") 
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q2","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k2 <- ggplot(df, aes(Sample,value, fill=key))
k2 <- k2 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_blank(),
        strip.text.x = element_text(angle=90, hjust=0),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k2
pdf("Plots/OnlyCTT_k2.pdf", width = 10, height = 3)
k2
dev.off()

# K = 3
file_Q=paste0("Files/Admixture/Rep8/CTT_onlyCTT_filter.3.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q2","Q1","Q3","Sample", "Site")
df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k3 <- ggplot(df, aes(Sample,value, fill=key))
k3 <- k3 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k3
pdf("Plots/OnlyCTT_k3.pdf", width = 10, height = 3)
k3
dev.off()

# K = 4
file_Q=paste0("Files/Admixture/Rep8/CTT_onlyCTT_filter.4.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q2","Q1","Q3","Q4","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k4 <- ggplot(df, aes(Sample,value, fill=key))
k4 <- k4 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
pdf("Plots/OnlyCTT_k4.pdf", width = 10, height = 3)
k4
dev.off()

# K = 5
file_Q=paste0("Files/Admixture/Rep4/CTT_onlyCTT_filter.5.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q1","Q5","Q3","Q2","Q4","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k5 <- ggplot(df, aes(Sample,value, fill=key))
k5 <- k5 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
pdf("Plots/OnlyCTT_k5.pdf", width = 10, height = 3)
k5
dev.off()

# K = 6
file_Q=paste0("Files/Admixture/Rep5/CTT_onlyCTT_filter.6.Q")
tbl=read.table(file_Q)

mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q2","Q4","Q3","Q1","Q6","Q5","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)

k6 <- ggplot(df, aes(Sample,value, fill=key))
k6 <- k6 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x =element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
pdf("Plots/OnlyCTT_k6.pdf", width = 10, height = 3)
k6
dev.off()

# K = 7
file_Q=paste0("Files/Admixture/Rep6/CTT_onlyCTT_filter.7.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl,  metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q4","Q1","Q2","Q5","Q3","Q7","Q6","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k7 <- ggplot(df, aes(Sample,value, fill=key))
k7 <- k7 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
pdf("Plots/OnlyCTT_k7.pdf", width = 10, height = 3)
k7
dev.off()

# K = 8
file_Q=paste0("Files/Admixture/Rep6/CTT_onlyCTT_filter.8.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl,  metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q2","Q4","Q3","Q6","Q8","Q1","Q5","Q7","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k8 <- ggplot(df, aes(Sample,value, fill=key))
k8 <- k8 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x =element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 
k8
pdf("Plots/OnlyCTT_k8.pdf", width = 10, height = 3)
k8
dev.off()


# K = 9
file_Q=paste0("Files/Admixture/Rep3/CTT_onlyCTT_filter.9.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl,  metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q5","Q8","Q3","Q1","Q2","Q7","Q6","Q4","Q9","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k9 <- ggplot(df, aes(Sample,value, fill=key))
k9 <- k9 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k9
pdf("Plots/OnlyCTT_k9.pdf", width = 10, height = 3)
k9
dev.off()

# K = 10
file_Q=paste0("Files/Admixture/Rep9/CTT_onlyCTT_filter.10.Q")
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl,  metadata$Full_ID, metadata$Site)
colnames(mergedAdmixtureTable) <- c("Q9","Q8","Q1","Q2","Q3","Q7","Q6","Q10","Q4","Q5","Sample", "Site")

df <- gather(mergedAdmixtureTable, key, value, -Sample, -Site)
df$Site  <- factor(df$Site,levels=order,ordered=TRUE)


k10 <- ggplot(df, aes(Sample,value, fill=key))
k10 <- k10 + geom_col(position = "fill",width = 1) + facet_grid(~Site, space="free", scales="free_x") + theme_minimal() + xlab("") + 
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = NA, color=NA))+
  scale_y_continuous(expand = c(0, 0)) + scale_fill_brewer(palette = "Paired") + ylab("")+
  scale_x_discrete(expand = c(0,0)) 

k10

pdf("Plots/OnlyCTT_k10df", width = 10, height = 3)
k10
dev.off()



pdf("Plots/Admixture_allK.pdf", height = 14, width = 8)
ggarrange(k2,k3,k4,k5,k6,k7,k8,k9,k10, ncol=1, heights = c(0.8,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.8))
dev.off()



# FIGURE 1 ----
library(png)
map <- readPNG("HCTT_mapEdited.png")
im_A <- ggplot() + 
  background_image(map) +coord_fixed(2/1.5)

library(cowplot)
my_legend <- get_legend(pc1)
legend <- as_ggplot(my_legend)

pc1 <- ggplot()
pc1 <- pc1 +   geom_point(data=pca_plink, aes(V3, V4, color=Site, shape=Type), size=4) + ylab("PC2 (7.65%)") + xlab("PC1 (12.32%)") + 
  theme_test() +coord_flip() +  
  geom_convexhull(data=pca_plink, aes(V3, V4,  fill=Site),alpha=.3, colour = NA)+
  geom_text_repel(data=centroids, aes(V3, V4,  label=Site), force = 4,max.overlaps = Inf,box.padding = 1.2)  +
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  theme(legend.position="none", legend.title = element_blank()) 
pc1


file_Q=paste0("Files/Admixture/Rep6/CTT_onlyCTT_filter.2.Q") 
tbl=read.table(file_Q)
mergedAdmixtureTable = cbind(tbl, metadata$Full_ID, metadata$Site)
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

pdf("Plots/Figure1_new.pdf", height = 10, width = 8)
ggarrange(ggarrange(im_A,pc1, ncol=2, labels = c("A","B")),ncol = 1, legend, k2,labels = c("","","C"), heights = c(4,0.5,2))
dev.off()


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

# F INBREEDING PLINK -----------

g <- ggplot(metadata, aes(Site, F_plink, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() +
  theme(legend.position = "none")
g

pdf("Plots/Inbreeding_F.pdf", height=4, width=11)
g
dev.off()

g <- ggplot(metadata[metadata$Coverage>5,], aes(Site, F_plink, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() +ggtitle("Only Samples > 5x")+
  theme(legend.position = "none")
g

pdf("Plots/Inbreeding_F_highCov.pdf", height=4, width=11)
g
dev.off()

g <- ggplot(metadata, aes(Coverage, F_plink, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() 
g

pdf("Plots/Inbreeding_FvsCoverage.pdf", height=5, width=6)
g
dev.off()

g <- ggplot(metadata, aes(Missingess, F_plink, color=Site, shape=Type))
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

Ç
# F INBREEDING NGSRELATE to do-----

# Heterozygosity with SNPs ----

het <- list.files("Files/Het/scaffolds/", pattern = ".het", full.names = TRUE)
ldf_het <- lapply(het, read.table)

samplesid <- unlist(lapply(1:length(basename(het)), function(i) strsplit(basename(het)[i],"_t")[[1]][1]))

for (i in 1:length(samplesid)){
  Full_ID <- samplesid[i]
  ldf_het[[i]]<- cbind.data.frame(ldf_het[[i]], Full_ID=as.character(Full_ID))
}


df_het <- do.call(rbind.data.frame,ldf_het)

df_het1 <- merge(df_het, metadata[ , 1:12], by="Full_ID")                            

g <- ggplot(df_het1, aes(Full_ID, V4, fill=Type))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Site, space="free", scales="free_x") + 
  theme_classic() +
  theme(legend.position = "none")
g


g <- ggplot(df_het1, aes(V1, V4))
g <- g+ geom_boxplot() + geom_jitter()+   
  theme_classic() +
  theme(legend.position = "none")
g

#Calculate averag per scaffold and per individual
average_het <- list()
for (i in 1:length(samplesid)){
  Full_ID <- samplesid[i]
  average_het[[i]]<- cbind.data.frame(Full_ID=as.character(Full_ID), 
                                      MeanHet=mean(ldf_het[[i]]$V4),
                                      MedianHet=median(ldf_het[[i]]$V4))
}
df_avhet <- do.call(rbind.data.frame,average_het)

#save the values in metadata

g <- ggplot(metadata, aes(Site, `Median Het`, fill=Site))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() + ylab("Median Genome-wide Heterozygosity (bp-1)")+
  theme(legend.position = "none")
g

pdf("Plots/Heterozygosity_new.pdf", height=4, width=11)
g
dev.off()

g <- ggplot(metadata[metadata$Coverage>5,], aes(Site, `Median Het`, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() + ggtitle("Only Samples > 5x")+
  theme(legend.position = "none")
g

pdf("Plots/Heterozygosity_highCov.pdf", height=4, width=11)
g
dev.off()

g <- ggplot(metadata, aes(Missingess, `Median Het`, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() + xlab("Missingness")
g

pdf("Plots/HeterozygosityvsMissingness.pdf", height=5, width=6)
g
dev.off()

g <- ggplot(metadata, aes(Coverage, `Median Het`, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() 
g

pdf("Plots/HeterozygosityvsCoverage.pdf", height=5, width=6)
g
dev.off()

# Het from ANGSD -----
fins <- list.files("Files/Het/angsd/",pattern=".ml",full.names=TRUE)
samples <- sub('\\.ml$', '', basename(fins))

fold <- 0 # Folded caculation, since we use the same ref and anc.
df_list <- lapply(fins,
                  FUN = function(files) {
                    scan(files)
                  })
df_list2 <- list()
for (i in 1:length(samples)){
  
  df_list2[[i]] <- cbind.data.frame(Sample=samples[i], het=df_list[[i]][2]/sum(df_list[[i]]))
}

het_df <- do.call(rbind,df_list2)
head(het_df)
#save the values in metadata
g <- ggplot(metadata, aes(Site, `Angsd Het`, fill=Site))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() + ylab("Heterozygosity (bp-1)")+
  theme(legend.position = "none")
g

pdf("Plots/Heterozygosity_angsd.pdf", height=4, width=11)
g
dev.off()


g <- ggplot(metadata, aes(Coverage, `Angsd Het`, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() 
g
pdf("Plots/Heterozygosity_angsd_coverage.pdf", height=4, width=11)
g
dev.off()

ctt <- metadata[metadata$Site!="Saguinus_geoffroyi",]
test <- t.test(ctt[ctt$Type=="Historical",]$`Angsd Het`, 
               ctt[ctt$Type=="Modern",]$`Angsd Het`)

test$p.value
library(ggpubr)
compare_means(formula, data, method = "wilcox.test", paired = FALSE,
              group.by = NULL, ref.group = NULL, ...)

g <- ggplot(ctt, aes(Type, `Angsd Het`))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+
  stat_compare_means(label.x = 1.5, label.y = 0.0014)+
  theme_classic() + xlab("") + ylab("Heterozygosity (bp-1)")
g

pdf("Plots/Heterozygosity_testHistModern.pdf", height=4, width=7)
g
dev.off()

g <- ggplot(ctt[ctt$Coverage > 5,], aes(Type, `Angsd Het`))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+
  stat_compare_means(label.x = 1.5, label.y = 0.0014)+
  theme_classic() + xlab("") + ylab("Heterozygosity (bp-1)")
g

pdf("Plots/Heterozygosity_testHistModern_highcov.pdf", height=4, width=7)
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
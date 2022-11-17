library("ggplot2")
library("tidyr")
library("ggConvexHull")
library("scales")
library("ggrepel")
library("RColorBrewer")
library("viridis")
library("readxl")

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
metadata <- read_excel("Shotgun_Metadata.xlsx", sheet = 1)
metadata <- metadata[-which(metadata$Full_ID=="37_AMNH_33174"|metadata$Full_ID=="CEI_060"),]

colors = c("#645244","#832232","grey","#70a37f","#b3dec1","#fbefa6","#fab2ea","#6d466b",
           "#508991","#0a4a33", "#eb9486","grey")
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


# F INBREEDING NGSRELATE to do-----

g <- ggplot(metadata[metadata$Full_ID!="CEI_060",], aes(Site, `F NgsRelate`, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() + ylab("Inbreeding (F)")+ scale_fill_manual(values=colors)+
  theme(legend.position = "none")
g

pdf("Plots/Inbreeding_F_ngsrelate.pdf", height=4, width=11)
g
dev.off()


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
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ 
  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() + ylab("Heterozygosity (bp-1)")+
  theme(legend.position = "none") + scale_fill_manual(values=colors)
g

pdf("Plots/Heterozygosity_angsd.pdf", height=4, width=11)
g
dev.off()


g <- ggplot(metadata, aes(Coverage, `Angsd Het`, color=Site, shape=Type))
g <- g+ geom_point() +
  theme_classic() + scale_color_manual(values=colors)
g
pdf("Plots/Heterozygosity_angsd_coverage.pdf", height=4, width=11)
g
dev.off()

ctt <- metadata
test <- t.test(ctt[ctt$Type=="Historical",]$`Angsd Het`, 
               ctt[ctt$Type=="Modern",]$`Angsd Het`)

test$p.value
library(ggpubr)
compare_means(formula, data, method = "wilcox.test", paired = FALSE,
              group.by = NULL, ref.group = NULL, ...)

g <- ggplot(ctt, aes(Type, `Angsd Het`, fill=Type))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ 
  stat_compare_means(label.x = 1.5, label.y = 0.0014)+
  theme_classic() + xlab("") + ylab("Heterozygosity (bp-1)") + 
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))
g

pdf("Plots/Heterozygosity_testHistModern.pdf", height=4, width=7)
g
dev.off()

mean(ctt[ctt$Group!="Unknown"&ctt$Type=="Historical",]$`Angsd Het`)/mean(ctt[ctt$Group!="Unknown"&ctt$Type=="Modern",]$`Angsd Het`)
mean(ctt[ctt$Group!="Unknown"&ctt$Type=="Historical"&ctt$Group=="Greater Northeast",]$`Angsd Het`)/mean(ctt[ctt$Group!="Unknown"&ctt$Type=="Modern"&ctt$Group=="Greater Northeast",]$`Angsd Het`)
mean(ctt[ctt$Group!="Unknown"&ctt$Type=="Historical"&ctt$Group=="Southwest",]$`Angsd Het`)/mean(ctt[ctt$Group!="Unknown"&ctt$Type=="Modern"&ctt$Group=="Southwest",]$`Angsd Het`)


g11 <- ggplot(ctt[ctt$Group!="Unknown",], aes(Type, `Angsd Het`, fill=Type))
g11 <- g11+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ facet_grid(.~Group)+
  stat_compare_means(label.x = 1.5, label.y = 0.0014)+
  theme_classic() + xlab("") + ylab("Heterozygosity (bp-1)")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))

g11

pdf("Plots/Heterozygosity_testHistModern_group.pdf", height=4, width=7)
g11
dev.off()


g <- ggplot(ctt[ctt$Coverage > 7&ctt$Group!="Unknown",], aes(Type, `Angsd Het`, fill=Type))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ facet_grid(.~Group)+
  stat_compare_means(label.x = 1.5, label.y = 0.0014)+
  theme_classic() + xlab("") + ylab("Heterozygosity (bp-1)")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))

g

pdf("Plots/Heterozygosity_testHistModern_highcov.pdf", height=4, width=7)
g
dev.off()

# RoH for the paper

library(ggplot2)
library(scales)
library("readxl")
library("ggpubr")
library(tidyr)
setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
metadata  <- read_excel("Shotgun_Metadata.xlsx")

metadata_8x <- metadata[metadata$Site!="Saguinus_geoffroyi"&metadata$Coverage>8&metadata$Full_ID!="CEI_060",]


# ROHan summary length RoHs -----
roh <- list.files("Files/ROHan/", pattern = "_2.5e4_allscaffolds_1Mb.mid.hmmrohl", full.names = TRUE)# This one
ldf_het <- lapply(roh, read.table)
V1 <- gsub("_2.5e4_allscaffolds_1Mb.mid.hmmrohl","", basename(roh)) # This one
finalV1 <- V1

for (i in 1:length(finalV1)){
  ind <- finalV1[i]
  ldf_het[[i]]<- cbind.data.frame(ldf_het[[i]], ind=as.character(ind))
}
df_roh <- do.call(rbind.data.frame,ldf_het)
colnames(df_roh) <- c("V1", "Scaffold","BEGIN","END" ,"ROH_LENGTH","VALIDATED_SITES","Sample")

assembly <- read.table("Files/assembly.fai")
colnames(assembly) <- c("Scaffold","Length","Pos","V4","V5")

roh_final <- merge(df_roh, assembly[,c("Scaffold","Length"),], by="Scaffold")
roh_final_rohan <- roh_final

list_roh <- list()
# per sample  -----
for (i in 1:length(finalV1)){  
  sample <- as.character(finalV1[i])
  df <- roh_final[roh_final$Sample==sample,]
  df_het_sample <- df_het[df_het$Sample == sample,]
  p1 <- ggplot(df)
  p1 <- p1 + facet_wrap(.~Scaffold, ncol=1,strip.position="right")+ ggtitle(sample) + 
    ylab("")+ xlab("Position (bp)")+
    scale_x_continuous(name = "", waiver(), labels=comma, expand = c(0.001,0), 
                       limits = c(0,224379228))+
    scale_y_continuous(name = "", waiver(), labels=comma,
                       limits = c(0,0.01))+
    geom_rect(data=df,
              mapping=aes(xmin = BEGIN, ymin = 0, xmax = END, ymax = 0.01), 
              fill="#913C5C", alpha=0.8) + 
    geom_rect(data=assembly,mapping = aes(xmin=0, xmax=Length, ymin=0, ymax=0.01), alpha = .3,
              color = "#2e2e2e", fill = NA) +
    theme_classic()+
    theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),panel.grid.minor = element_blank(),
          legend.position = "none", axis.ticks.y = element_blank(), panel.border=element_blank(),
          panel.background = element_blank(),axis.line = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank())
  print(p1)
  list_roh[[i]] <- p1
  ggsave(paste0("Plots/RoHs/RoHans_", sample, "2.5e-4_mid_1Mb.pdf"), width = 6, height = 6)
}

ggarrange(plotlist=list_roh, labels = c("A", "B", "C","D","E"))
ggsave("Plots/RoH_final_rohan.pdf", height = 15, width = 13)




# summary ROH
samples_new <- read.table("Files/Samples", header = FALSE)

roh_perc_bcftools <- lapply(1:length(samples_new$V1), function(i) {
  df2 <- roh_final[roh_final$Sample == samples_new$V1[i],]
  data.frame(sample=samples_new$V1[i],
             Total=sum(df2$ROH_LENGTH)/2779462738,  
             Average=mean(df2$ROH_LENGTH),
             Intermediate=sum(df2[df2$ROH_LENGTH>=1000000 & df2$ROH_LENGTH<2500000,]$ROH_LENGTH)/2779462738, 
             Intermediate2=sum(df2[df2$ROH_LENGTH>=2500000 & df2$ROH_LENGTH<5000000,]$ROH_LENGTH)/2779462738, 
             Long=sum(df2[df2$ROH_LENGTH>=5000000 & df2$ROH_LENGTH<8000000,]$ROH_LENGTH)/2779462738, 
             ExtraLong=sum(df2[df2$ROH_LENGTH>=8000000,]$ROH_LENGTH)/2779462738, 
             count=length(df2$ROH_LENGTH),
             sumMb=sum(df2$ROH_LENGTH))
})
roh_perc_df_bcftools <- do.call(rbind, roh_perc_bcftools)

roh_perc_df_gather <- gather(roh_perc_df_bcftools, key, value, -sample, -Average, -count,-sumMb)
roh_perc_df_gather$key <- gsub("Intermediate2","1-2.5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate","2.5-5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("ExtraLong",">8Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long","5-8Mb", roh_perc_df_gather$key )
order <- c("Total", "1-2.5Mb","2.5-5Mb", "5-8Mb",">8Mb")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$Full_ID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="Full_ID")



roh_perc_df_gather3 <- roh_perc_df_gather2[which(roh_perc_df_gather2$Coverage>8&
                                                   roh_perc_df_gather2$key!="Total"&
                                                   roh_perc_df_gather2$key!="Average"&
                                                   roh_perc_df_gather2$sample!="CEI_060"),]

pa <- ggplot(roh_perc_df_gather3[roh_perc_df_gather3$value>0,], 
             aes(sample,value*100,fill=key))
pa <- pa + geom_col()+ theme_classic() +  xlab("")+ facet_grid(.~Group, space = "free", scales="free")+
  ylab("% of Genome in RoH") + 
  theme(legend.position = "right", axis.text.y = element_text(hjust=0.5)) + 
  scale_fill_manual(values=rev(c("#fcd021","#234f96","#5e97f2","#b6d1fc")), name="RoH Size") 
pa

pdf("Plots/Summary_rohs.pdf", width = 6, height = 5)
p
dev.off()


# Final plot for figure 3
ggarrange(g11,pa,nrow = 2, align = "v",labels=c("A","B"))
ggsave("Plots/Figure3.pdf", height = 8, width = 8)

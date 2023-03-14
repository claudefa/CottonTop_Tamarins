library("ggplot2")
library("tidyr")
library("ggConvexHull")
library("scales")
library("ggrepel")
library("RColorBrewer")
library("readxl")
library("ggpubr")

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
metadata <- read_excel("Paper/SupplementaryTables_v1.1.xlsx", sheet = 4)

metadata <- metadata[-which(metadata$Full_ID=="CEI_060"),]

colors = c("#645244","#832232","grey","#70a37f","#b3dec1","#fbefa6","#fab2ea","#6d466b",
           "#508991","#0a4a33", "#eb9486","grey")

# F INBREEDING NGSRELATE -----

g <- ggplot(metadata, aes(Site, `Fi (NGSrelate)`, fill=Site))
g <- g+ geom_boxplot() + geom_jitter()+  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() + ylab("Inbreeding (Fi)")+ scale_fill_manual(values=colors)+
  theme(legend.position = "none")
g

pdf("Plots/Inbreeding_F_ngsrelate.pdf", height=4, width=11)
g
dev.off()

# HETEROZYGOSITY from ANGSD -----
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
g <- ggplot(metadata, aes(Site, `Heterozygosity (ANGSD)`, fill=Site))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ 
  facet_grid(.~Type, space="free", scales="free_x") + 
  theme_classic() + ylab("Heterozygosity (bp-1)")+
  theme(legend.position = "none") + scale_fill_manual(values=colors)
g

pdf("Plots/Heterozygosity_angsd.pdf", height=4, width=11)
g
dev.off()


cor.test(metadata$Coverage,metadata$`Heterozygosity (ANGSD)`)
summary(lm(metadata$Coverage~metadata$`Heterozygosity (ANGSD)`))

g <- ggplot(metadata, aes(Coverage, `Heterozygosity (ANGSD)`))
g <- g+ geom_point() + geom_smooth(method = "lm")+ ylab("Heterozygosity (bp-1)")+
  theme_classic() + scale_color_manual(values=colors)
g

pdf("Plots/Heterozygosity_angsd_coverage.pdf", height=4, width=6)
g
dev.off()

ctt <- metadata
test <- t.test(ctt[ctt$Type=="Historical",]$`Heterozygosity (ANGSD)`, 
               ctt[ctt$Type=="Modern",]$`Heterozygosity (ANGSD)`)
test$p.value

g <- ggplot(ctt, aes(Type, `Heterozygosity (ANGSD)`, fill=Type))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ 
  stat_compare_means(label.x = 1.5, label.y = 0.0014)+
  theme_classic() + xlab("") + ylab("Heterozygosity (bp-1)") + 
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))
g

pdf("Plots/Heterozygosity_testHistModern.pdf", height=4, width=7)
g
dev.off()

groups <- c("Southwest", "Greater Northeast")
mean(ctt[which(ctt$`Grouping (2 Regions)`%in%groups&ctt$Type=="Historical"),]$`Heterozygosity (ANGSD)`)/
  mean(ctt[which(ctt$`Grouping (2 Regions)`%in%groups&ctt$Type=="Modern"),]$`Heterozygosity (ANGSD)`)
mean(ctt[ctt$Type=="Historical",]$`Heterozygosity (ANGSD)`)/mean(ctt[ctt$Type=="Modern",]$`Heterozygosity (ANGSD)`)

mean(ctt[ctt$`Grouping (2 Regions)`!="Unknown"&ctt$Type=="Historical"&ctt$`Grouping (2 Regions)`=="Greater Northeast",]$`Heterozygosity (ANGSD)`)/mean(ctt[ctt$`Grouping (2 Regions)`!="Unknown"&ctt$Type=="Modern"&ctt$`Grouping (2 Regions)`=="Greater Northeast",]$`Heterozygosity (ANGSD)`)
mean(ctt[ctt$`Grouping (2 Regions)`!="Unknown"&ctt$Type=="Historical"&ctt$`Grouping (2 Regions)`=="Southwest",]$`Heterozygosity (ANGSD)`)/mean(ctt[ctt$`Grouping (2 Regions)`!="Unknown"&ctt$Type=="Modern"&ctt$`Grouping (2 Regions)`=="Southwest",]$`Heterozygosity (ANGSD)`)


g11 <- ggplot(ctt[which(ctt$`Grouping (2 Regions)`%in%groups),], aes(Type, `Heterozygosity (ANGSD)`, fill=Type))
g11 <- g11+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ facet_grid(.~`Grouping (2 Regions)`)+
  stat_compare_means(label.x = 1.5, label.y = 0.0014)+
  theme_classic() + xlab("") + ylab("Heterozygosity (bp-1)")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))

g11

pdf("Plots/Heterozygosity_testHistModern_group.pdf", height=4, width=7)
g11
dev.off()


g <- ggplot(ctt[which(ctt$`Grouping (2 Regions)`%in%groups & ctt$Coverage > 7),], aes(Type, `Heterozygosity (ANGSD)`, fill=Type))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ facet_grid(.~`Grouping (2 Regions)`)+
  stat_compare_means(label.x = 1.5, label.y = 0.0014)+
  theme_classic() + xlab("") + ylab("Heterozygosity (bp-1)")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))

g

pdf("Plots/Heterozygosity_testHistModern_highcov.pdf", height=4, width=7)
g
dev.off()

# RUNS OF HOMOGYZOSITY form ROHan -------
#only samples with coverage >7x with Rohan
samples <- c("24_FMNH_69286","36_AMNH_32702","39_AMNH_130390","41_FMNH_69938","CEI_051","SJN_001","SJN_011","TLP2_03","TLP9_01")
metadata_7x <- metadata[which(metadata$Full_ID%in%samples),]

# ROHan summary length RoHs -----
roh <- list.files("Files/Het/rohan/", pattern = "_aDNA_2e5_allscaffolds_1Mb.mid.hmmrohl", full.names = TRUE)
ldf_het <- lapply(roh, read.table)
V1 <- gsub("_aDNA_2e5_allscaffolds_1Mb.mid.hmmrohl","", basename(roh)) # This one
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

# per sample  -----
list_roh<- list()
for (i in 1:length(finalV1)){  
  sample <- as.character(finalV1[i])
  df <- roh_final_rohan[roh_final_rohan$Sample==sample,]
  p1 <- ggplot(df)
  p1 <- p1 + facet_wrap(.~Scaffold, ncol=1,strip.position="right")+ ggtitle(sample) + 
    ylab("")+ xlab("Position (bp)")+
    scale_x_continuous(name = "", waiver(), labels=comma, expand = c(0.001,0), 
                       limits = c(0,224379228))+
    scale_y_continuous(name = "", waiver(), labels=comma,
                       limits = c(0,0.01))+
    geom_rect(data=df,
              mapping=aes(xmin = BEGIN, ymin = 0, xmax = END, ymax = 0.01), 
              fill="#fc7f03", alpha=0.8) + 
    geom_rect(data=assembly,mapping = aes(xmin=0, xmax=Length, ymin=0, ymax=0.01), alpha = .3,
              color = "#2e2e2e", fill = NA) +
    theme_classic()+
    theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),panel.grid.minor = element_blank(),
          legend.position = "none", axis.ticks.y = element_blank(), panel.border=element_blank(),
          panel.background = element_blank(),axis.line = element_blank(), 
          strip.background = element_blank(), strip.text.y = element_text(angle=0, hjust=1))
  print(p1)
  list_roh[[i]] <- p1
}


list_roh_plot <- c(list_roh[4],list_roh[6],list_roh[7:13])

ggarrange(plotlist=list_roh_plot, labels = c("A", "B", "C","D","E", "F","G","H","I"))
ggsave("Plots/RoH_final_rohan.pdf", height = 15, width = 13) 

# Age of fragments ------
roh_final_rohan$Full_ID <- roh_final_rohan$Sample
roh_final_rohan_metadata <- merge(roh_final_rohan, metadata, by="Full_ID")
roh_final_rohan_metadata$Age <- 100/(2*(roh_final_rohan_metadata$ROH_LENGTH/1000000))*6
roh_final_rohan_metadata$Year <- 2020-roh_final_rohan_metadata$Age

roh_final_rohan_metadata <- roh_final_rohan_metadata[roh_final_rohan_metadata$Coverage > 7,] 

roh_final_rohan_metadata$LengthMb <- roh_final_rohan_metadata$ROH_LENGTH/1000000

mean(roh_final_rohan_metadata[roh_final_rohan_metadata$Full_ID=="CEI_051",]$Age)
min(roh_final_rohan_metadata[roh_final_rohan_metadata$Full_ID=="CEI_051",]$Age)

median(roh_final_rohan_metadata[roh_final_rohan_metadata$Full_ID=="CEI_051",]$Age)
mean(roh_final_rohan_metadata[roh_final_rohan_metadata$Full_ID=="CEI_051",]$ROH_LENGTH)/1000000


# summary ROH ------
samples_new <- read.table("Files/Samples", header = FALSE)

roh_perc_bcftools <- lapply(1:length(samples_new$V1), function(i) {
  df2 <- roh_final_rohan[roh_final_rohan$Sample == samples_new$V1[i],]
  data.frame(sample=samples_new$V1[i],
             Total=sum(df2$ROH_LENGTH)/2779462738,  
             Average=mean(df2$ROH_LENGTH),
             Median=median(df2$ROH_LENGTH),
             Intermediate=sum(df2[df2$ROH_LENGTH>=1000000 & df2$ROH_LENGTH<2500000,]$ROH_LENGTH)/2779462738, 
             Intermediate2=sum(df2[df2$ROH_LENGTH>=2500000 & df2$ROH_LENGTH<5000000,]$ROH_LENGTH)/2779462738, 
             Long=sum(df2[df2$ROH_LENGTH>=5000000 & df2$ROH_LENGTH<8000000,]$ROH_LENGTH)/2779462738, 
             ExtraLong=sum(df2[df2$ROH_LENGTH>=8000000,]$ROH_LENGTH)/2779462738, 
             count=length(df2$ROH_LENGTH),
             sumMb=sum(df2$ROH_LENGTH),
             maxMb=max(df2$ROH_LENGTH))
})
roh_perc_df_bcftools <- do.call(rbind, roh_perc_bcftools)

roh_perc_df_gather <- gather(roh_perc_df_bcftools, key, value, -sample, -Median,-Average, -count,-sumMb,-maxMb)
roh_perc_df_gather$key <- gsub("Intermediate2","1-2.5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate","2.5-5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("ExtraLong",">8Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long","5-8Mb", roh_perc_df_gather$key )
order <- c("Total", "1-2.5Mb","2.5-5Mb", "5-8Mb",">8Mb")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$Full_ID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="Full_ID")



roh_perc_df_gather3 <- roh_perc_df_gather2[which(roh_perc_df_gather2$Coverage>7&
                                                   roh_perc_df_gather2$key!="Total"&
                                                   roh_perc_df_gather2$key!="Average"&
                                                   roh_perc_df_gather2$sample!="CEI_060"),]

pa <- ggplot(roh_perc_df_gather3, 
             aes(sample,value*100,fill=key))
pa <- pa + geom_col()+ theme_classic() +  xlab("")+ facet_grid(.~Type,  scales="free")+
  ylab("% of Genome in RoH") + 
  theme(legend.position = "right", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=rev(c("#fcd021","#234f96","#5e97f2","#b6d1fc")), name="RoH Size") 
pa

pdf("Plots/Summary_rohs.pdf", width = 6, height = 5)
pa
dev.off()

# Final plot for figure 3
ggarrange(g11,pa,nrow = 2, align = "v",labels=c("A","B"))
ggsave("Plots/Figure3.pdf", height = 8, width = 8)

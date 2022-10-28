# Runs of homozygosity
library(ggplot2)
library(scales)
library("readxl")
library("ggpubr")
setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
metadata  <- read_excel("Shotgun_Metadata.xlsx")

#Coverage ----

g <- ggplot(metadata, aes(Type,Coverage))
g <- g+ geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_classic()
g

# Het per window ---------
het <- list.files("Files/Het/window/", pattern = ".het", full.names = TRUE)
ldf_het <- lapply(het, read.table)
V1 <- read.table("Files/Samples", header = FALSE)
for (i in 1:length(samples$V1)){
  ind <- V1$V1[i]
  ldf_het[[i]]<- cbind.data.frame(ldf_het[[i]], ind=as.character(ind))
}

df_het <- do.call(rbind.data.frame,ldf_het)
colnames(df_het) <- c("Scaffold", "Pos", "NumHet", "Callable", "Heterozygosity","Sample")

scaffolds_fin <- c( "CM038391.1", "CM038392.1" ,"CM038393.1", "CM038394.1", "CM038395.1")
df_het_f <- df_het[(df_het$Scaffold%in%scaffolds_fin),]

# Density heterozygosity
#Remove all the windows with two standard deviation below the median
nocall <- median(df_het_f$Callable) - 2 *sd(df_het_f$Callable)
df_het_f[df_het_f$Callable<nocall,]$Heterozygosity <- NA
df_het_2 <- df_het_f[complete.cases(df_het_f),]

g <- ggplot(df_het_2, aes(Heterozygosity, color=Sample))+
  geom_density(alpha=0.2)+ xlim(c(0,0.004))+ facet_wrap(~Sample)+
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_blank(),axis.line = element_line(colour = "black")) 

g

pdf("Plots/Density_Heterozygosity_windows.pdf",height = 8, width = 12)
g
dev.off()

# Median global heterozygosity for each window
list_het_mean <- list()
list_het_mean <- lapply(1:length(as.character(V1$V1)), function(i){
  df_sample<- df_het_2[df_het_2$Sample == as.character(V1$V1[i]),]
  cbind.data.frame(median=median(df_sample$Heterozygosity),sd=sd(df_sample$Heterozygosity),
                   Full_ID=samples$V1[i])
  
})

df_het_median <- do.call(rbind.data.frame,list_het_mean)
head(df_het_median)

df_het_median <- merge(df_het_median, metadata[ , 1:12], by="Full_ID")                            



p1 <- ggplot(df_het_median, aes(Type, median, fill=Type, color=Type)) +
  geom_violin(position = position_dodge(width = 1),width=1,color=NA, alpha=0.6)+ 
  geom_boxplot(width=0.3, position=position_dodge(width =1), color="black", outlier.shape=NA,lwd=0.3)+
  geom_jitter(width=0.1, size=0.8, color="black")  + ylab("Heterozygosity (bp-1)") + 
  theme_classic() + 
  theme(legend.position="none", 
        panel.border = element_blank(),#strip.text = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=13),
        axis.line = element_line(colour = "black")) + xlab("")

p1

mean(df_het_median$median) - 1*sd(df_het_median$median)


# ROH ---------
j=9
#for (j in 1:9){
j=9
roh <- read.table(paste0("Files/RoHs/BcfToolsRoh_alt04recomb_1e-",j,"_final_AF.txt")) 

roh$V9 <- roh$V2
samples <- read.table("Files/Samples", header = FALSE)

unique(roh$V2)
#size of roh
head(roh)
g <- ggplot(roh, aes(as.numeric(V6), color=V9))
g <- g +geom_density()+ geom_vline(xintercept = 10000)+
  theme(panel.grid.major = element_line(colour = "grey"),
        legend.position = "right",
        panel.background = element_blank(),axis.line = element_line(colour = "black")) + xlim(c(0,1000000))
g

assembly <- read.table("Files/assembly.fai")
colnames(assembly) <- c("Scaffold","Length","Pos","V4","V5")

assembly$V3 <- assembly$Scaffold
roh_final <- merge(roh, assembly[,c("V3","Scaffold","Length"),], by="V3")

chrom_want <- c("CM038391.1","CM038392.1","CM038393.1", "CM038394.1",
"CM038395.1", "CM038396.1" ,"CM038397.1", "CM038399.1","CM038400.1",
"CM038401.1", "CM038402.1", "CM038403.1","CM038404.1" ,"CM038405.1" ,"CM038406.1",
"CM038407.1","CM038408.1" ,"CM038409.1", "CM038410.1", "CM038411.1","CM038412.1", "CM038413.1")

roh_final$Scaffold <- factor(roh_final$Scaffold,levels=chrom_want, ordered = TRUE)

# density
g <- ggplot(roh_final, aes(V8, color=V2))
g <- g +geom_density()+
  theme(panel.grid.major = element_line(colour = "grey"),
        legend.position = "right",
        panel.background = element_blank(),axis.line = element_line(colour = "black")) 
g

roh_final <- roh_final[roh_final$V8 > 90,]
roh_final_SAVE <- roh_final
colnames(roh_final_SAVE) <- c("scaffold", "RG", "Sample", "Start","End","Length","Number of markers", "Quality","Library","Chrom")

#write.csv(roh_final_SAVE,"Files/FinalRoHs_newBed.csv", quote = FALSE)
samples_new <- unique(roh$V9)

# per sample  -----
j=9
for (i in 1:length(samples_new)){  
  sample <- as.character(samples_new[i])
  df <- roh_final[roh_final$V9==sample,]
  #df_het_sample <- df_het_2[df_het_2$Sample == sample,]
  p1 <- ggplot(df)
  p1 <- p1 + facet_wrap(.~Scaffold, ncol=1,strip.position="right")+ ggtitle(sample) + 
    ylab("")+ xlab("Position (bp)")+
    scale_x_continuous(name = "", waiver(), labels=comma, expand = c(0.001,0), 
                       limits = c(0,224379228))+
    scale_y_continuous(name = "", waiver(), labels=comma,
                       limits = c(0,0.008))+
    geom_rect(data=df,
              mapping=aes(xmin = V4, ymin = 0, xmax = V5, ymax = 0.008), 
              fill="#f56c42", alpha=0.8) + 
    geom_rect(data=assembly,mapping = aes(xmin=0, xmax=Length, ymin=0, ymax=0.008), alpha = .3,
              color = "#2e2e2e", fill = NA) +
   # geom_point(data= df_het_sample[df_het_sample$Heterozygosity>1.790398e-05,], aes(Pos,Heterozygosity), size=0.05, color="#b0b0b0", alpha=0.6) + 
    #geom_point(data=df_het_sample[df_het_sample$Heterozygosity<=1.790398e-05,], size=0.05, color="#afff6e", alpha=0.9, aes(Pos,Heterozygosity)) +
    
    theme_classic()+
    theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),panel.grid.minor = element_blank(),
          legend.position = "none", axis.ticks.y = element_blank(), panel.border=element_blank(),
          panel.background = element_blank(),axis.line = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank())
  print(p1)
  ggsave(paste0("Plots/RoHs/Final_GenomeinRoHs_", sample, "final90_rec",j,"_AF.pdf"), width = 6, height = 6)
}


#}
roh_final$Sample <- roh_final$V2

# per chrom  -----

j=9
for (i in 1:length(scaffolds_fin)){  
  chrom <- as.character(scaffolds_fin[i])
  df <- roh_final[roh_final$Scaffold==chrom,]
  df_het_sample <- df_het_2[df_het_2$Scaffold == chrom,]
  p1 <- ggplot(df)
  p1 <- p1 + facet_wrap(.~Sample, ncol=1,strip.position="top")+ ggtitle(chrom) + 
    ylab("")+ xlab("Position (bp)")+
    scale_x_continuous(name = "", waiver(), labels=comma, expand = c(0.001,0), 
                       limits = c(0,224379228))+
    scale_y_continuous(name = "", waiver(), labels=comma,
                       limits = c(0,0.008))+
    geom_rect(data=df,
              mapping=aes(xmin = V4, ymin = 0, xmax = V5, ymax = 0.008), 
              fill="#f56c42", alpha=0.8) + 
    geom_rect(data=assembly,mapping = aes(xmin=0, xmax=Length, ymin=0, ymax=0.008), alpha = .3,
              color = "#2e2e2e", fill = NA) +
    geom_point(data= df_het_sample[df_het_sample$Heterozygosity>1.790398e-05,], aes(Pos,Heterozygosity), size=0.05, color="#b0b0b0", alpha=0.6) + 
    geom_point(data=df_het_sample[df_het_sample$Heterozygosity<=1.790398e-05,], size=0.05, color="#afff6e", alpha=0.9, aes(Pos,Heterozygosity)) +
    
    theme_classic()+
    theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),panel.grid.minor = element_blank(),
          legend.position = "none", axis.ticks.y = element_blank(), panel.border=element_blank(),
          panel.background = element_blank(),axis.line = element_blank(), 
          strip.background = element_blank())
  print(p1)
  ggsave(paste0("Plots/RoHs/AllSampes_scaffold_", chrom, "final90_rec",j,"_alt0001_noAF.pdf"), width = 8, height = 16)
}

# % of genome in roh and separating by different sizes----
chrom_want
sum(assembly[which(assembly$Scaffold%in%chrom_want),]$Length)


roh_perc_bcftools <- lapply(1:length(samples_new), function(i) {
  df2 <- roh_final[roh_final$V2 == samples_new[i],]
  data.frame(sample=samples_new[i],
             Total=sum(df2$V6)/2779462738,  
             Short=sum(df2[df2$V6<=100000,]$V6)/2779462738, 
             Intermediate=sum(df2[df2$V6>100000 & df2$V6<=500000,]$V6)/2779462738, 
             Intermediate2=sum(df2[df2$V6>500000 & df2$V6<=1000000,]$V6)/2779462738, 
             Long=sum(df2[df2$V6>1000000 & df2$V6<=5000000,]$V6)/2779462738, 
             ExtraLong=sum(df2[df2$V6>5000000,]$V6)/2779462738, 
             count=length(df2$V6),
             sumMb=sum(df2$V6))
})
roh_perc_df_bcftools <- do.call(rbind, roh_perc_bcftools)


roh_perc_df_gather <- gather(roh_perc_df_bcftools, key, value, -sample, -count,-sumMb)
roh_perc_df_gather$key <- gsub("Short","<100Kb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate2","500Kb-1Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate","100-500Kb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("ExtraLong",">5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long","1-5Mb", roh_perc_df_gather$key )
order <- c("Total", "<100Kb","100-500Kb","500Kb-1Mb","1-5Mb",">5Mb")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$Full_ID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="Full_ID")

p <- ggplot(roh_perc_df_gather2, aes(key, value*100, fill=Site))
p <- p + geom_col(position = "dodge")+ theme_classic() + 
  facet_wrap(Site~Full_ID, labeller = label_wrap_gen(multi_line=FALSE)) + xlab("")+
  ylab("% of genome in RoHs") + theme(axis.text.x = element_text(angle=45, hjust=1),
                                      legend.position = "none") 
p

pdf("Plots/Genome_roh_size_indv_final.pdf", width = 6, height = 5)
p
dev.off()

# FROH ----
roh_perc_df_bcftools$Full_ID <- roh_perc_df_bcftools$sample
roh_perc_df_bcftools <- merge(roh_perc_df_bcftools, metadata, by="Full_ID")


p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$key=="Total",], aes(sample, value, color=Site))
p <- p + geom_segment( color="grey",aes(x=sample, xend=sample, y=0, yend=value)) +
  geom_point( size=5,alpha=0.7,stroke=2) +
  theme_classic() + facet_grid(~Type, space="free", scales="free") + xlab("")+
  ylab(expression(F[ROH])) + theme(axis.text.x = element_text(angle=45, hjust=1),
                                   legend.position = "none") 
p

pdf("Plots/FROH_indv_final.pdf", width = 6, height = 5)
p
dev.off()

p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$key=="Total" & roh_perc_df_gather2$Coverage>5 ,], aes(sample, value, color=Site))
p <- p + geom_segment( color="grey",aes(x=sample, xend=sample, y=0, yend=value)) +
  geom_point( size=5,alpha=0.7,stroke=2) +
  theme_classic() + facet_grid(~Type, space="free", scales="free") + xlab("")+
  ylab(expression(F[ROH])) + theme(axis.text.x = element_text(angle=45, hjust=1),
                                   legend.position = "none") 
p

pdf("Plots/FROH_indv_final_5x.pdf", width = 6, height = 5)
p
dev.off()

p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$key=="Total",], aes(F_plink, value, color=Site))
p <- p + 
  geom_point( size=5,alpha=0.7,stroke=2) +
  theme_classic() + 
  ylab(expression(F[ROH]))
p

pdf("Plots/FROH_indv_vs_Fplink.pdf", width = 6, height = 5)
p
dev.off()

p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$key=="Total"& roh_perc_df_gather2$Coverage>5 ,], aes(F_plink, value, color=Site))
p <- p + 
  geom_point( size=5,alpha=0.7,stroke=2) +
  theme_classic() + 
  ylab(expression(F[ROH]))
p

pdf("Plots/FROH_indv_vs_Fplink_5x.pdf", width = 6, height = 5)
p
dev.off()

library(ggpubr)
p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$Full_ID!="CEI_060",], aes(key, value, color=Type))
p <- p + geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge())+
  theme_classic() + xlab("")+
  ylab(expression(F[ROH])) + stat_compare_means(aes(group = Type))
p

pdf("Plots/FROH_comparison_means.pdf", width = 12, height = 5)
p
dev.off()

p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$Full_ID!="CEI_060"&roh_perc_df_gather2$Coverage>5,], aes(key, value, color=Type))
p <- p + geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge())+
  theme_classic() +  xlab("")+
  ylab(expression(F[ROH])) + stat_compare_means(aes(group = Type))
p

pdf("Plots/FROH_comparison_means_5x.pdf", width = 12, height = 5)
p
dev.off()

p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$Full_ID!="CEI_060"&roh_perc_df_gather2$Coverage>9,], aes(key, value, color=Type))
p <- p + geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge())+
  theme_classic() +  xlab("")+
  ylab(expression(F[ROH])) + stat_compare_means(aes(group = Type))
p

pdf("Plots/FROH_comparison_means_10x.pdf", width = 12, height = 5)
p
dev.off()

p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$Site=="San Juán",], aes(key, value, color=Type))
p <- p + geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge())+
  theme_classic() +  xlab("")+
  ylab(expression(F[ROH])) + stat_compare_means(aes(group = Type))
p

pdf("Plots/FROH_comparison_means_Sanjuan.pdf", width = 12, height = 5)
p
dev.off()

p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$Full_ID!="CEI_060"&roh_perc_df_gather2$Coverage>5,], aes(Site, value, color=Type))
p <- p + geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge())+
  theme_classic() +  xlab("")+ facet_grid(.~key, space="free", scales="free_x")+
  ylab(expression(F[ROH]))
p



p <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$key=="Total",], aes(Coverage, value, color=Full_ID))
p <- p +
  geom_point( size=5,alpha=0.7,stroke=2) +
  theme_classic() +
  ylab(expression(F[ROH])) + theme(axis.text.x = element_text(angle=45, hjust=1),
                                   legend.position = "none") 
p
pdf("Plots/FROH_indv_coverage.pdf", width = 6, height = 5)
p
dev.off()

# plot cumulative rohs versos counts
library(ggConvexHull)
library(ggpubr)
p <- ggplot(roh_perc_df_bcftools[roh_perc_df_bcftools$Full_ID!="CEI_060",], aes(sumMb, count, color=Type, fill=Type))
p <- p + geom_point(size=3)+ theme_classic()  + xlab("")+ ylab(expression(F[ROH]))+
  geom_convexhull(alpha = 0.3)

p

pdf("Plots/Genome_Roh_TestType.pdf", width = 5, height = 4)
p
dev.off()

p <- ggplot(roh_perc_df_bcftools[roh_perc_df_bcftools$Coverage>5&roh_perc_df_bcftools$Full_ID!="CEI_060",], aes(sumMb, count, color=Type,
                                                                                                                fill=Type))
p <- p + geom_point(size=3)+ theme_classic()  + xlab("sROH")+ ylab("nROH") + 
  geom_convexhull(alpha = 0.3)
p

pdf("Plots/Genome_sroh_roh_final_5x.pdf", width = 5, height = 4)
p
dev.off()


# Test between type modern and historical and by site
# average per type
df_new_historical <- roh_perc_df_bcftools[roh_perc_df_bcftools$Type=="Historical",]
df_new_modern <- roh_perc_df_bcftools[roh_perc_df_bcftools$Type=="Modern"&roh_perc_df_bcftools$Full_ID!="CEI_060",]

df_new_historical_average<-cbind.data.frame(Type="Historical",Mean=colMeans(x=df_new_historical[,3:8]),SD=apply(df_new_historical[,3:8],2,sd))
df_new_historical_average$RoHs <- row.names(df_new_historical_average)
df_new_modern_average<-cbind.data.frame(Type="Modern",Mean=colMeans(x=df_new_modern[,3:8]),SD=apply(df_new_modern[,3:8],2,sd))
df_new_modern_average$RoHs <- row.names(df_new_modern_average)

df_new_average <- rbind.data.frame(df_new_historical_average,df_new_modern_average)

df_new_average$RoHs <- gsub("Short","<100Kb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Intermediate2","500Kb-1Mb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Intermediate","100-500Kb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("ExtraLong",">5Mb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Long","1-5Mb", df_new_average$RoHs )

order <- c("Total", "<100Kb","100-500Kb","500Kb-1Mb","1-5Mb",">5Mb")

df_new_average$RoHs <- factor(df_new_average$RoHs, levels = order, ordered = TRUE)

p <- ggplot(df_new_average, aes(RoHs, Mean*100, fill=Type))
p <- p + geom_col(position = "dodge")+ geom_errorbar(aes(ymin= ifelse(Mean*100 - SD*100 < 0, 0,  Mean*100 - SD*100),
                                                         ymax=Mean*100+SD*100), width=.2,
                                                     position=position_dodge(.9)) +
  theme_classic()  + xlab("")+ 
  ylab("% of genome in RoHs") + theme(legend.position = "top", legend.title = element_blank()) + 
  scale_fill_manual(values=c("lightgrey","#9a0026","#7293e8"))
p

pdf("Plots/Genome_roh_size_average_final.pdf", width = 6, height = 5)
p
dev.off()

#only >5x
df_new_historical <- roh_perc_df_bcftools[roh_perc_df_bcftools$Type=="Historical"&roh_perc_df_bcftools$Coverage>5,]
df_new_modern <- roh_perc_df_bcftools[roh_perc_df_bcftools$Type=="Modern"&roh_perc_df_bcftools$Full_ID!="CEI_060"&roh_perc_df_bcftools$Coverage>5,]

df_new_historical_average<-cbind.data.frame(Type="Historical",Mean=colMeans(x=df_new_historical[,3:8]),SD=apply(df_new_historical[,3:8],2,sd))
df_new_historical_average$RoHs <- row.names(df_new_historical_average)
df_new_modern_average<-cbind.data.frame(Type="Modern",Mean=colMeans(x=df_new_modern[,3:8]),SD=apply(df_new_modern[,3:8],2,sd))
df_new_modern_average$RoHs <- row.names(df_new_modern_average)

df_new_average <- rbind.data.frame(df_new_historical_average,df_new_modern_average)

df_new_average$RoHs <- gsub("Short","<100Kb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Intermediate2","500Kb-1Mb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Intermediate","100-500Kb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("ExtraLong",">5Mb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Long","1-5Mb", df_new_average$RoHs )

order <- c("Total", "<100Kb","100-500Kb","500Kb-1Mb","1-5Mb",">5Mb")

df_new_average$RoHs <- factor(df_new_average$RoHs, levels = order, ordered = TRUE)


p <- ggplot(df_new_average, aes(RoHs, Mean*100, fill=Type))
p <- p + geom_col(position = "dodge")+ geom_errorbar(aes(ymin= ifelse(Mean*100 - SD*100 < 0, 0,  Mean*100 - SD*100),
                                                         ymax=Mean*100+SD*100), width=.2,
                                                     position=position_dodge(.9)) +
  theme_classic()  + xlab("")+ 
  ylab("% of genome in RoHs") + theme(legend.position = "top", legend.title = element_blank()) + 
  scale_fill_manual(values=c("lightgrey","#9a0026","#7293e8"))
p

pdf("Plots/Genome_roh_size_average_final_5x.pdf", width = 6, height = 5)
p
dev.off()

#only San Juan
df_new_historical <- roh_perc_df_bcftools[roh_perc_df_bcftools$Type=="Historical"&roh_perc_df_bcftools$Site>"San Juan",]
df_new_modern <- roh_perc_df_bcftools[roh_perc_df_bcftools$Type=="Modern"&roh_perc_df_bcftools$Site>"San Juan",]

df_new_historical_average<-cbind.data.frame(Type="Historical",Mean=colMeans(x=df_new_historical[,3:8]),SD=apply(df_new_historical[,3:8],2,sd))
df_new_historical_average$RoHs <- row.names(df_new_historical_average)
df_new_modern_average<-cbind.data.frame(Type="Modern",Mean=colMeans(x=df_new_modern[,3:8]),SD=apply(df_new_modern[,3:8],2,sd))
df_new_modern_average$RoHs <- row.names(df_new_modern_average)

df_new_average <- rbind.data.frame(df_new_historical_average,df_new_modern_average)

df_new_average$RoHs <- gsub("Short","<100Kb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Intermediate2","500Kb-1Mb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Intermediate","100-500Kb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("ExtraLong",">5Mb", df_new_average$RoHs )
df_new_average$RoHs <- gsub("Long","1-5Mb", df_new_average$RoHs )

order <- c("Total", "<100Kb","100-500Kb","500Kb-1Mb","1-5Mb",">5Mb")

df_new_average$RoHs <- factor(df_new_average$RoHs, levels = order, ordered = TRUE)


p <- ggplot(df_new_average, aes(RoHs, Mean*100, fill=Type))
p <- p + geom_col(position = "dodge")+ geom_errorbar(aes(ymin= ifelse(Mean*100 - SD*100 < 0, 0,  Mean*100 - SD*100),
                                                         ymax=Mean*100+SD*100), width=.2,
                                                     position=position_dodge(.9)) +
  theme_classic()  + xlab("")+ 
  ylab("% of genome in RoHs") + theme(legend.position = "top", legend.title = element_blank()) + 
  scale_fill_manual(values=c("lightgrey","#9a0026","#7293e8"))
p

pdf("Plots/Genome_roh_size_average_final_SanJuan.pdf", width = 6, height = 5)
p
dev.off()


# RoHs per individual
roh_perc_df_bcftools_subset <- cbind.data.frame(roh_perc_df_bcftools$Full_ID,roh_perc_df_bcftools$Short,
                                                roh_perc_df_bcftools$Intermediate, roh_perc_df_bcftools$Intermediate2,
                                                roh_perc_df_bcftools$Long, roh_perc_df_bcftools$ExtraLong, roh_perc_df_bcftools$Type,
                                                roh_perc_df_bcftools$Site, roh_perc_df_bcftools$Coverage)

colnames(roh_perc_df_bcftools_subset) <- c("Full_ID","Short","Intermediate", "Intermediate2","Long","ExtraLong","Type","Site", "Coverage")
roh_new2 <-  gather(roh_perc_df_bcftools_subset, key, value, -Full_ID, -Site, -Type, -Coverage)

roh_new2$key <- gsub("Short","<100Kb", roh_new2$key )
roh_new2$key <- gsub("Intermediate2","500Kb-1Mb", roh_new2$key )
roh_new2$key <- gsub("Intermediate","100-500Kb", roh_new2$key )
roh_new2$key <- gsub("ExtraLong",">5Mb", roh_new2$key )
roh_new2$key <- gsub("Long","1-5Mb", roh_new2$key )

order <- c("<100Kb","100-500Kb","500Kb-1Mb","1-5Mb",">5Mb")
roh_new2$key <- factor(roh_new2$key, levels=order, ordered = TRUE)

pa <- ggplot(roh_new2, aes(value*100, Full_ID,fill=key))
pa <- pa + geom_col()+ theme_classic() +  ylab("")+ facet_grid(Type~.,space="free", scales="free")+
  xlab("% of Genome in RoH") + 
  theme(legend.position = "right", strip.background = element_blank()) + 
  scale_fill_manual(values=rev(c("#fcd021","#234f96","#5e97f2","#b6d1fc","grey")), name="RoH Size") 
pa
pdf("Plots/ROHs_individual_size.pdf", height = 8, width = 6)
pa
dev.off()
pa <- ggplot(roh_new2[roh_new2$Coverage>5,], aes(value*100, Full_ID,fill=key))
pa <- pa + geom_col()+ theme_classic() +  ylab("")+ facet_grid(Type~.,space="free", scales="free")+
  xlab("% of Genome in RoH") + 
  theme(legend.position = "right", strip.background = element_blank()) + 
  scale_fill_manual(values=rev(c("#fcd021","#234f96","#5e97f2","#b6d1fc","grey")), name="RoH Size") 
pa
pdf("Plots/ROHs_individual_size_5x.pdf", height = 8, width = 6)
pa
dev.off()
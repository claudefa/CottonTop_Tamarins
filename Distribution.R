#Distribution of quality and reath depth for genetic load analysis
library("ggplot2")
library("readxl")
library(ggpubr)

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Revision/")
metadata <- read_excel("../Paper/RevisionMolEcol/SupplementaryTables_v2.xlsx", sheet = 4)


# annotated positions High moderate and synonymous
fins <- list.files("Files/Category_v2/",pattern=".txt",full.names=TRUE)
samples <- sub('\\_v2.txt$', '', basename(fins))
samples <- gsub("DP_GQ_","", samples)
samples <- gsub("_HIGH","", samples)
samples <- gsub("_MODERATE","", samples)
samples <- gsub("_SYNONYM","", samples)

type <- rep(c("High","Moderate","Synonym"), times=15)

df_list <- lapply(fins,FUN = function(files) {read.table(files)})
df_list2 <- list()
for (i in 1:length(samples)){
  df_list2[[i]] <- cbind.data.frame(Sample=samples[i], Data=df_list[[i]],
                                    Type=type[i])
}

qual_df <- do.call(rbind,df_list2)
colnames(qual_df) <- c("Full_ID","Genotype","DP","GQ","Type")

qual_df_M <- merge(qual_df, metadata[ , 1:13], by="Full_ID")                            
colnames(qual_df_M) <- c("Full_ID"  
                         ,"Genotype"                         
                         ,"DP"                               
                         ,"GQ"                               
                         ,"Type_Mut"                           
                         ,"Type_Sample"                           
                         ,"Collection Year"                  
                         ,"Sample Material"                  
                         ,"Department"                       
                         , "Municipality"                     
                         , "Site"                             
                         , "Latitude"                         
                         , "Longitude"                        
                         , "Sex"                              
                         , "Coverage"                         
                        , "Coverage Subset ~5x (downsampled)","Group")


qual_df_M$Group <- gsub(" - Uncertain","",qual_df_M$Group  )

g <- ggplot(qual_df_M[qual_df_M$Genotype!="./.",],aes(DP, color=Full_ID))
g <- g+ geom_density() + facet_grid(Type_Mut~Type_Sample) + theme_classic() + ggtitle("Depth per site") +
  xlim(c(0,100))+ geom_vline(xintercept = 5)
g


pdf("DP_min5.pdf", height = 8, width = 10)
g
dev.off()

g <- ggplot(qual_df_M[qual_df_M$Genotype=="0/0",],aes(DP, color=Full_ID))
g <- g+ geom_density() + facet_grid(Type_Mut~Type_Sample) + theme_classic() + ggtitle("Depth per site") +
  xlim(c(0,100)) + geom_vline(xintercept = 5)

g

g <- ggplot(qual_df_M[qual_df_M$Genotype=="1/1",],aes(DP, color=Full_ID))
g <- g+ geom_density() + facet_grid(Type_Mut~Type_Sample) + theme_classic() + ggtitle("Depth per site") +
  xlim(c(0,100)) + geom_vline(xintercept = 5)

g
pdf("DP_all.pdf", height = 8, width = 10)
g
dev.off()

g <- ggplot(qual_df_M[qual_df_M$Genotype)="./.",],aes(GQ, color=Full_ID))
g <- g+ geom_density() + facet_grid(Type_Mut~Type_Sample) + theme_classic() + ggtitle("GQ per site") 

g
pdf("GQ_min5.pdf", height = 8, width = 10)
g
dev.off()

g <- ggplot(qual_df_M,aes(GQ, color=Full_ID))
g <- g+ geom_density() + facet_grid(Type_Mut~Type_Sample) + theme_classic() + ggtitle("GQ per site") 

g
pdf("GQ_all.pdf", height = 8, width = 10)
g
dev.off()

#Boxplots

#Boxplots DP----
g <- ggplot(qual_df_M,aes(Type_Mut, DP, group= interaction(Type_Mut,Type_Sample)))
g <- g+ geom_violin()+ 
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  +
  theme_classic() + 
  ggtitle("Depth per site (all sites)")+ #ylim(c(0,50))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test") + xlab("")+ theme(legend.position = "top")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))
g

g1 <- ggplot(qual_df_M[qual_df_M$Genotype!="./.",],aes(Type_Mut, DP, group= interaction(Type_Mut,Type_Sample)))
g1 <- g1+geom_violin()+
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  + 
  theme_classic() + ggtitle("Depth per site (DP >= 5 & GQ >=30)")+
 # ylim(c(0,50))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test")+  theme(legend.position = "top")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))+ xlab("")

g1

g11 <- ggplot(qual_df_M[qual_df_M$Genotype=="0/0",],aes(Type_Mut, DP, group= interaction(Type_Mut,Type_Sample)))
g11 <- g11+geom_violin()+
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  + 
  theme_classic() + ggtitle("Depth per site 0/0 (DP >= 5 & GQ >=30)")+
  # ylim(c(0,50))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test")+  theme(legend.position = "top")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))+ xlab("")

g11

g2 <- ggplot(qual_df_M[qual_df_M$Genotype=="0/1",],aes(Type_Mut, DP, group= interaction(Type_Mut,Type_Sample)))
g2 <- g2+geom_violin()+ 
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  + 
  theme_classic() + ggtitle("Depth per site 0/1 (DP > 5 & GQ >=30)")+
#  ylim(c(0,50))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test")+ theme(legend.position = "top")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF")) + xlab("")

g2
g3 <- ggplot(qual_df_M[qual_df_M$Genotype=="1/1",],aes(Type_Mut, DP, group= interaction(Type_Mut,Type_Sample)))
g3 <- g3+geom_violin()+
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  + 
  theme_classic() + ggtitle("Depth per site 1/1 (DP > 5  & GQ >=30)")+ theme(legend.position = "top")+
#  ylim(c(0,50))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF")) + xlab("")

g3

ggarrange(g11,g2,g3 , labels = c("A", "B","C","D"), ncol=2, nrow = 2)
ggsave("Distribution_DPt.pdf", height = 13, width = 23) 

# Boxplots genotype quality -----
g <- ggplot(qual_df_M,aes(Type_Mut, GQ, group= interaction(Type_Mut,Type_Sample)))
g <- g+ geom_violin()+ 
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  +
  theme_classic() + 
  ggtitle("Genotype Quality per site (all sites)")+ 
  #ylim(c(0,100))+ #facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test") + xlab("")+ theme(legend.position = "top")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))
g

g1 <- ggplot(qual_df_M[qual_df_M$Genotype!="./.",],aes(Type_Mut, GQ, group= interaction(Type_Mut,Type_Sample)))
g1 <- g1+geom_violin()+
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  + 
  theme_classic() + ggtitle("Genotype Quality per site (DP > 5 & GQ >=30)")+
 # ylim(c(25,100))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test")+  theme(legend.position = "top")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))+ xlab("")

g1
p11 <- ggplot(qual_df_M[qual_df_M$Genotype=="0/0",],aes(Type_Mut, GQ, group= interaction(Type_Mut,Type_Sample)))
p11 <- p11+geom_violin()+
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  + 
  theme_classic() + ggtitle("Genotype Quality per site 0/0 (DP >= 5 & GQ >=30)")+
  # ylim(c(0,50))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test")+  theme(legend.position = "top")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF"))+ xlab("")

p11

p2 <- ggplot(qual_df_M[qual_df_M$Genotype=="0/1",],aes(Type_Mut, GQ, group= interaction(Type_Mut,Type_Sample)))
p2 <- p2+geom_violin()+ 
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  + 
  theme_classic() + ggtitle("Genotype Quality per site 0/1 (DP > 5 & GQ >=30)")+
 # ylim(c(25,100))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test")+ theme(legend.position = "top")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF")) + xlab("")

p2
p3 <- ggplot(qual_df_M[qual_df_M$Genotype=="1/1",],aes(Type_Mut, GQ, group= interaction(Type_Mut,Type_Sample)))
p3 <- p3+geom_violin()+
  geom_boxplot(width=0.5,position = position_dodge(width=0.8),outlier.shape = NA,aes(fill=Type_Sample))  + 
  theme_classic() + ggtitle("Genotype Quality per site 1/1 (DP > 5 & GQ >=30)")+ theme(legend.position = "top")+
  #ylim(c(25,100))+ 
  facet_grid(.~Group)+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values=c("#913C5C","#FEDFFF")) + xlab("")

p3

ggarrange(g11,g2,g3 , labels = c("D", "E","F"), ncol=2, nrow = 2)
ggsave("Distribution_GQ.pdf", height = 13, width = 23) 

ggarrange(g11,g2,g3,p11,p2,p3, labels = c("A","B","C","D", "E","F"), ncol=2, nrow = 3)
ggsave("Revision/Distribution_DP_GQ.pdf", height = 18, width = 21) 


mean(qual_df_M[qual_df_M$Genotype=="0/0"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Historical",]$GQ)
mean(qual_df_M[qual_df_M$Genotype=="0/0"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Modern",]$GQ)

mean(qual_df_M[qual_df_M$Genotype=="0/0"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Historical",]$DP)
mean(qual_df_M[qual_df_M$Genotype=="0/0"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Modern",]$DP)

mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Historical",]$GQ)
mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Modern",]$GQ)

mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Historical",]$DP)
mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Modern",]$DP)

mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Historical",]$GQ)
mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Modern",]$GQ)

mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Historical",]$DP)
mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="High"&qual_df_M$Type_Sample=="Modern",]$DP)


mean(qual_df_M[qual_df_M$Genotype=="0/0"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Historical",]$GQ)
mean(qual_df_M[qual_df_M$Genotype=="0/0"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Modern",]$GQ)

mean(qual_df_M[qual_df_M$Genotype=="0/0"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Historical",]$DP)
mean(qual_df_M[qual_df_M$Genotype=="0/0"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Modern",]$DP)

mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Historical",]$GQ)
mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Modern",]$GQ)

mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Historical",]$DP)
mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Modern",]$DP)

mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Historical",]$GQ)
mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Modern",]$GQ)

mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Historical",]$DP)
mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="Moderate"&qual_df_M$Type_Sample=="Modern",]$DP)



mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="Synonym"&qual_df_M$Type_Sample=="Historical",]$GQ)
mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="Synonym"&qual_df_M$Type_Sample=="Modern",]$GQ)

mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="Synonym"&qual_df_M$Type_Sample=="Historical",]$DP)
mean(qual_df_M[qual_df_M$Genotype=="1/1"&qual_df_M$Type_Mut=="Synonym"&qual_df_M$Type_Sample=="Modern",]$DP)

mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="Synonym"&qual_df_M$Type_Sample=="Historical",]$GQ)
mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="Synonym"&qual_df_M$Type_Sample=="Modern",]$GQ)

mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="Synonym"&qual_df_M$Type_Sample=="Historical",]$DP)
mean(qual_df_M[qual_df_M$Genotype=="0/1"&qual_df_M$Type_Mut=="Synonym"&qual_df_M$Type_Sample=="Modern",]$DP)

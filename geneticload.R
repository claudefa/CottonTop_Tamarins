#Genetic Load test
library(ggplot2)
library(readxl)
library(tidyr)
library(ggpubr)

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
df<- read_excel("Paper/SupplementaryTables_v1.1.xlsx", sheet = 7) 

colnames(df) <- c("Sample","Type","Grouping","Load","RealizedCount",
                  "MaskedCount","TotalCount","SynonymousCount" ,
                  "Total Load","Realized Load","Masked Load")
                  
df_tidy <- gather(data = df, key = GeneticLoad, value = RelativeCount, -Sample, -Type, -Grouping,
                  -RealizedCount, -TotalCount, -MaskedCount, -Load, -SynonymousCount)

# Samples > 5x --------
g1 <- ggplot(df_tidy, aes(GeneticLoad,RelativeCount, fill=Type))
g1 <- g1+ geom_boxplot(outlier.shape = NA) +   geom_point(position=position_jitterdodge(), size=1, color="grey9")+
  facet_grid(Load~Grouping, scales="free_y")+   xlab("") + 
  scale_fill_manual(values=c("#913C5C","#FEDFFF")) +   theme_classic() +
   ylab("Deleterious Allele Count / Synonymous")+ 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        legend.title = element_blank(),
        panel.grid.major = element_line(colour="grey", size=0.1),
        panel.border = element_rect(colour="grey", size=0.5, fill = NA),
        legend.position = "top")
g1
pdf("Plots/GeneticLoad_bySampleType.pdf", height = 7, width = 10)
g1
dev.off()


# comparisons for Ratio Genetic Load
df_comp <- data.frame(Comparisons=rep(c("Northeast: Modern/Historical",
                                    "Southwest: Modern/Historical",
                                    "Historical: Northeast/Southwest",
                                    "Modern: Northeast/Southwest"), times=2),
                      Type=rep(c("High","Moderate"),each=4))

df_tidy_final <- df_tidy[df_tidy$GeneticLoad=="Total Load",]
x11 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="High"&
                df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="High"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)
x12 <- mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="High"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="High"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)

x13 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="High"&
                            df_tidy_final$Type=="Historical",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="High"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)

x14 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="High"&
                              df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="High"&
                       df_tidy_final$Type=="Modern",]$RelativeCount)


x21 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="Moderate"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="Moderate"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)
x22 <- mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="Moderate"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="Moderate"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)

x23 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="Moderate"&
                            df_tidy_final$Type=="Historical",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="Moderate"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)

x24 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="Moderate"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="Moderate"&
                       df_tidy_final$Type=="Modern",]$RelativeCount)

df_comp$Values <- c(x11,x12,x13,x14,x21,x22,x23,x24)

g <- ggplot(df_comp, aes(Comparisons, Values-1, fill=Type))
g <- g + geom_col(position = "dodge") + theme_classic() + ylab("Total Load Ratio") +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.title = element_blank(),
        legend.position = c(0.9,0.9))+ xlab("")+
  scale_y_continuous(labels = function(y) y + 1) +
  scale_fill_manual(values=c("#ff875c","#1a659e"))  
g

pdf("Plots/RatioLoad.pdf", width = 6, height = 4)
g
dev.off()

ggarrange(ggarrange(g1, labels=c("A"), ncol=1),
          ggarrange(NULL,g,NULL, ncol=3, widths = c(1,2,1), 
                                     labels = c("","B","")), heights = c(1.5,1), ncol=1)
ggsave("Plots/Figure4.pdf", width = 8, height = 10)


# Only samples with > 7x coverage ----
highcov <- c("24_FMNH_69286","36_AMNH_32702","39_AMNH_130390","41_FMNH_69938","CEI_051","SJN_001","SJN_011","TLP2_03","TLP9_01")
df_tidy <- df_tidy[which(df_tidy$Sample%in%highcov),]

g1 <- ggplot(df_tidy, aes(GeneticLoad,RelativeCount, fill=Type))
g1 <- g1+ geom_boxplot(outlier.shape = NA) +   geom_point(position=position_jitterdodge(), size=1, color="grey9")+
  facet_grid(Load~Grouping, scales="free_y")+   xlab("") + 
  scale_fill_manual(values=c("#913C5C","#FEDFFF")) +   theme_classic() +
  ylab("Deleterious Allele Count / Synonymous")+ 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        legend.title = element_blank(),
        panel.grid.major = element_line(colour="grey", size=0.1),
        panel.border = element_rect(colour="grey", size=0.5, fill = NA),
        legend.position = "top")
g1
pdf("Plots/GeneticLoad_bySampleType_7x.pdf", height = 7, width = 10)
g1
dev.off()


df_tidy_final <- df_tidy[df_tidy$GeneticLoad=="Total Load",]
x11 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="High"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="High"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)
x12 <- mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="High"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="High"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)

x13 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="High"&
                            df_tidy_final$Type=="Historical",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="High"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)

x14 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="High"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="High"&
                       df_tidy_final$Type=="Modern",]$RelativeCount)


x21 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="Moderate"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="Moderate"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)
x22 <- mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="Moderate"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="Moderate"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)

x23 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="Moderate"&
                            df_tidy_final$Type=="Historical",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="Moderate"&
                       df_tidy_final$Type=="Historical",]$RelativeCount)

x24 <- mean(df_tidy_final[df_tidy_final$Grouping=="Greater Northeast"&df_tidy_final$Load=="Moderate"&
                            df_tidy_final$Type=="Modern",]$RelativeCount)/
  mean(df_tidy_final[df_tidy_final$Grouping=="Southwest"&df_tidy_final$Load=="Moderate"&
                       df_tidy_final$Type=="Modern",]$RelativeCount)

df_comp$Values <- c(x11,x12,x13,x14,x21,x22,x23,x24)

g <- ggplot(df_comp, aes(Comparisons, Values-1, fill=Type))
g <- g + geom_col(position = "dodge") + theme_classic() + ylab("Total Load Ratio") +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.title = element_blank(),
        legend.position = c(0.9,0.9))+ xlab("")+
  scale_y_continuous(labels = function(y) y + 1) +
  scale_fill_manual(values=c("#ff875c","#1a659e"))  
g

pdf("Plots/RatioLoad_7x.pdf", width = 6, height = 4)
g
dev.off()

ggarrange(ggarrange(g1, labels=c("A"), ncol=1),
          ggarrange(NULL,g,NULL, ncol=3, widths = c(1,2,1), 
                    labels = c("","B","")), heights = c(1.5,1), ncol=1)
ggsave("Plots/FigureS12.pdf", width = 8, height = 10)

#Genetic Load test
library(ggplot2)
library(readxl)
library(tidyr)
library(ggpubr)

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
df<- read_excel("CountsGeneticLoad_paper.xlsx", sheet = 1)


df_tidy <- gather(data = df, key = GeneticLoad, value = RelativeCount, -Sample, -Type, -Grouping,
                 -Realized, -Total, -Masked, -Load,-TotalSynonymous)

colnames(df_tidy)
df_tidy$GeneticLoad <- gsub("Masked/Synonymous", "Masked Load",df_tidy$GeneticLoad )
df_tidy$GeneticLoad <- gsub("Realized/Synonymous", "Realized Load",df_tidy$GeneticLoad )
df_tidy$GeneticLoad <- gsub("Total/Synonymous", "Total Load",df_tidy$GeneticLoad )

g <- ggplot(df_tidy, aes(GeneticLoad,RelativeCount, fill=Grouping))
g <- g+ geom_boxplot(outlier.shape = NA) +   geom_point(position=position_jitterdodge())+
  facet_grid(Load~Type, scales="free_y")+   xlab("") + 
   scale_fill_manual(values=c("#a6cee3","#1f78b4")) + 
  theme_light()+ ylab("Relative Count of Genetic Load") + theme(axis.text.x = element_text(angle=45, hjust=1))
g
pdf("Plots/GeneticLoad_bygrouping.pdf", height = 7, width = 10)
g
dev.off()

# Or this one
g1 <- ggplot(df_tidy, aes(GeneticLoad,RelativeCount, fill=Type))
g1 <- g1+ geom_boxplot(outlier.shape = NA) +   geom_point(position=position_jitterdodge())+
  facet_grid(Load~Grouping, scales="free_y")+   xlab("") + 
  scale_fill_manual(values=c("#913C5C","#FEDFFF")) + 
  theme_light()+ ylab("Relative Count of Genetic Load")+ theme(axis.text.x = element_text(angle=45, hjust=1),
                                                               legend.title = element_blank(),
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
        legend.position = c(0.9,0.8))+ xlab("")+
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

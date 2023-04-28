# AB per sample
library(tidyr)
library(dplyr)
library(ggplot2)

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")

fins <- list.files("Files/AB/",pattern=".txt",full.names=TRUE)
samples <- sub('\\.txt$', '', basename(fins))

df_list <- lapply(fins,
                  FUN = function(files) {
                    read.table(files)
                  })
df_list2 <- list()
for (i in 1:length(samples)){
  df1 <- df_list[[i]]
  df_list[[i]] <- cbind.data.frame(Sample=samples[i], 
                                   new1=sapply(strsplit(df1$V2, "[ ,]+"), function(i) sum(as.numeric(i)))/df1$V1,
                                   new2=sapply(strsplit(df1$V3, "[ ,]+"), function(i) sum(as.numeric(i)))/df1$V1,
                                   new3=sapply(strsplit(df1$V4, "[ ,]+"), function(i) sum(as.numeric(i)))/df1$V1,
                                   new4=sapply(strsplit(df1$V5, "[ ,]+"), function(i) sum(as.numeric(i)))/df1$V1
                                 )
}

df_list2[[i]] <- df_list[[i]] %>%
  mutate(across(-1, na_if, 0)) %>%
  tidyr::unite(final, new1:new4, na.rm = TRUE, sep = ';', remove = FALSE)


ab_df <- do.call(rbind,df_list)


df2 <- ab_df %>%
  mutate(across(-1, na_if, 0)) %>%
  tidyr::unite(final, new1:new4, na.rm = TRUE, sep = ';', remove = FALSE)

df2$AB<- as.numeric(unlist(lapply(1:length(df2$final), function(i) 
                    strsplit(df2$final[i],";")[[1]][1])))



#write.csv2(df2,file = "Files/AB.csv", quote = FALSE, row.names = TRUE)
#df2 <- read.csv2(file = "Files/AB.csv", sep=";")


df2_final <- cbind.data.frame(df2$Sample, df2$AB)
colnames(df2_final) <- c("Sample", "AB")
df2_final$Sample <- gsub("AB_", "", df2_final$Sample)

g <- ggplot(df2_final, aes(AB))
g <- g+ geom_density(adjust = 2) + facet_wrap(.~Sample) + theme_minimal()
g

pdf("Plots/AB_allsamples.pdf", height = 10, width = 10)
g
dev.off()



# Dissimilarity tree
library(Imap)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(ape)
library(corrplot)
setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")

# matrix
gen_dist_all <- read.table("Files/Distance_allsamples.txt", row.names = 1) 
colnames(gen_dist_all) <-  row.names(gen_dist_all)

indTable <- row.names(gen_dist_all)

gen_dist_all2 <- gen_dist_all
gen_dist_all2$Samples <- row.names(gen_dist_all2)

gen_dist_all3 <- gather(gen_dist_all2, Samples2, Dissimilarity, -Samples)

g <- ggplot(gen_dist_all3, aes(Samples2, Samples, fill=Dissimilarity))
g <- g + geom_tile() +   scale_fill_gradientn(colours = colorRampPalette(c(brewer.pal( 11, "RdYlBu")))(200)) +   
  ylab("") + xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1))
g

# remove CEI_060
gen_dist_all3_no <- gen_dist_all3[gen_dist_all3$Samples!= "CEI_060"& gen_dist_all3$Samples2!= "CEI_060",]


x <- spread(gen_dist_all3_no, Samples2, Dissimilarity, fill=NA,convert = FALSE)
names <- x$Samples
x <- x[,2:ncol(x)]
rownames(x) <- names
y <- data.matrix(x)
y

clustering <- hclust(as.dist(y))

pdf("Plots/Dendogram_noCEI_060.pdf", width = 10, height = 8)
plot(clustering)
dev.off()

library(ape)
class(clustering) # must be hclust class
my_tree <- as.phylo(clustering) 
write.tree(phy=my_tree, file="Dendogram.newick") # look for the file in your working directory



pdf("Plots/Dissimilarity.pdf", height=8, width = 10 )
corrplot(y,addrect = 4, na.label = "NA",
         tl.col = "black", type="full",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90, 
         rect.col="black",
         pch.cex = 1.0, method="color",is.corr = FALSE, col=colorRampPalette(c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(200),
         pch.col = "purple")

dev.off()


#without outgroup
gen_dist_all3_no <- gen_dist_all3[gen_dist_all3$Samples!= "37_AMNH_33174"& gen_dist_all3$Samples2!= "37_AMNH_33174",]
x <- spread(gen_dist_all3_no, Samples2, Dissimilarity, fill=NA,convert = FALSE)
names <- x$Samples
x <- x[,2:ncol(x)]
rownames(x) <- names
y <- data.matrix(x)
y
clustering <- hclust(as.dist(y))
pdf("Plots/Dendogram_no.pdf", width = 10, height = 8)
plot(clustering)
dev.off()


pdf("Plots/Dissimilarity_no.pdf", height=8, width = 10 )
corrplot(y,addrect = 4, na.label = "NA",
         tl.col = "black", type="full",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90, 
         rect.col="black",
         pch.cex = 1.0, method="color",is.corr = FALSE, col=colorRampPalette(c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(200),
         pch.col = "purple")

dev.off()


# make tree 
tree <- read.tree("Files/all.raxml.support")

tree <- root(tree, outgroup = "37_AMNH_33174",resolve.root = TRUE)

#pdf("Tree_distance.pdf", height = 20, width = 20)
plot(tree, type = "r", show.tip.lab = TRUE, show.node.labe=TRUE,font=1, cex=0.9)
plot(tree, type = "phylogram", show.tip.lab = TRUE, show.node.labe=TRUE,font=1, cex=0.9,no.margin=T)
plot(tree, type = "tidy", show.tip.lab = TRUE, show.node.labe=TRUE,font=1, cex=0.9,no.margin=T)

#dev.off()

#COLS = rep(c("#0085DD","#DB3670","#FF7F00","#33A02C","black"), times=c(109,16,60,26,1))
#plotTree(tree,tip.color=COLS)

#pdf("Tree_all_nooutliers_genotype2.pdf", height = 18, width = 15)
#plot(tree, type = "r", show.tip.lab = TRUE, show.node.labe=TRUE,font=1, cex=0.7, tip.color=COLS)
#plot(tree, type = "phylogram", show.tip.lab = TRUE,show.node.labe=TRUE, font=1, cex=0.7,tip.color=COLS)
#dev.off()
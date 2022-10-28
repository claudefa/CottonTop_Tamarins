# EEMS, F3 and maps
library(reemsplots2)
library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)
library(RColorBrewer)
library(ggplot2)
library("tidyr")
library(corrplot)

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")

# EEMS CTT

## Part 1: Install rEEMSplots
## Check that the current directory contains the rEEMSplots source directory
#if (file.exists("/Applications/eems/plotting/rEEMSplots/")) {
#  install.packages("/Applications/eems/plotting/rEEMSplots", repos = NULL, type = "source")
#} else {
#  stop("Move to the directory that contains the rEEMSplots source to install the package.")
#}

## Part 2: Generate graphics -------

# all sites
#for (i in 1:10){
  i=1
  print(i)
  mcmcpath = paste0("./Files/CTT-EEMS-nDemes1000-chain",i,"/")
  plotpath = paste0("./Plots/NewIter_", i,"/")
  dir.create(plotpath)
  
  projection_none <- "+proj=longlat +datum=WGS84" 
  projection_mercator <- "+proj=merc +datum=WGS84"
  
  map_world <- getMap()
  map_africa <- map_world[which(map_world@data$continent == "South America"), ]
  map_africa <- spTransform(map_africa, CRSobj = CRS(projection_mercator))
  plots <- make_eems_plots(mcmcpath, longlat = TRUE)
  
  plots <- make_eems_plots(mcmcpath, longlat = TRUE, dpi = 250, add_grid = FALSE,
                  col_grid = "#BBBBBB", add_demes = TRUE, col_demes = "#000000",
                  add_outline = FALSE, col_outline = "#FFFFFF", eems_colors = NULL,
                  prob_level = 0.9, m_colscale = NULL, q_colscale = NULL,
                  add_abline = FALSE)
  
  plots$mrates02
  
  eems.plots(mcmcpath, plotpath, longlat = TRUE, 
             add.grid=FALSE, add.outline = TRUE, add.demes = TRUE,
             add.map = TRUE,projection.in = projection_none,
             projection.out = projection_mercator,
             m.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             q.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             eems.colors = colorRampPalette(brewer.pal(11, "RdBu"))(200),
             out.png=FALSE, min.cex.demes = 0.5, max.cex.demes = 1.5 )
  
  
  
  
  # no maf -----
  i=1
  mcmcpath = paste0("./Files/CTT-EEMS-nDemes1000-chain",i,"_nomaf/")
  plotpath = paste0("./Plots/Iter_", i,"_nomaf/")
  dir.create(plotpath)
  projection_none <- "+proj=longlat +datum=WGS84" 
  projection_mercator <- "+proj=merc +datum=WGS84"
  
  map_world <- getMap()
  map_africa <- map_world[which(map_world@data$continent == "Africa"), ]
  
  map_africa <- spTransform(map_africa, CRSobj = CRS(projection_mercator))
  #Plot all subspecies
  
  eems.plots(mcmcpath, plotpath, longlat = TRUE, 
             add.grid=FALSE, add.outline = TRUE, add.demes = TRUE,
             add.map = TRUE,projection.in = projection_none,
             projection.out = projection_mercator,
             m.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             q.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             eems.colors = colorRampPalette(brewer.pal(11, "RdBu"))(200),
             out.png=FALSE, min.cex.demes = 0.5, max.cex.demes = 1.5 )
  
  
  
  # chr20-------
  i=3
  mcmcpath = paste0("./Files/Baboons-EEMS-nDemes2000-chain",i,"/")
  plotpath = paste0("./Plots/Iter_", 4,"/")
  dir.create(plotpath)
  
  projection_none <- "+proj=longlat +datum=WGS84" 
  projection_mercator <- "+proj=merc +datum=WGS84"
  
  map_world <- getMap()
  map_africa <- map_world[which(map_world@data$continent == "Africa"), ]
  
  map_africa <- spTransform(map_africa, CRSobj = CRS(projection_mercator))
  #Plot
  
  eems.plots(mcmcpath, plotpath, longlat = TRUE, 
             add.grid=FALSE, add.outline = TRUE, add.demes = TRUE,
             add.map = TRUE,projection.in = projection_none,
             projection.out = projection_mercator,
             m.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             q.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             min.cex.demes = 0.5, max.cex.demes = 1.5 ,
             eems.colors = colorRampPalette(brewer.pal(11, "RdBu"))(200), out.png=TRUE)
  
# F3--------------
#F3
#onlyCTT  
f3dat_sites <- read.table("Files/all_sites_onlyCTT.qp3Pop.out",
                            col.names=c("PopA", "PopB", "Outgroup", "F3", "StdErr", "Z", "SNPs"))
  
f3dat_sites$low <- f3dat_sites$F3 - f3dat_sites$StdErr
f3dat_sites$high <- f3dat_sites$F3 + f3dat_sites$StdErr

#f3dat_sites$PopA <- factor(f3dat_sites$PopA, levels=order, ordered=TRUE)
#f3dat_sites$PopB <- factor(f3dat_sites$PopB, levels=order, ordered=TRUE)

f3dat_sites_final <- cbind.data.frame(f3dat_sites$PopA,f3dat_sites$PopB, f3dat_sites$F3)
colnames(f3dat_sites_final) <- c("PopA","PopB","F3")

f3<- ggplot(f3dat_sites_final, aes(PopA, PopB, fill=F3))
f3 <- f3 + geom_tile() +   scale_fill_viridis_c()+
  #scale_fill_gradientn(colours = colorRampPalette(c(brewer.pal( 11, "RdYlBu")))(200)) +   
  ylab("") + xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1))
f3


pdf("Plots/F3_ggplot2.pdf", height = 8, width = 10)
f3
dev.off()

# spread
colnames(f3dat_sites_final) <- c("PopA","PopB","F3")
x <- spread(f3dat_sites_final, PopA, F3, fill=NA,convert = FALSE)
names <- x$PopB
x <- x[,2:ncol(x)]
rownames(x) <- names
y <- data.matrix(x)
y

pdf("Plots/F3_sites_onlyCTT_newVCF.pdf", height=8, width = 10 )
corrplot(y,addrect = 4, na.label = "NA",
         tl.col = "black", type="full",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90, 
         rect.col="black",
         pch.cex = 1.0, method="color",is.corr = FALSE, col=colorRampPalette(c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(200),
         pch.col = "purple")

dev.off()


#only historical 
historical <- c("Arjona_Historical","Caracas_Historical","Cauca_Historical",
                "Coloso_Historical","Mutata_Historical",
                "PlanetaRica_Historical","SanJuan_Historical",
                "Tierralta_Historical","Turbo_Historical","Unknown_Historical")
f3dat_sites_final_hist <- f3dat_sites_final[which(f3dat_sites_final$PopA%in%historical & f3dat_sites_final$PopB%in% historical),]

colnames(f3dat_sites_final_hist) <- c("PopA","PopB","F3")
x <- spread(f3dat_sites_final_hist, PopA, F3, fill=NA,convert = FALSE)
names <- x$PopB
x <- x[,2:ncol(x)]
rownames(x) <- names
y <- data.matrix(x)
y

pdf("Plots/F3_sites_onlyCTT_newVCF_historical.pdf", height=8, width = 10 )
corrplot(y,addrect = 3, na.label = "NA",
         tl.col = "black", type="full",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90, 
         rect.col="black",
         pch.cex = 1.0, method="color",is.corr = FALSE, 
         col=colorRampPalette(c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(200),
         pch.col = "purple")

dev.off()


# with S. geoffroyi
f3dat_sites <- read.table("Files/all_sites_newvcf.qp3Pop.out",
                          col.names=c("PopA", "PopB", "Outgroup", "F3", "StdErr", "Z", "SNPs"))

f3dat_sites$low <- f3dat_sites$F3 - f3dat_sites$StdErr
f3dat_sites$high <- f3dat_sites$F3 + f3dat_sites$StdErr

f3dat_sites_final <- cbind.data.frame(f3dat_sites$PopA,f3dat_sites$PopB, f3dat_sites$F3)
colnames(f3dat_sites_final) <- c("PopA","PopB","F3")

f3<- ggplot(f3dat_sites_final, aes(PopA, PopB, fill=F3))
f3 <- f3 + geom_tile() +   scale_fill_viridis_c()+
  #scale_fill_gradientn(colours = colorRampPalette(c(brewer.pal( 11, "RdYlBu")))(200)) +   
  ylab("") + xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1))
f3

# spread
colnames(f3dat_sites_final) <- c("PopA","PopB","F3")
x <- spread(f3dat_sites_final, PopA, F3, fill=NA,convert = FALSE)
names <- x$PopB
x <- x[,2:ncol(x)]
rownames(x) <- names
y <- data.matrix(x)
y

pdf("Plots/F3_sites_allSamples_newVCF.pdf", height=8, width = 10 )
corrplot(y,addrect = 4, na.label = "NA",
         tl.col = "black", type="full",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90, 
         rect.col="black",
         pch.cex = 1.0, method="color",is.corr = FALSE, col=colorRampPalette(c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(200),
         pch.col = "purple")

dev.off()

# EEMS CTT

## Part 1: Install rEEMSplots
## Check that the current directory contains the rEEMSplots source directory
#if (file.exists("~/Applications/eems/plotting/rEEMSplots/")) {
#  install.packages("~/Applications/eems/plotting/rEEMSplots", repos = NULL, type = "source")
#} else {
#  stop("Move to the directory that contains the rEEMSplots source to install the package.")
#}

## Part 2: Generate graphics
library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)
library(RColorBrewer)
library(ggplot2)
library("tidyr")
library(corrplot)



#library(ghostscript)
setwd("~/Documents/OneDrive - upf.edu/Doctorat/Cotton-top//")

# all sites
#for (i in 1:10){
  i=1
  print(i)
  mcmcpath = paste0("./Files/CTT-EEMS-nDemes1000-chain",i,"/")
  plotpath = paste0("./Plots/Iter_", i,"/")
  dir.create(plotpath)
  
  projection_none <- "+proj=longlat +datum=WGS84" 
  projection_mercator <- "+proj=merc +datum=WGS84"
  
  map_world <- getMap()
  map_africa <- map_world[which(map_world@data$continent == "South America"), ]
  map_africa <- map_world[which(map_world@data$continent == "South America"), ]
  
  map_africa <- spTransform(map_africa, CRSobj = CRS(projection_mercator))
  
  eems.plots(mcmcpath, plotpath, longlat = TRUE, 
             add.grid=FALSE, add.outline = TRUE, add.demes = TRUE,
             add.map = TRUE,projection.in = projection_none,
             projection.out = projection_mercator,
             m.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             q.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             eems.colors = colorRampPalette(brewer.pal(11, "RdBu"))(200),
             out.png=TRUE, min.cex.demes = 0.5, max.cex.demes = 1.5 )
  
  
  
  
  # prunned -----
  i=3
  mcmcpath = paste0("./Files/Baboons-EEMS-nDemes2000-chain",i,"/")
  plotpath = paste0("./Plots/Iter_", i,"/")
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
             out.png=TRUE, min.cex.demes = 0.5, max.cex.demes = 1.5 )
  
  
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
  
  # chrX -------
  i=5
  mcmcpath = paste0("./Files/Baboons-EEMS-nDemes2000-chain",i,"/")
  plotpath = paste0("./Plots/Iter_", i,"/")
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
             min.cex.demes = 0.5, max.cex.demes = 1.5 ,
             eems.colors = colorRampPalette(brewer.pal(11, "RdBu"))(200), out.png=TRUE)
  
  


# F3--------------
#F3
f3dat_sites <- read.table("Files/all_sites.qp3Pop.out",
                          col.names=c("PopA", "PopB", "Outgroup", "F3", "StdErr", "Z", "SNPs"))
f3dat_sites$low <- f3dat_sites$F3 - f3dat_sites$StdErr
f3dat_sites$high <- f3dat_sites$F3 + f3dat_sites$StdErr

order <- c("Papio_cynocephalus_IssaValley,Tanzania","Papio_hamadryas_Filoha,Ethiopia", "Papio_papio_Niokolo-Koba,Senegal",
"Papio_anubis_Tarangire,Tanzania","Papio_anubis_Arusha,Tanzania","Papio_anubis_LakeManyara,Tanzania",
"Papio_anubis_Ngorongoro,Tanzania","Papio_anubis_GogWoreda,Ethiopia","Papio_anubis_Gombe,Tanzania","Papio_anubis_Serengeti,Tanzania",
"Papio_cynocephalus_Ruaha,Tanzania","Papio_cynocephalus_Udzungwa,Tanzania","Papio_cynocephalus_Mikumi,Tanzania",
"Papio_cynocephalus_Selous,Tanzania","Papio_ursinus_(grayfoot)DendroPark,Zambia","Papio_kindae_Chunga,Zambia",
"Papio_cynocephalus_Katavi,Tanzania","Papio_cynocephalus_Mahale,Tanzania")

          

f3dat_sites$PopA <- factor(f3dat_sites$PopA, levels=order, ordered=TRUE)
f3dat_sites$PopB <- factor(f3dat_sites$PopB, levels=order, ordered=TRUE)

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
hola <- as.character(names)
hola1 <- unlist(lapply(1:length(hola), function(i) strsplit(hola[i],"_")[[1]][3]))
hola2 <- unlist(lapply(1:length(hola1), function(i) strsplit(hola1[i],",")[[1]][1]))



pdf("Plots/F3_sites.pdf", height=8, width = 10 )
corrplot(y,addrect = 4, na.label = "NA",
         tl.col = "black", type="full",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90, 
         rect.col="black",
         pch.cex = 1.0, method="color",is.corr = FALSE, col=colorRampPalette(c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(200),
         pch.col = "purple")

dev.off()



# Maps -------------
library('maps')
library('mapdata')
library('maptools')
library("ggplot2")
library("readxl")
library("rgdal") 
library(svglite)
library("tidyr")
library("ggConvexHull")
library("scales")
library("ggrepel")
library("RColorBrewer")
library("viridis")
library("ggpubr")
library("cowplot")
# Construct Map
map <-readOGR("~/Documents/OneDrive - upf.edu/Doctorat/PanAf/TechPaper/Files/Maps/ne_50m_admin_0_countries/","ne_50m_admin_0_countries")
rivers <-readOGR("~/Documents/OneDrive - upf.edu/Doctorat/PanAf/TechPaper/Files/Maps/ne_110m_rivers_lake_centerlines/ne_110m_rivers_lake_centerlines.shp")
lakes <-readOGR("~/Documents/OneDrive - upf.edu/Doctorat/PanAf/TechPaper/Files/Maps/ne_110m_lakes/ne_110m_lakes.shp")
ocean<-readOGR("~/Documents/OneDrive - upf.edu/Doctorat/PanAf/TechPaper/Files/Maps/ne_50m_ocean/ne_50m_ocean.shp")

# CTT
area.1 <- readOGR("Files/redlist_species_data_94c6c081-f927-416b-8ea9-ef72a2300dc4/")
area.1.points <- fortify(area.1)
area.all <- rbind.data.frame(area.1.points)

# Plain Map #c2d9ff #accae5
k <-ggplot() + geom_polygon(data=ocean,aes(x = long,y = lat,group = group),fill = "#a3c2f7",size=0.3) + 
  geom_path(data=ocean,aes(x = long,y = lat,group = group),colour = "#a3c2f7",size=0.7)
k <- k +  geom_polygon(data = map, aes(x = long,y = lat, group = group),fill = "#fdfdfb",size=0.3) + 
  geom_path(data = map,aes(x = long,y = lat,group = group),colour = "grey50",size=0.4) + 
  labs(x = "Longitude", y = "Latitude") + coord_fixed(xlim = c(-79, -74), ylim = c(6, 12), ratio=1) + 
  theme(panel.background = element_rect(fill = '#e3f0ff'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

k
# Subspecies Map

area.all$species <- area.all$group



colors_habitat=rep(c("#e31a1c"),  times=c(15))

#'kindae', 'west yellow', 'eash yellow', 'papio', 'north anubis', 'south anubis', 'hamadryas', 'chacma'
#'#1f78b4', '#FFD87E',      '#FFB100',  '#e31a1c', '#b2df8a',      '#33a02c',      "#8d9da3", '#b15928'

p <- k
p <- p + geom_polygon(data = area.all,aes(x = long,y = lat, fill = group, color=species), color=NA, size=0.01, alpha = 0.6)  +
  ylab("") +xlab("")+
  geom_point(data = area.all,aes(x = long,y = lat, color=species), fill=NA, size=0, alpha = 0)  +
  scale_fill_manual(values=colors_habitat) +   
  scale_color_manual(values=colors_habitat,  guide = guide_legend(override.aes = list(size = 4, alpha=0.6) ) ) +
  geom_path(data=rivers,aes(x = long,y = lat,group = group),colour = "#badaff",size=0.4) + 
  geom_polygon(data=lakes,aes(x = long,y = lat,group = group),fill = "#badaff",alpha=0.4) + 
  guides(fill="none" )+
  #geom_polygon(data = area.ref.points,aes(x = long,y = lat, color=group), linetype="dashed", size=0.8,fill=NA)+
  # scale_color_manual(values=rep(c("#4d4d4d"), times=c(55))) +
  theme(panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.8), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.box.background = element_blank(),
        legend.background = element_blank(),legend.title =element_blank(),
        legend.key=element_rect(fill="white"), legend.text = element_text(size=15,face = "italic"),
        legend.box = "horizontal") 
p

Site_coordinates <- read.csv2("metadata_sites.csv", header = FALSE)

LABELS <- Site_coordinates$V1
LON <- Site_coordinates$V2
LAT <- Site_coordinates$V3
coord.SITES <- cbind.data.frame(LABELS, LON, LAT)
coord.SITES$LABELS <- as.factor(coord.SITES$LABELS)
coord.SITES$LAT <- as.numeric(coord.SITES$LAT)
coord.SITES$LON <- as.numeric(coord.SITES$LON)
source("~/Documents/OneDrive - upf.edu/Doctorat/PanAf/PhaseII/Genotype Analysis/Admixture/scalebar_function.R")

s <- p
s <- s +scaleBar(lon = -75, lat = 6, distanceLon = 50,
                 distanceLat = 25, distanceLegend = 52, dist.unit = "km", arrow.length = 52, arrow.distance = 52) + 
  geom_text_repel(data=coord.SITES,aes(LON, LAT, label=LABELS),color="gray2", size=6, max.overlaps = Inf)+
  geom_point(data=coord.SITES,aes(LON, LAT), alpha=0.8, color="grey2",size=2) + ylab("") +xlab("")+
  geom_rect(aes(xmin = 27, xmax = 41, ymin = -10, ymax = 0),
           fill = "transparent", color = "black", size = 1)+
  theme(panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.8), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.background = element_blank(), legend.position = "none",
        legend.key=element_rect(fill="white"),
        strip.background = element_rect(colour="black", fill="#c7dbf9"),
         legend.box = "horizontal") 
s

pdf("Plots/PlainMap_CTT.pdf", height = 11, width = 11)
s
dev.off()


## Zoom in in the plot 

r <- p
r <- r + 
  geom_text_repel(data=coord.SITES[1:14,][-2,],aes(LON, LAT, label=LABELS),color="gray2", size=6,
                  max.overlaps = Inf)+theme_minimal_hgrid(color = NA)+
 # scaleBar(lon = 38, lat = 0, distanceLon = 100,
    #       distanceLat = 25, distanceLegend = 50, dist.unit = "km", arrow.length = 100, arrow.distance = 50) + 
  geom_point(data=coord.SITES[1:14,],aes(LON, LAT), alpha=0.8, color="grey2",size=2) + ylab("") +xlab("")+
  theme(panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=2), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.box.background = element_rect(fill="white",size = 0.3),
        legend.background = element_blank(), legend.position = "none",
        legend.key=element_rect(fill="white"), 
        strip.background = element_rect(colour="black"),
        legend.box = "horizontal") + 
   coord_fixed(xlim = c(27, 41), ylim = c(-10, 0), ratio=1) 
r

pdf("Plots/PlainMap_baboons_zoomin.pdf", height = 6.5, width = 6.5) 
r
dev.off()

gg_inset_map2 = ggdraw() +
  draw_plot(s) +
  draw_plot(r, x = 0.029, y = 0.26, width = 0.45, height = 0.45)

pdf("Plots/PlainMap_baboons_comb.pdf", height = 11, width = 12) 
gg_inset_map2
dev.off()




#all samples plotted

Site_coordinates <- read.csv2("metadata.csv", header = FALSE)

LABELS <- Site_coordinates$V9
SAMPLES <- Site_coordinates$V1
LON <- Site_coordinates$V7
LAT <- Site_coordinates$V8
coord.SITES <- cbind.data.frame(LABELS, LON, LAT, SAMPLES)
coord.SITES$LABELS <- as.factor(coord.SITES$LABELS)
coord.SITES$LAT <- as.numeric(coord.SITES$LAT)
coord.SITES$LON <- as.numeric(coord.SITES$LON)

s <- p
s <- s + 
  geom_point(data=coord.SITES,aes(LON, LAT, color=LABELS), alpha=0.8, size=3) + ylab("") +xlab("")+
  theme(panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.8), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.box.background = element_rect(fill="white",size = 0.3),
        legend.background = element_blank(),
        legend.key=element_rect(fill="white"),
        strip.background = element_rect(colour="black", fill="#c7dbf9"),
        legend.position = "right", legend.box = "horizontal") + scale_size(range = c(1, 5))
s


pdf("Plots/PlainMap_baboons_ALLPOS.pdf", height = 12, width = 18)
s
dev.off()

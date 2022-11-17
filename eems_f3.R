# EEMS, F3 and maps
library(reemsplots2)
#library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)
library(RColorBrewer)
library(ggplot2)
library("tidyr")
library(corrplot)

setwd("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/")
metadata  <- read_excel("Shotgun_Metadata.xlsx")

# EEMS CTT ------

# map specs

map <-readOGR("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Files/ne_50m_admin_0_countries/","ne_50m_admin_0_countries")
# all sites
for (i in 1:10){
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

  plots <- make_eems_plots(mcmcpath, longlat = TRUE, dpi = 250, add_grid = FALSE,
                  col_grid = "#BBBBBB", add_demes = TRUE, col_demes = "#000000",
                  add_outline = FALSE, col_outline = "#FFFFFF", eems_colors = NULL,
                  prob_level = 0.9, m_colscale = NULL, q_colscale = NULL,
                  add_abline = FALSE)
  
  eems.plots(mcmcpath, plotpath, longlat = TRUE, 
             add.grid=FALSE, add.outline = TRUE, add.demes = TRUE,
             add.map = TRUE,projection.in = projection_none,
             projection.out = projection_mercator,
             m.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             q.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
             eems.colors = colorRampPalette(brewer.pal(11, "RdBu"))(200),
             out.png=FALSE, min.cex.demes = 0.5, max.cex.demes = 1.5 )
  
 new_plots <-plots$mrates01 +
   geom_path(data = map,aes(x = long,y = lat,group = group),colour = "black",size=0.4) + 
    xlab("Longitude") + ylab("Latitude") + theme_minimal() +
   theme(panel.grid.minor = element_blank(),panel.background = element_blank(),
         panel.border = element_rect(color="black", fill=NA, size=0.8), 
         legend.box.background = element_blank(),
         legend.box = "horizontal")  
 ggsave(paste0(plotpath, "mrates01.pdf"),new_plots, width = 6, height = 4)

  new_plots2 <-plots$mrates02 + 
   geom_path(data = map,aes(x = long,y = lat,group = group),colour = "black",size=0.4) + 
   xlab("Longitude") + ylab("Latitude") + theme_minimal() +
   theme(panel.grid.minor = element_blank(),panel.background = element_blank(),
         panel.border = element_rect(color="black", fill=NA, size=0.8), 
         legend.box.background = element_blank(),
         legend.box = "horizontal")  
 ggsave(paste0(plotpath, "mrates02.pdf"), new_plots2, width = 6, height = 4) 

  ggarrange(new_plots, new_plots2, labels = c("A", "B"), widths = c(1,1.255))
  ggsave("Plots/EEMS_final.pdf", width = 8, height = 4) 
  
 } 
  
  
  

# F3--------------
#F3
#onlyCTT  
f3dat_sites <- read.table("Files/all_sites_onlyCTT_final.qp3Pop.out",
                            col.names=c("PopA", "PopB", "Outgroup", "F3", "StdErr", "Z", "SNPs"))
f3dat_sites$low <- f3dat_sites$F3 - f3dat_sites$StdErr
f3dat_sites$high <- f3dat_sites$F3 + f3dat_sites$StdErr




f3<- ggplot(f3dat_sites_final, aes(PopA, F3))
f3 <- f3 + geom_point() +   
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

pdf("Plots/F3_sites_onlyCTT_newVCF.pdf", height=8, width = 10 )
corrplot(y,addrect = 4, na.label = "NA",
         tl.col = "black", type="full",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90, 
         rect.col="black",
         pch.cex = 1.0, method="color",is.corr = FALSE, col=colorRampPalette(c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(200),
         pch.col = "purple")

dev.off()

# only ones with coordinates
order <- c("Tulenapa_Modern", "Mutata_Historical", "Turbo_Historical",
           "Tierralta_Historical",
           "PlanetaRica_Historical",
           "Caracas_Historical", 
           "Coloso_Historical",
           "SanJuan_Historical",
           "SanJuan_Modern","Ceibal_Modern")
f3dat_sites_onlyGPS <- f3dat_sites[which(f3dat_sites$PopA%in%order& f3dat_sites$PopB%in%order),]

f3dat_sites_onlyGPS$PopA <- factor(f3dat_sites_onlyGPS$PopA, levels=order, ordered=TRUE)
f3dat_sites_onlyGPS$PopB <- factor(f3dat_sites_onlyGPS$PopB, levels=order, ordered=TRUE)

f3dat_sites_final <- cbind.data.frame(f3dat_sites_onlyGPS$PopA,f3dat_sites_onlyGPS$PopB, f3dat_sites_onlyGPS$F3)
colnames(f3dat_sites_final) <- c("PopA","PopB","F3")

f3<- ggplot(f3dat_sites_final, aes(PopA, PopB, fill=F3))
f3 <- f3 + geom_tile() +   scale_fill_viridis_c()+
  #scale_fill_gradientn(colours = colorRampPalette(c(brewer.pal( 11, "RdYlBu")))(200)) +   
  ylab("") + xlab("")+ ggtitle("f3(S.midas;PopA,PopB)")+ theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.line = element_blank())
f3
pdf("Plots/F3_onlyGPS.pdf", height = 8, width = 10)
f3
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



# Map per population only CTT all with coordinates-----
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
map <-readOGR("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Files/ne_50m_admin_0_countries/","ne_50m_admin_0_countries")
rivers <-readOGR("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Files/ne_110m_rivers_lake_centerlines/ne_110m_rivers_lake_centerlines.shp")
ocean<-readOGR("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Files/ne_10m_ocean/ne_10m_ocean.shp")

# CTT
area.1 <- readOGR("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Files/redlist_species_data_94c6c081-f927-416b-8ea9-ef72a2300dc4/")
area.1.points <- fortify(area.1)
area.all <- rbind.data.frame(area.1.points)

# Plain Map #c2d9ff #accae5
k <-ggplot() + geom_polygon(data=ocean,aes(x = long,y = lat,group = group),fill = "#d9d9d9",size=0.3) + 
  geom_path(data=ocean,aes(x = long,y = lat,group = group),colour = "#d9d9d9",size=0.7)
k <- k +  geom_polygon(data = map, aes(x = long,y = lat, group = group),fill = "#fdfdfb",size=0.3) + 
  geom_path(data = map,aes(x = long,y = lat,group = group),colour = "grey50",size=0.4) + 
  labs(x = "Longitude", y = "Latitude") + coord_fixed(xlim = c(-78, -74), ylim = c(7, 11), ratio=1) + 
  theme(panel.background = element_rect(fill = '#e3f0ff'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

k
# Subspecies Map
area.all$species <- area.all$group
colors_habitat=rep(c("#f2cc9b"),  times=c(15))


coordinates <- read_excel("coordinates.xlsx")
coordinates$PopB <- coordinates$Site
f3dat_sites <- merge(f3dat_sites, coordinates[ , c("Latitude","Longitude","Site","PopB")], by="PopB") 


for (i in 1:length(coordinates$Site)){
site <- coordinates$Site[i]
p <- k
p <- p + geom_polygon(data = area.all,aes(x = long,y = lat, fill = group), color=NA, size=0.01, alpha = 0.6)  +
  ylab("") +xlab("")+
  geom_point(data = f3dat_sites[f3dat_sites$PopA==site,],
             aes(x = Longitude,y = Latitude, color=F3), size=3)  +
  geom_point(data = coordinates[coordinates$PopB==site,],
             aes(x = Longitude,y = Latitude), color="grey20", size=3, shape=4)  +
  geom_text_repel(data = coordinates,
                  aes(x = Longitude,y = Latitude, label=PopB), color="grey20",size=3)+
  scale_fill_manual(values=colors_habitat) +   
  scale_color_viridis_c()+
  geom_path(data=rivers,aes(x = long,y = lat,group = group),colour = "#badaff",size=0.4) + 
  #geom_polygon(data=lakes,aes(x = long,y = lat,group = group),fill = "#badaff",alpha=0.4) + 
  guides(fill="none" )+
  theme(panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.8), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key=element_rect(fill="white"),
        legend.box = "horizontal") +
  ggtitle(paste0("f3(S.midas;X;",site,")"))
p

pdf(paste0("Plots/Map_f3_",site,".pdf"), width = 6, height = 5)
print(p)
dev.off()
}

i=6
site <- coordinates$Site[i]
p <- k
p <- p + geom_polygon(data = area.all,aes(x = long,y = lat, fill = group), color=NA, size=0.01, alpha = 0.6)  +
  ylab("") +xlab("")+
  geom_point(data = f3dat_sites[f3dat_sites$PopA==site,],
             aes(x = Longitude,y = Latitude, color=F3), size=3)  +
  geom_point(data = coordinates[coordinates$PopB==site,],
             aes(x = Longitude,y = Latitude), color="grey20", size=3, shape=4)  +
  geom_text_repel(data = coordinates,
                  aes(x = Longitude,y = Latitude, label=PopB), color="grey20",size=3)+
  scale_fill_manual(values=colors_habitat) +   
  scale_color_viridis_c()+
  geom_path(data=rivers,aes(x = long,y = lat,group = group),colour = "#badaff",size=0.4) + 
  #geom_polygon(data=lakes,aes(x = long,y = lat,group = group),fill = "#badaff",alpha=0.4) + 
  guides(fill="none" )+
  theme(panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.8), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key=element_rect(fill="white"),
        legend.box = "horizontal") +
  ggtitle(paste0("f3(S.midas;X;",site,")"))
p

i=8
site <- coordinates$Site[i]
q <- k
q <- q + geom_polygon(data = area.all,aes(x = long,y = lat, fill = group), color=NA, size=0.01, alpha = 0.6)  +
  ylab("") +xlab("")+
  geom_point(data = f3dat_sites[f3dat_sites$PopA==site,],
             aes(x = Longitude,y = Latitude, color=F3), size=3)  +
  geom_point(data = coordinates[coordinates$PopB==site,],
             aes(x = Longitude,y = Latitude), color="grey20", size=3, shape=4)  +
  geom_text_repel(data = coordinates,
                  aes(x = Longitude,y = Latitude, label=PopB), color="grey20",size=3)+
  scale_fill_manual(values=colors_habitat) +   
  scale_color_viridis_c()+
  geom_path(data=rivers,aes(x = long,y = lat,group = group),colour = "#badaff",size=0.4) + 
  #geom_polygon(data=lakes,aes(x = long,y = lat,group = group),fill = "#badaff",alpha=0.4) + 
  guides(fill="none" )+
  theme(panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.8), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key=element_rect(fill="white"),
        legend.box = "horizontal") +
  ggtitle(paste0("f3(S.midas;X;",site,")"))
q

ggarrange(p,q,ncol=2)
ggsave("Plots/F3_maps.pdf", width = 10, height = 6)

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

# F3 per sample
f3dat_indv <- read.table("Files/allIndv.qp3Pop.out",
                          col.names=c("PopA", "PopB", "Outgroup", "F3", "StdErr", "Z", "SNPs"))

f3dat_indv$low <- f3dat_indv$F3 - f3dat_indv$StdErr
f3dat_indv$high <- f3dat_indv$F3 + f3dat_indv$StdErr

metadata$PopB <- metadata$Full_ID
f3dat_indv <- merge(f3dat_indv, metadata[ , c("Site","PopB","Coverage","Type")], by="PopB")                            


# spread
f3dat_sites_final <- cbind.data.frame(f3dat_indv$PopA,f3dat_indv$PopB, f3dat_indv$F3)
colnames(f3dat_sites_final) <- c("PopA","PopB","F3")
colnames(f3dat_sites_final) <- c("PopA","PopB","F3")
x <- spread(f3dat_sites_final, PopA, F3, fill=NA,convert = FALSE)
names <- x$PopB
x <- x[,2:ncol(x)]
rownames(x) <- names
y <- data.matrix(x)

corrplot(y,addrect = 4, na.label = "NA",
         tl.col = "black", type="full",order="hclust",na.label.col="white",
         tl.cex = 0.95,  
         tl.srt = 90, 
         rect.col="black",
         pch.cex = 1.0, method="color",is.corr = FALSE, col=colorRampPalette(c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(200),
         pch.col = "purple")


f3<- ggplot(f3dat_indv[f3dat_indv$PopB=="36_AMNH_32702",], aes(PopA, F3, fill=F3))
f3 <- f3 + geom_point() + geom_errorbar(aes(ymin=low, ymax=high), width=.2,
                                          position=position_dodge(0.05))  +
  scale_fill_viridis_c()+ facet_wrap(.~PopB)+
  #scale_fill_gradientn(colours = colorRampPalette(c(brewer.pal( 11, "RdYlBu")))(200)) +   
  ylab("") + xlab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1))
f3


cor.test(f3dat_indv[f3dat_indv$PopA=="36_AMNH_32702",]$Coverage, f3dat_indv[f3dat_indv$PopA=="36_AMNH_32702",]$F3)

f3<- ggplot(f3dat_indv[f3dat_indv$PopA=="36_AMNH_32702",], aes(Coverage, F3))
f3 <- f3 + geom_point() + geom_smooth(method = "lm")+
  scale_fill_viridis_c()
f3

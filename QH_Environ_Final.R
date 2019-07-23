#*******************************************#
#Pulling Bioclim data and pop coordinates together
#*******************************************#

library(stats); library(dismo); library(maptools); library(rgdal)

#For PCA and correlations

#install_github("vqv/ggbiplot")
library(tibble);library(scales);library(ggbiplot);library(ggplot2)
#install_github("vqv/ggbiplot") # must be installed after devtools is loaded in
library("FactoMineR")
library("factoextra")
library(ggrepel)


#*******************************************#
#Read in data and building matrix
#*******************************************#
setwd("E:/Zumwalde/havardii_environmental/QH_EnvironmentalAnalyses")
#setwd("~/Desktop/QH_Environ")
# http://worldclim.org/version2 is where data was pulled from
files <- list.files("WORLDCLIM", pattern='tif', full.names=TRUE) #Load climate files
bioclim2.5 <- stack(files) #Create a raster stack
#dim(bioclim30s) #Looks at data structure

#Extract data for QH Pops
QHloc_Final <- read.csv("QH_Pops_Final.csv", sep=",", header=T, stringsAsFactors = FALSE)
data(wrld_simpl)
climateFinal <- climate <- extract(bioclim2.5, QHloc_Final[,1:2]) #Extract climate data for locations
#write.csv(climateFinal, file = "climate_Final.csv")

#Renaming columns and rows 
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
colnames(climateFinal) <- bioclim_names
rownames(climateFinal) <-  c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E10", "E13", "E14", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12", "WAUX3") 
pops <-  c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E10", "E13", "E14", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12", "WAUX3") 



#*******************************************#

#PCA

#*******************************************#

climate.region <- c(rep("E",10), rep("W",13))
climate.region_Final <- c(rep("E",10), rep("W",13))

#Getting rid of highly correlated variables. See QH_Environ_Normality_Correlations.R script for correlation codes
climateFinal2 <- as.data.frame(climateFinal, row.names = TRUE)
climateFinal2 <- climateFinal2[,-c(1,5,6,13,16,18)]
climateFinal.pca2 <- prcomp(climateFinal2, center = TRUE,scale. = TRUE)
climateFinal.pca2  #loadings

#Adding column for regions
climateFinal2 <- cbind(climate.region,pops,climateFinal2)
colnames(climateFinal2) <- c("reg", "pops", "BIO2","BIO3","BIO4","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO14","BIO15","BIO17","BIO19")
dim(climateFinal2)


#****With Correlated variables removed Final!****#
#res.pca <- PCA(climateFinal2, graph = FALSE)
#var <- get_pca_var(res.pca)
#write.csv(var$contrib, file = "contributions_Final.csv")
#fviz_pca_var(res.pca, col.var = "black")


#Run this first to get % variation explained on x and y axes for next plot
ggbiplot(climateFinal.pca2, ellipse = FALSE, obs.scale = 1.5, var.scale = 1.5, varname.adjust = 3.2, 
         var.axes = TRUE, labels= NULL, groups=climate.region_Final)+
  geom_text_repel(aes(label = rownames(climateFinal)), nudge_x = 0.1, size=4)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_point(aes(shape=climate.region_Final, color=climate.region_Final, size=3))+
  theme_light(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "bottom") + 
  theme(legend.position = "none")

#Run script Arrow_edits.R to change color of arrows for this graph.

#Final Figure with Axes changed
#if color changes to the arrows should want to be made, run R script Arrow_edits.R before running this plot
e <- ggbiplot(climateFinal.pca2, ellipse = FALSE, obs.scale = 1.25, var.scale = 1.5, varname.adjust = 1.5, 
              var.axes = FALSE, labels= NULL, groups=climate.region_Final)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_point(aes(shape=climate.region_Final, color=climate.region_Final, size=1, labels=FALSE))+
  geom_text_repel(aes(label = rownames(climateFinal)), color = "black", nudge_x = 0.2, nudge_y = 0.2, size=4.5)+
  xlab(expression("Axis 1 (42.4%)")) + 
  ylab(expression("Axis 2 (21.5%)"))+
  theme_classic(base_size = 10) +
  ggtitle("(C) Environmental") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black")) +
  theme(legend.position = "none") +
  scale_x_reverse()
e

e$layers <- c(e$layers, e$layers[[3]])
e$layers <- c(e$layers, e$layers[[1]])
e



# Run script graphing_3panelfig.R to combine environmental PCA with genetic and morphological PCAs


#*******************************************#
#Pulling Bioclim data and pop coordinates together
#*******************************************#

library(stats); library(dismo); library(maptools); library(rgdal)

#For PCA and correlations

#install_github("vqv/ggbiplot")
library(tibble);library(scales);library(devtools);library(ggbiplot);library(ggplot2)
#install_github("vqv/ggbiplot") # must be installed after devtools is loaded in
library("FactoMineR")
library("factoextra")
library(ggrepel)

setwd("E:/Zumwalde/havardii_environmental/data")

#*******************************************#
#Read in data and building matrix
#*******************************************#

# http://worldclim.org/version2 is where data was pulled from
files <- list.files("WORLDCLIM", pattern='tif', full.names=TRUE) #Load climate files
bioclim30s <- stack(files) #Create a raster stack
#dim(bioclim30s) #Looks at data structure

#Extract data for QH Pops
QHloc_E15Removed <- read.csv("QH_Pops_E15Removed.csv", sep=",", header=T, stringsAsFactors = FALSE)
data(wrld_simpl)
climateE15Removed <- climate <- extract(bioclim30s, QHloc_E15Removed[,1:2]) #Extract climate data for locations

#Renaming columns and rows 
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
colnames(climateE15Removed) <- bioclim_names
rownames(climateE15Removed) <-  c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E10", "E13", "E14", "EAUX8", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12", "WAUX3") 


#*******************************************#

# Testing for Correlation

#*******************************************#

library("Hmisc")
corrclimE15Removed <- cor(climateE15Removed)
round(corrclimE15Removed,3)
cor(climate, use = "complete.obs")
corrclimE15Removed2 <- rcorr(as.matrix(corrclimE15Removed))
corrclimE15Removed2
corrclimE15Removed2$r
corrclimE15Removed2$P
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
QHcorrmat_E15Removed <- flattenCorrMatrix(corrclimE15Removed2$r, corrclimE15Removed2$P)
#write.csv(QHcorrmat_E15Removed, file = "corrmat_E15Removed.csv")
symbols_E15Removed <- symnum(corrclimE15Removed, abbr.colnames = FALSE)
#write.csv(symbols_E15Removed, file = "corrmatsymbols_E15Removed.csv")




#*******************************************#

#PCA

#*******************************************#



#****With Correlated variables Removed and E15 Removed****#
climateE15Removed2 <- climateE15Removed[,-c(1,5,6,13,16,18)] #biological 
climateE15Removed.pca2 <- prcomp(climateE15Removed2, center = TRUE,scale. = TRUE)
climateE15Removed.pca2  #loadings
climate.region_E15Removed <- c(rep("E",12), rep("W",13))

ggbiplot(climateE15Removed.pca2, ellipse = FALSE, obs.scale = 1.5, var.scale = 1.5, varname.adjust = 3.25, 
         var.axes = TRUE, labels= NULL, groups=climate.region_E15Removed)+
  geom_text_repel(aes(label = rownames(climateE15Removed)), nudge_x = 0.1, size=4)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_point(aes(shape=climate.region_E15Removed, color=climate.region_E15Removed, size=3))+
  theme_light(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "bottom") + 
  theme(legend.position = "none")

#Final Figure with Axes changed
ggbiplot(climateE15Removed.pca2, ellipse = FALSE, obs.scale = 1.5, var.scale = 1.5, varname.adjust = 3.3, 
         var.axes = TRUE, labels= NULL, groups=climate.region_E15Removed)+
  geom_text_repel(aes(label = rownames(climateE15Removed)), nudge_x = 0.1, size=4)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_point(aes(shape=climate.region_E15Removed, color=climate.region_E15Removed, size=3))+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  xlab(expression("Axis 1 (40.2%)")) + 
  ylab(expression("Axis 2 (24.6%)"))+
  ggtitle("Environmental") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black")) +
  theme(legend.position = "none")

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
setwd("E:/Zumwalde/havardii_environmental/data")
#setwd("~/Desktop/QH_Environ")
# http://worldclim.org/version2 is where data was pulled from
files <- list.files("WORLDCLIM", pattern='tif', full.names=TRUE) #Load climate files
bioclim2.5 <- stack(files) #Create a raster stack
#dim(bioclim30s) #Looks at data structure

#Extract data for QH Pops
setwd("~/Documents/GitHub/QH_EnvironmentalAnalyses")
QHloc_Final <- read.csv("QH_Pops_Final.csv", sep=",", header=T, stringsAsFactors = FALSE)
data(wrld_simpl)
climateFinal <- climate <- extract(bioclim2.5, QHloc_Final[,1:2]) #Extract climate data for locations
#write.csv(climateFinal, file = "climate_Final.csv")


#Renaming columns and rows 
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
colnames(climateFinal) <- bioclim_names
rownames(climateFinal) <-  c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E10", "E13", "E14", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12", "WAUX3") 
pops <-  c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E10", "E13", "E14", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12", "WAUX3") 

library(dplyr)

climate.region <- c(rep("E",10), rep("W",13))
climateFinal2 <- as.data.frame(climateFinal, row.names = TRUE)
colnames(climateFinal2) <- c("reg", "pops", "BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
climateFinal2 <- cbind(climate.region,pops,climateFinal2)
dim(climateFinal2)

library("ggpubr")
ggboxplot(climateFinal2, x = "reg", y = "BIO10", 
          color = "reg",
          order = c("E", "W"),
          ylab = "Region", xlab = "BIO10")

res.aov <- aov(BIO19 ~ reg, data = climateFinal2)
summary(res.aov)


#tapply(climateFinal2$BIO1, climate.region, mean)
#tapply(climateFinal2$BIO1, climate.region, sd)


#*******************************************#

# Testing for Correlation

#*******************************************#

library("Hmisc")
corrclimFinal <- cor(climateFinal)
round(corrclimFinal,3)
cor(climate, use = "complete.obs")
corrclimFinal2 <- rcorr(as.matrix(corrclimFinal))
corrclimFinal2
corrclimFinal2$r
corrclimFinal2$P
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
QHcorrmat_Final <- flattenCorrMatrix(corrclimFinal2$r, corrclimFinal2$P)
#write.csv(QHcorrmat_Final, file = "corrmat_Final.csv")
symbols_Final <- symnum(corrclimFinal, abbr.colnames = FALSE)
#write.csv(symbols_Final, file = "corrmatsymbols_Final.csv")




#*******************************************#

#PCA

#*******************************************#



#****With Correlated variables Final!****#
climateFinal2 <- climateFinal[,-c(1,5,6,13,16,18)] #biological 
climateFinal.pca2 <- prcomp(climateFinal2, center = TRUE,scale. = TRUE)
climateFinal.pca2  #loadings
climate.region_Final <- c(rep("E",10), rep("W",13))

res.pca <- PCA(climateFinal2, graph = FALSE)
var <- get_pca_var(res.pca)
write.csv(var$contrib, file = "contributions_Final.csv")
fviz_pca_var(res.pca, col.var = "black")

ggbiplot(climateFinal.pca2, ellipse = FALSE, obs.scale = 1.5, var.scale = 1.5, varname.adjust = 3.2, 
         var.axes = TRUE, labels= NULL, groups=climate.region_Final)+
  geom_text_repel(aes(label = rownames(climateFinal)), nudge_x = 0.1, size=4)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_point(aes(shape=climate.region_Final, color=climate.region_Final, size=3))+
  theme_light(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "bottom") + 
  theme(legend.position = "none")

#Final Figure with Axes changed

#if color changes to the arrows should want to be made, run R script Arrow_edits.R before running this plot
g <- ggbiplot(climateFinal.pca2, ellipse = FALSE, obs.scale = 1.25, var.scale = 1.5, varname.adjust = 1.5, 
              var.axes = TRUE, labels= NULL, groups=climate.region_Final)+
  geom_text_repel(aes(label = rownames(climateFinal)), color = "black", nudge_x = 0.25, nudge_y = 0.1, size=4.5)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_point(aes(shape=climate.region_Final, color=climate.region_Final, size=3))+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  xlab(expression("Axis 1 (42.4%)")) + 
  ylab(expression("Axis 2 (21.5%)"))+
  ggtitle("Environmental") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black")) +
  theme(legend.position = "none") +
  scale_x_reverse()
g

g$layers <- c(g$layers, g$layers[[3]])
g$layers <- c(g$layers, g$layers[[1]])
g

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
QHloc <- read.csv("QH_Pops.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
#QHloc <- read.csv("BeckHob_QHOccur_Vetted.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
summary(QHloc) #Checks that file loaded properly and looks at columns

data(wrld_simpl)

climate <- extract(bioclim30s, QHloc[,1:2]) #Extract climate data for locations

#dim(climate) #Checking dimensions

#Pulling in elevation data
#library(elevatr)
#prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#df_elev_epqs <- get_elev_point(QHloc, prj = prj_dd, src = "epqs")
#data.frame(df_elev_epqs)

#QHsoils <- read.csv("QH_Soils.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
#QHsoils <- read.csv("Soils_Data_AllPops_Quha_Vetted.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
#QHsoils[QHsoils == 999] <- NA

#climate <- cbind(clim, df_elev_epqs$elevation)

#Renaming columns and rows 
#bioclim_names <- c("1_Ann_T","2_Diurnal_Range","3_Isothermality","4_T_Seasonality","5_MaxT_Wrmst_Month","6_MinT_Cldst_Month","7_T_Ann_Range","8_T_Wettest_Qtr","9_T_Driest_Qtr","10_T_Wrmst_Qtr","11_T_Cldst_Qtr","12_Ann_Precip","13_P_Wettest_Month","14_P_Driest_Month","15_P_Seasonality","16_P_Wettest_Qtr","17_P_Driest_Qtr","18_P_Wrmst_Qtr","19_P_Cldst_Qtr")
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
colnames(climate) <- bioclim_names
rownames(climate) <-  c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E10", "E13", "E14", "E15", "EAUX8", "W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12", "WAUX3") 

##PCA and checking contributions
#climate.pca <- prcomp(climate, center = TRUE,scale. = TRUE)
#climate.region <- c(rep("E",13), rep("W",13))
#climate.region <- c(rep("EB_E",221),rep("SH_E",22), rep("SH_W",17),rep("EB_W",43), rep("Unknown",4))
#climate.pca  #loadings
#res.pca <- PCA(climate, graph = FALSE)
#var <- get_pca_var(res.pca)
#write.csv(var$contrib, file = "contributions.csv")
#fviz_pca_var(res.pca, col.var = "black")



#*******************************************#

# Testing for Correlation

#*******************************************#

library("Hmisc")
corrclim <- cor(climate)
round(corrclim,3)
cor(climate, use = "complete.obs")
corrclim2 <- rcorr(as.matrix(climate))
corrclim2
#The output of the function rcorr() is a list containing the following elements : 
#- r : the correlation matrix 
#- n : the matrix of the number of observations used in analyzing each pair of variables 
#- P : the p-values corresponding to the significance levels of correlations.
# Extract the correlation coefficients
corrclim2$r
# Extract p-values
corrclim2$P
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
QHcorrmat <- flattenCorrMatrix(corrclim2$r, corrclim2$P)
#write.csv(QHcorrmat, file = "corrmat.csv")
symbols <- symnum(corrclim, abbr.colnames = FALSE)
#write.csv(symbols, file = "corrmatsymbols_26Pops.csv")

#*******************************************#

#Removing Correlated Variables Above 0.9

#*******************************************#
library("Hmisc")

climate2 <- climate[,-c(1,5,6,13,16,18)] #biological 



#*******************************************#

#PCA

#*******************************************#


# ***** PCA without removing correlated variables ****#

climate.pca <- prcomp(climate, center = TRUE,scale. = TRUE)
climate.region <- c(rep("E",13), rep("W",13))
#climate.region <- c(rep("EB_E",221),rep("SH_E",22), rep("SH_W",17),rep("EB_W",43), rep("Unknown",4))
climate.pca  #loadings
#res.pca <- PCA(climate, graph = FALSE)
#var <- get_pca_var(res.pca)
#write.csv(var$contrib, file = "contributions.csv")
#fviz_pca_var(res.pca, col.var = "black")


#****With Correlated variables Removed****#
climate2.pca2 <- prcomp(climate2, center = TRUE,scale. = TRUE)
climate.region <- c(rep("E",13), rep("W",13))

ggbiplot(climate2.pca2, ellipse = FALSE, obs.scale = 1.5, var.scale = 1.5, varname.adjust = 3.25, 
    var.axes = TRUE, labels= NULL, groups=climate.region)+
  geom_text_repel(aes(label = rownames(climate2)), nudge_x = 0.1, size=4)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_point(aes(shape=climate.region, color=climate.region, size=3))+
  theme_light(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "bottom") + 
  theme(legend.position = "none")

#Final Figure with Axes changed to PCA1 and PCA2 values from previous graph
ggbiplot(climate2.pca2, ellipse = FALSE, obs.scale = 1.5, var.scale = 1.5, varname.adjust = 3.1, 
  var.axes = TRUE, labels= NULL, groups=climate.region)+
  geom_text_repel(aes(label = rownames(climate2)), nudge_x = 0.1, size=4)+
  scale_colour_manual(values = c("red", "blue"))+
  geom_point(aes(shape=climate.region, color=climate.region, size=3))+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  xlab(expression("Axis 1 (41.4%)")) + 
  ylab(expression("Axis 2 (24.1%)"))+
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black")) +
  theme(legend.position = "none")


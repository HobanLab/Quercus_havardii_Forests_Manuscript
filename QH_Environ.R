#*******************************************#
#Pulling Bioclim data and pop coordinates together
#*******************************************#


library(stats); library(dismo); library(maptools); library(rgdal)

setwd("E:/Zumwalde/havardii_environmental/data")


#*******************************************#

#Read in data and building matrix

#*******************************************#

# http://worldclim.org/version2 is where data was pulled from
files <- list.files("WORLDCLIM", pattern='tif', full.names=TRUE) #Load climate files
bioclim2.5 <- stack(files) #Create a raster stack
#dim(bioclim2.5) #Looks at data structure

#Extract data for QH Pops
QHloc <- read.csv("BeckHob_QHOccur_NoUnknown.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
summary(QHloc) #Checks that file loaded properly and looks at columns

data(wrld_simpl)
#plot(wrld_simpl, add=T) #can look at how samples are distributed
#points(QHloc[,1:2], col="red", cex=0.5)

clim <- extract(bioclim2.5, QHloc[,1:2]) #Extract climate data for locations
#dim(climate) #Checking dimensions

#Pulling in elevation data
library(elevatr)
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
df_elev_epqs <- get_elev_point(QHloc, prj = prj_dd, src = "epqs")
# data.frame(df_elev_epqs)

climate <- cbind(clim, df_elev_epqs$elevation)

View(climate)
#Renaming columns and rows 
#bioclim_names <- c("1_Ann_T","2_Diurnal_Range","3_Isothermality","4_T_Seasonality","5_MaxT_Wrmst_Month","6_MinT_Cldst_Month","7_T_Ann_Range","8_T_Wettest_Qtr","9_T_Driest_Qtr","10_T_Wrmst_Qtr","11_T_Cldst_Qtr","12_Ann_Precip","13_P_Wettest_Month","14_P_Driest_Month","15_P_Seasonality","16_P_Wettest_Qtr","17_P_Driest_Qtr","18_P_Wrmst_Qtr","19_P_Cldst_Qtr")
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19", "Elev")
colnames(climate) <- bioclim_names
rownames(climate) <-  QHloc$Pop



#*******************************************#

# Non-metric Multidimensional Scaling

#*******************************************#

library(vegan)
library(ape)
library(dplyr)

climate %>%
  metaMDS(trace = F) %>%
  ordiplot(type = "none") %>%
  text("sites")

dist <- vegdist(climate,  method = "euc")

NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

NMDS.scree(dist)

#  A plot of stress (a measure of goodness-of-fit) vs. dimensionality can be used to assess the proper choice of dimensions. 
#  The stress values themselves can be used as an indicator. 
#  Stress values >0.2 are generally poor and potentially uninterpretable, 
#  whereas values <0.1 are good and <0.05 are excellent, leaving little danger of misinterpretation. 
#  Stress values between 0.1 and 0.2 are useable but some of the distances will be misleading. 
#  Finding the inflexion point can instruct the selection of a minimum number of dimension





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
#write.csv(symbols, file = "corrmatsymbols.csv")


#*******************************************#

#Removing Correlated Variables Above 0.9

#*******************************************#

climate2 <- climate[,-c(1,5,6,13,14,16,18)]


res <- cor(climate2)
round(res,3)
cor(climate2, use = "complete.obs")

res2 <- rcorr(as.matrix(climate2))
res2
res2$r
res2$P

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

QHcorrmat_varsremoved <- flattenCorrMatrix(res2$r, res2$P)
#write.csv(QHcorrmat_varsremoved, file = "corrmat.csv")

symbols_varsremoved <- symnum(res, abbr.colnames = FALSE)
#write.csv(symbols_varsremoved, file = "corrmatsymbols_varsremoved.csv")


#*******************************************#

#PCA

#*******************************************#

#install.packages("tibble")
#install.packages("scales")
#install_github("vqv/ggbiplot")

library(tibble);library(scales);library(devtools);library(ggbiplot);library(ggplot2)


# ***** PCA without removing correlated variables ****#


climate.pca <- prcomp(climate, center = TRUE,scale. = TRUE)
climate.region <- c(rep("E",243), rep("W",60))
#climate.region <- c(rep("EB_E",221),rep("SH_E",22), rep("SH_W",17),rep("EB_W",43), rep("Unknown",4))

climate.pca  #loadings
#With arrows
ggbiplot(climate.pca,ellipse=TRUE,obs.scale = 1.5, var.scale = 1.5, var.axes = TRUE, groups=climate.region)+
  scale_colour_manual(name="Key", values = c("gray", "black")) +
  ggtitle("PCA of Bioclim Variables")+
  theme_minimal()+
  theme(legend.position = "bottom")


#****With Correlated variables Removed****#

climate.pca2 <- prcomp(climate2, center = TRUE,scale. = TRUE)
climate.pca2  #loadings

ggbiplot(climate.pca2,ellipse=TRUE,obs.scale = 1.5, var.scale = 1.5, var.axes = TRUE, groups=climate.region)+
  scale_colour_manual(name="Key", values = c("gray", "black")) +
  ggtitle("PCA of Bioclim Variables")+
  theme_minimal()+
  theme(legend.position = "bottom")







#*****************Miscellaneous code**************************#

#Without arrows
#ggbiplot(climate.pca,ellipse=TRUE,obs.scale = 1.5, var.scale = 1.5, var.axes = FALSE,  labels=rownames(climate), groups=climate.region)+
#  scale_colour_manual(name="Key", values = c("red", "purple", "dark blue")) +
#  ggtitle("PCA of Bioclim Variables")+
#  theme_minimal()+
#  theme(legend.position = "bottom")

ggbiplot(climate.pca,ellipse=TRUE,obs.scale = 1.5, var.scale = 1.5, var.axes = FALSE, labels=rownames(climate),groups=climate.region)+
  scale_colour_manual(name="Key", values = c("gray", "lime green", "black", "purple", "orange")) +
  ggtitle("PCA of Bioclim Variables")+
  theme_minimal()+
  theme(legend.position = "bottom")

#Axes 1 and 3
#ggbiplot(climate.pca,ellipse=TRUE,choices=c(1,3),   labels=rownames(climate), groups=climate.region)


#-------------------------------------------#
#-------------------------------------------#

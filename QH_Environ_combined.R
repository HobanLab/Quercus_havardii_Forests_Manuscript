#-------------------------------------------#
#Pulling Bioclim data and pop coordinates together


#install.packages("dismo")
#install.packages("maptools")
library(stats); library(dismo); library(maptools); library(rgdal)

setwd("E:/Zumwalde/havardii_environmental/data")

# http://worldclim.org/version2
files <- list.files("WORLDCLIM", pattern='tif', full.names=TRUE) #Load climate files
bioclim2.5 <- stack(files) #Create a raster stack
#dim(bioclim2.5) #Looks at data structure
#bioclim2.5 #Looks at data structure
#plot(bioclim2.5[[1:9]]) #Looks at data structure

#Extract data for QH Pops
QHloc <- read.csv("BeckHob_QHOccur.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
summary(QHloc) #Checks that file loaded properly and looks at columns

data(wrld_simpl)
#plot(wrld_simpl, add=T) #can look at how samples are distributed
#points(QHloc[,1:2], col="red", cex=0.5)

climate <- extract(bioclim2.5, QHloc[,1:2]) #Extract climate data for locations
#dim(climate) #Checking dimensions

#Renaming columns and rows 
#bioclim_names <- c("1_Ann_T","2_Diurnal_Range","3_Isothermality","4_T_Seasonality","5_MaxT_Wrmst_Month","6_MinT_Cldst_Month","7_T_Ann_Range","8_T_Wettest_Qtr","9_T_Driest_Qtr","10_T_Wrmst_Qtr","11_T_Cldst_Qtr","12_Ann_Precip","13_P_Wettest_Month","14_P_Driest_Month","15_P_Seasonality","16_P_Wettest_Qtr","17_P_Driest_Qtr","18_P_Wrmst_Qtr","19_P_Cldst_Qtr")
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
colnames(climate) <- bioclim_names
rownames(climate) <-  QHloc$Pop


#-------------------------------------------#
#Correlation


res <- cor(climate)
round(res,3)
cor(climate, use = "complete.obs")

library("Hmisc")
res2 <- rcorr(as.matrix(climate))
res2

#The output of the function rcorr() is a list containing the following elements : 
#- r : the correlation matrix 
#- n : the matrix of the number of observations used in analyzing each pair of variables 
#- P : the p-values corresponding to the significance levels of correlations.

# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P


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

QHcorrmat <- flattenCorrMatrix(res2$r, res2$P)
write.csv(QHcorrmat, file = "corrmat.csv")

symbols <- symnum(res, abbr.colnames = FALSE)
write.csv(symbols, file = "corrmatsymbols.csv")





#-------------------------------------------#
#PCA

#install.packages("tibble")
#install.packages("scales")
#install_github("vqv/ggbiplot")

library(tibble);library(scales);library(devtools);library(ggbiplot);library(ggplot2)

climate.pca <- prcomp(climate, center = TRUE,scale. = TRUE)
climate.region <- c(rep("EB_E",220),rep("SH_E",22), rep("SH_W",17),rep("EB_W",42), rep("Unknown",4))

summary(climate.pca)

climate.pca  #loadings

#With arrows
ggbiplot(climate.pca,ellipse=TRUE,obs.scale = 1.5, var.scale = 1.5, var.axes = TRUE, groups=climate.region)+
  scale_colour_manual(name="Key", values = c("red", "purple", "dark blue")) +
  ggtitle("PCA of Bioclim Variables")+
  theme_minimal()+
  theme(legend.position = "bottom")


#Without arrows
ggbiplot(climate.pca,ellipse=TRUE,obs.scale = 1.5, var.scale = 1.5, var.axes = FALSE,  labels=rownames(climate), groups=climate.region)+
  scale_colour_manual(name="Key", values = c("red", "purple", "dark blue")) +
  ggtitle("PCA of Bioclim Variables")+
  theme_minimal()+
  theme(legend.position = "bottom")

ggbiplot(climate.pca,ellipse=TRUE,obs.scale = 1.5, var.scale = 1.5, var.axes = FALSE, labels=rownames(climate),groups=climate.region)+
  scale_colour_manual(name="Key", values = c("gray", "lime green", "black", "purple", "orange")) +
  ggtitle("PCA of Bioclim Variables")+
  theme_minimal()+
  theme(legend.position = "bottom")



#Axes 1 and 3
#ggbiplot(climate.pca,ellipse=TRUE,choices=c(1,3),   labels=rownames(climate), groups=climate.region)


#-------------------------------------------#
#-------------------------------------------#

#Determining and removing correlated variables

#install.packages("caret")
#library(caret)
#cor_mat = cor(climate)
#highlyCorrelated = findCorrelation(cor_mat, cutoff = 0.1)
#names(climate)[highlyCorrelated]
#NULL = no variables are highly correlated

#Testing another version
#cors <- cor(as.matrix(climate), method="pearson", use="everything")

#highlyCorrelated2 = findCorrelation(cors, cutoff = 0.1)
#names(climate)[highlyCorrelated2]

#install.packages("corrplot")
#library(corrplot)
#cor.climate <-cor(climate)
#colnames(cor.climate) <- bioclim_names
#head(round(cor.climate,5))
#corrplot(cor.climate, method = "color", type = "lower")
#Reorder
#corrplot(cor.climate, method = "color", type = "lower", order = "hclust")
#cor.climate
#Computing p-value of correlations
#cor.mtest <- function(mat, ...) {
mat <- as.matrix(mat)
n <- ncol(mat)
p.mat<- matrix(NA, n, n)
diag(p.mat) <- 0
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    tmp <- cor.test(mat[, i], mat[, j], ...)
    p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
  }
}
colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
p.mat
}
#p.mat <- cor.mtest(M)
#corrplot(M, type="upper", order="hclust", 
#         p.mat = p.mat, sig.level = 0.1, insig = "blank")

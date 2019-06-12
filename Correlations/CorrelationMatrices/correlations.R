#*******************************************#
Checking Correlations
#*******************************************#

##PCA and checking contributions
climate.pca <- prcomp(climate, center = TRUE,scale. = TRUE)
climate.region <- c(rep("E",13), rep("W",13))
#climate.region <- c(rep("EB_E",221),rep("SH_E",22), rep("SH_W",17),rep("EB_W",43), rep("Unknown",4))
climate.pca  #loadings
res.pca <- PCA(climate, graph = FALSE)
var <- get_pca_var(res.pca)
#write.csv(var$contrib, file = "contributions.csv")
fviz_pca_var(res.pca, col.var = "black")




#*******************************************#

# Testing for Correlation

#*******************************************#
library("Hmisc")

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





##EAST ONLY
QHlocEast <- read.csv("QH_Pops_East.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
data(wrld_simpl)
clim_East <- extract(bioclim2.5, QHlocEast[,1:2]) #Extract climate data for locations
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
df_elev_epqs_East <- get_elev_point(QHlocEast, prj = prj_dd, src = "epqs")
data.frame(df_elev_epqs_East)
climate_East <- cbind(clim_East, df_elev_epqs_East$elevation)
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19", "Elev")
colnames(climate_East) <- bioclim_names
rownames(climate_East) <-  QHlocWEast$Pop


climate.pca.east <- prcomp(climate_East, center = TRUE,scale. = TRUE)
climate.region.east <- c(rep("E",13))
climate.pca.east  #loadings
res.pca.east <- PCA(climate_East, graph = FALSE)
var.east <- get_pca_var(res.pca.east)
write.csv(var.east$contrib, file = "contributionsEast.csv")
fviz_pca_var(res.pca.east, col.var = "black")

corrclim.east <- cor(climate_East)
round(corrclim.east,3)
cor(climate_East, use = "complete.obs")
corrclim.east2 <- rcorr(as.matrix(climate_East))
corrclim.east2
corrclim.east2$r
corrclim.east2$P
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
QHcorrmat.east <- flattenCorrMatrix(corrclim.east2$r, corrclim.east2$P)
#write.csv(QHcorrmat, file = "corrmat_east.csv")
symbols.east <- symnum(corrclim.east, abbr.colnames = FALSE)
write.csv(symbols.east, file = "corrmatsymbols_east.csv")





##WEST ONLY
QHlocWest <- read.csv("QH_Pops_West.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
clim_West <- extract(bioclim2.5, QHlocWest[,1:2]) #Extract climate data for locations
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
df_elev_epqs_West <- get_elev_point(QHlocWest, prj = prj_dd, src = "epqs")
data.frame(df_elev_epqs_West)
climate_West <- cbind(clim_West, df_elev_epqs_East$elevation)
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19", "Elev")
colnames(climate_West) <- bioclim_names
rownames(climate_West) <-  QHlocWest$Pop

climate.pca.west <- prcomp(climate_West, center = TRUE,scale. = TRUE)
climate.region.west <- c(rep("E",13))
climate.pca.west  #loadings
res.pca.west <- PCA(climate_West, graph = FALSE)
var.west <- get_pca_var(res.pca.west)
write.csv(var.west$contrib, file = "contributionsWest.csv")
fviz_pca_var(res.pca.west, col.var = "black")

corrclim.west <- cor(climate_West)
round(corrclim.west,3)
cor(climate_West, use = "complete.obs")
corrclim.west2 <- rcorr(as.matrix(climate_West))
corrclim.west2
corrclim.west2$r
corrclim.west2$P
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
QHcorrmat.west <- flattenCorrMatrix(corrclim.west2$r, corrclim.west2$P)
#write.csv(QHcorrmat, file = "corrmat_west.csv")
symbols.west <- symnum(corrclim.west, abbr.colnames = FALSE)
write.csv(symbols.west, file = "corrmatsymbols_west.csv")


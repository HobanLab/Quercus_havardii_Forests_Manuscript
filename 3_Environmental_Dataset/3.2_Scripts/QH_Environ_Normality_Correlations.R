

#Checking for Normality and MANOVA


library(mvnormtest); library("ggpubr"); library(dplyr)


# From the output, the p-value > 0.05 implying that the distribution of the data 
# are not significantly different from normal distribution. 
# In other words, we can assume the normality.

#Bioclim variables not used in PCA
#shapiro.test(climateFinal2$BIO1)
#shapiro.test(climateFinal2$BIO5) #not normal p = 0.02368
#shapiro.test(climateFinal2$BIO6)
#shapiro.test(climateFinal2$BIO13) #not normal p = 0.0121
#shapiro.test(climateFinal2$BIO16) #not normal p = 0.001442
#shapiro.test(climateFinal2$BIO18) #not normal p = 0.0002872


#Bioclim variables used in PCA

shapiro.test(climateFinal2$BIO2) #p = 0.5012
shapiro.test(climateFinal2$BIO3) #p = 0.1046
shapiro.test(climateFinal2$BIO4) #p = 0.0604
shapiro.test(climateFinal2$BIO7) #p = 0.6214
shapiro.test(climateFinal2$BIO11)#p = 0.1942

shapiro.test(climateFinal2$BIO8) #not normal p = 0.009768
hist(climateFinal2$BIO8)
BIO8_Trans <- sqrt(climateFinal2$BIO8)
hist(BIO8_Trans)
shapiro.test(BIO8_Trans) #can't get data to normalize 

shapiro.test(climateFinal2$BIO9) #not normal p = 0.0153
hist(climateFinal2$BIO9)
BIO9_Trans <- sqrt(climateFinal2$BIO9)
hist(BIO9_Trans)
shapiro.test(BIO9_Trans) # p normal = 0.0886


shapiro.test(climateFinal2$BIO10) #not normal p = 0.01895
hist(climateFinal2$BIO10)
BIO10_Trans <- (climateFinal2$BIO10)^3
hist(BIO10_Trans)
shapiro.test(BIO10_Trans) # p normal = 0.1904


shapiro.test(climateFinal2$BIO12) #not normal p = 0.008731
hist(climateFinal2$BIO12)


shapiro.test(climateFinal2$BIO14) #not normal p = 0.004716
hist(climateFinal2$BIO14)

shapiro.test(climateFinal2$BIO15) #not normal p = 0.009577
hist(climateFinal2$BIO15)


shapiro.test(climateFinal2$BIO17) #not normal p = 0.0008614
hist(climateFinal2$BIO17)


shapiro.test(climateFinal2$BIO19) #not normal p = 0.001092
hist(climateFinal2$BIO19)



ggboxplot(climateFinal2, x = "climate.region", y = "BIO1", 
          color = "climate.region",
          order = c("E", "W"),
          ylab = "Region", xlab = "BIO10")

res.aov <- aov(BIO19 ~ climate.region, data = climateFinal2)
summary(res.aov)



res.man <- manova(climateFinal2)


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



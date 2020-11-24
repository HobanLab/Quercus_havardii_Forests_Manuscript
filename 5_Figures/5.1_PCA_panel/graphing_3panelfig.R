library(ggbiplot)
library(ggrepel)

#-----------------#
#Environmental PCA
#-----------------#

#Run script QH_Environ_Final.R first to get Environmental PCA


#-----------------#
#Morphological PCA
#-----------------#

setwd("E:/Zumwalde/havardii_environmental/QH_EnvironmentalAnalyses")
#Final morph pca with row names
#PCA using all variables, all variables not correlating, dataset imported from xlsx in rstudio
reduced_means_July_15_1 <- read.csv("reduced_means_July_15_1.csv", sep=",", header=T, stringsAsFactors = FALSE)
rownames(reduced_means_July_15_1) <-  c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10", "W11", "W12", "WAUX3", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E10") 

reduced_means_July_15_1.pca <- prcomp(reduced_means_July_15_1[,c(3:19)], center = TRUE,scale. = TRUE)

#grouping for west and east pops
two_groups <- c(rep("West",13), rep("East",8))

m<- ggbiplot(reduced_means_July_15_1.pca, ellipse = FALSE, obs.scale = 1.25, var.scale = 1.5, varname.adjust = 1.5, 
             var.axes = FALSE, labels= NULL, groups=two_groups)+
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_point(aes(shape=two_groups, color=two_groups, size=2))+
  geom_text_repel(aes(label = rownames(reduced_means_July_15_1)), color = "black", size=4.5)+
  scale_colour_manual(values = c("red", "blue"))+
  xlab(expression("Axis 1 (22.9%)"))+
  ylab(expression("Axis 2 (17.9%)"))+
  theme_classic(base_size = 10) +
  ggtitle("(B) Morphological") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black")) +
  theme(legend.position = "none") +
  scale_x_reverse()
m

m$layers <- c(m$layers, m$layers[[3]])
m$layers <- c(m$layers, m$layers[[1]])
m


#----------#
#Genetic PCA
#----------#

QH_eigen <- read.csv("QH_eigen.csv", sep=",", header=T, stringsAsFactors = FALSE)
rownames(QH_eigen)<- c("E1", "E10", "E13", "E15", "E2", "E3", "E4", "E5", "E6", "E7",
                       "W1", "W10", "W11", "W12", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "WAUX3")
colnames(QH_eigen) <- c("Pop","1", "2","3","4","5","6","7", "8","9","10","11","12","13", "14","15")
point_colors<-c(rep("red",10),rep("blue",13))
point_shapes<-c(rep(19,10),rep(24,13))
plot_col<-4
g <- ggplot(QH_eigen) +
  geom_hline(yintercept = 0, color = "grey")+
  geom_vline(xintercept = 0, color = "grey")+
  geom_point(aes(x =QH_eigen$`1`, y = QH_eigen$`2`),
             shape=point_shapes, color=point_colors, fill=point_colors, size=4.5)+
  xlab(expression("Axis 1 (63.3%)")) + 
  ylab(expression("Axis 2 (13.9%)"))+
  theme_classic(base_size = 10) +
  ggtitle("(A) Genetic") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black")) +
  theme(legend.position = "none")
g<- g + geom_text_repel(aes(x =QH_eigen$`1`, y = QH_eigen$`2`),
                        label = rownames(QH_eigen), color = "black",nudge_x = 0.02, nudge_y = 0.04, size=4.5)
g


#----------#
#Combining all PCA's
#----------#


library(egg)

ggarrange(g, m, e, widths = c(1.5,1.5,1.5))


setwd("C:/Users/shoban/Downloads/wc2.0_2.5m_bio")
 files <- list.files(pattern='tif', full.names=TRUE) #Load climate files
 bioclim2.5 <- stack(files) #Create a raster stack
 setwd("C:/Users/shoban/Documents/GitHub/QH_EnvironmentalAnalyses/")
 QHloc <- read.csv("Hob_QHOccur_gen.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations
data(wrld_simpl)
clim <- extract(bioclim2.5, QHloc[,1:2]) 
clim<-clim[,-c(1,5,6,13,14,16,18)]

 setwd("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Qhavardii/7indiv_11loci/")
gen_sum_stats<-read.csv("Qha_summ_stats.csv") 

#summary(lm(gen_sum_stats[,5]~clim[,1]))
reg_res<-matrix(nrow=dim(clim)[2],ncol=5)
stat_list<-c(5,6,8,9,10)

for (ss in 1:length(stat_list)){
	for (c in 1:dim(clim)[2]){
		p_val<-summary(lm(gen_sum_stats[,stat_list[ss]]~clim[,c]))[4][[1]][8]
		reg_res[c,ss]<-p_val
	}
}
 matrix(p.adjust(reg_res,"BH"),nrow=12,ncol=5)
 
 which(p.adjust(reg_res,"BH")<0.05)%%12
 [1]  2  3 10  2  3  8  9 10  1  2  8 10
 
 
  setwd("C:/Users/shoban/Documents/GitHub/QH_EnvironmentalAnalyses/")
list_bioc<-read.csv("list_bioclim_var.csv")
 list_bioc<-list_bioc[list_bioc[,1]=="x"]
 list_bioc<-list_bioc[-which(list_bioc[,1]=="x"),]
 
 par(mfrow=c(4,3),mar=c(3,3,2,2),oma=c(2,3,1,1))
 plot(gen_sum_stats[,6]~clim[,2])
 plot(gen_sum_stats[,6]~clim[,3])
 plot(gen_sum_stats[,6]~clim[,10])
 plot(gen_sum_stats[,8]~clim[,2])
 plot(gen_sum_stats[,8]~clim[,3])
 plot(gen_sum_stats[,8]~clim[,8])
 plot(gen_sum_stats[,8]~clim[,9])
 plot(gen_sum_stats[,8]~clim[,10])
 plot(gen_sum_stats[,9]~clim[,1])
 plot(gen_sum_stats[,9]~clim[,2])
 plot(gen_sum_stats[,10]~clim[,8])
 plot(gen_sum_stats[,10]~clim[,10])
 
 matrix(p.adjust(reg_res,"BY"),nrow=12,ncol=5)
 pdf(width=10, height=4,file="regr_gen_clim.pdf")
 par(mfrow=c(1,3))

 plot(gen_sum_stats[,6]~clim[,10],col=c(rep("red",13),rep("blue",13)),pch=19, xlab=list_bioc[10,3], ylab=colnames(gen_sum_stats)[6])
	text((gen_sum_stats[,6]-.005)~clim[,10], labels = gen_sum_stats[,2])
 abline(lm(gen_sum_stats[,6]~clim[,10]))
 	  text(45,.72,paste("R2= ",round(unlist(summary(lm(gen_sum_stats[,6]~clim[,10]))[8]),3),sep=""))
plot(gen_sum_stats[,8]~clim[,10],col=c(rep("red",13),rep("blue",13)),pch=19, xlab=list_bioc[10,3], ylab=colnames(gen_sum_stats)[8])
 	text((gen_sum_stats[,8]-.5)~clim[,10], labels = gen_sum_stats[,2])
	abline(lm(gen_sum_stats[,8]~clim[,10]))
 	  text(50,58,paste("R2= ",round(unlist(summary(lm(gen_sum_stats[,8]~clim[,10]))[8]),3),sep=""))
plot(gen_sum_stats[,9]~clim[,2],col=c(rep("red",13),rep("blue",13)),pch=19, xlab=list_bioc[2,3], ylab=colnames(gen_sum_stats)[9])
	text((gen_sum_stats[,9]+.001)~clim[,2], labels = gen_sum_stats[,2])
	abline(lm(gen_sum_stats[,9]~clim[,2]))
	  text(40,.13,paste("R2= ",round(unlist(summary(lm(gen_sum_stats[,9]~clim[,2]))[8]),3),sep=""))
	dev.off()
	
	

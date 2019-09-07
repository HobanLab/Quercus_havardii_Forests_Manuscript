##################
# Get climate variables 
##################
library(stats); library(dismo); library(maptools); library(rgdal)

setwd("C:/Users/shoban/Downloads/wc2.0_2.5m_bio")
setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Downloads/wc2.0_2.5m_bio")

 files <- list.files(pattern='tif', full.names=TRUE) #Load climate files
 bioclim2.5 <- stack(files) #Create a raster stack

setwd("C:/Users/shoban/Documents/GitHub/QH_EnvironmentalAnalyses/Hoban_work/")
setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Documents/git/QH_EnvironmentalAnalyses/Hoban_work/")

 QHloc <- read.csv("QH_Pops_Final.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations

  clim <- extract(bioclim2.5, QHloc[,1:2]) 

#add in lat/long- optional... I would argue that lat long are not biologically relevant here
#clim_and_loc<-cbind(clim, QHloc[,1:2])
 
 #read in the names of the bioclim variables for the axis names of the plots
 #the x is for those variables to be removed due to correlation
 list_bioc<-read.csv("list_bioclim_var.csv",header=F)

 
#################
# Get genetic data
##################
gen_sum_stats<-read.csv("Qha_summ_stats.csv") 
#these are the genetic statistics: number clones (5). heterozygosity (6), allelic richness (8), pw FST (9), relatedness (10)
stat_list<-c(6,8,10)



##################
# Do Regressions
################

pop_set<-list(1:10,11:23,1:23,12:23)
point_colors<-c(rep("red",10),rep("blue",13))
point_shapes<-c(rep(21,10),rep(24,13))
plot_col<-4

#Note: The following analysis assumes we ran correlations for populations within each region and picked out
	#the least correlated variables within each region to include
	#We originally looked at removing correlated variables based on all the populations 
	#but ran into the problem of two perfectly correlated variables in east
	#to do this, subtract the correlated variables clim<-clim[,-c(1,5,6,13,14,16,18)]
	#which would mean c(2:4,7:10,12,15,17,19); I did this and saved such results under folder "old" 
	
#Loop to go over East, West, East+West Together, and West without W-1
for (p in 1:4){
	#store name of region we are in for writing file names
	
	if (p==1) { region<-"E"; clim_to_keep<-c(1,3,4,6,11,15,18); clim_reg<-clim[,clim_to_keep] }
	if (p==2) { region<-"W"; clim_to_keep<-c(2,3,7,9,10,14,17,18,19); clim_reg<-clim[,clim_to_keep] }
	if (p==3) { region<-"overall"; clim_to_keep<-c(2:4,7:10,12,15,17,19); clim_reg<-clim[,clim_to_keep] }
	if (p==4) { region<-"W_sub"; clim_to_keep<-c(2,3,7,9,10,14,17,18,19); clim_reg<-clim[,clim_to_keep] }
	
	#place to keep regression results (p values) of all regressions run
	reg_pval<-matrix(nrow=dim(clim_reg)[2],ncol=3)
	reg_r2<-matrix(nrow=dim(clim_reg)[2],ncol=3)

	#linear models
	for (ss in 1:length(stat_list)){
		for (c in 1:dim(clim_reg)[2]){
			reg_pval[c,ss]<-summary(lm(gen_sum_stats[pop_set[[p]],stat_list[ss]]~clim_reg[pop_set[[p]],c]))[4][[1]][8]
			reg_r2[c,ss]<-unlist(summary(lm(gen_sum_stats[pop_set[[p]],stat_list[ss]]~clim_reg[pop_set[[p]],c]))[8]	)
		}
	}
	
	#identify rows and columns with significant p values 
	 	row_col_signif<-which( matrix(p.adjust(reg_pval,"none"),nrow=length(clim_to_keep),ncol=3)<=0.05,arr.ind=T)
	if (length(row_col_signif)!=0){
	#store the list of which variables are significant and what genetic statistic
		variables_id<-matrix(nrow=length(row_col_signif[,1]),ncol=2)
		for (l in 1:length(row_col_signif[,1])) variables_id[l,1]<-as.character(list_bioc[row_col_signif[l,1],3])
		for (l in 1:length(row_col_signif[,1])) variables_id[l,2]<-colnames(gen_sum_stats)[stat_list[row_col_signif[l,2]]]
		write.csv(variables_id[order(variables_id[,1]),],file=paste0("clim_variables_signif_",region,"loose.csv"))
	 
	 sort(reg_r2)
		 
	 #this will cycle through all signification associations (row_col_signif) and plot that regression, pulling from row_col_signif
		if (p==1) { pdf(width=10, height=5,file=paste0("regr_gen_clim_none_",region,".pdf")); par(mfrow=c(1,3))}
		if (p==2) { pdf(width=10, height=8,file=paste0("regr_gen_clim_none_",region,".pdf")); par(mfrow=c(2,3))}
		if (p==3) { pdf(width=12, height=9,file=paste0("regr_gen_clim_none_",region,".pdf")); par(mfrow=c(4,4),mar=c(4,4,2,2),oma=c(2,3,1,1))}
		if (p==4) { pdf(width=10, height=8,file=paste0("regr_gen_clim_none_",region,".pdf")); par(mfrow=c(1,2))}

		for (i in 1:nrow(row_col_signif)){
			eq1<-gen_sum_stats[pop_set[[p]],stat_list[row_col_signif[i,2]]]~clim_reg[pop_set[[p]],row_col_signif[i,1]]
			 plot(eq1,col=point_colors[pop_set[[p]]],pch=point_shapes[pop_set[[p]]], 
			 xlab=list_bioc[clim_to_keep[row_col_signif[i,1]],3], ylab=colnames(gen_sum_stats)[stat_list[row_col_signif[i,2]]])
			 abline(lm(eq1))
			mtext(paste(paste("R2= ",round(unlist(summary(lm(eq1))[8]),3),sep=""),paste("p= ",round(summary(lm(eq1))[4][[1]][8],3),sep=""),sep="  "),side=3)
			  #text(eq1, labels = gen_sum_stats[pop_set[[p]],2]) #optional to add labels for pop names
		 }
		dev.off() 
		
	}
	
	
	
	#identify rows and columns with significant p values 
	 	row_col_signif<-which( matrix(p.adjust(reg_pval,"BH"),nrow=length(clim_to_keep),ncol=3)<=0.054,arr.ind=T)
	if (length(row_col_signif)!=0){
	#store the list of which variables are significant and what genetic statistic
		variables_id<-matrix(nrow=length(row_col_signif[,1]),ncol=2)
		for (l in 1:length(row_col_signif[,1])) variables_id[l,1]<-as.character(list_bioc[row_col_signif[l,1],3])
		for (l in 1:length(row_col_signif[,1])) variables_id[l,2]<-colnames(gen_sum_stats)[stat_list[row_col_signif[l,2]]]
		write.csv(variables_id[order(variables_id[,1]),],file=paste0("clim_variables_signif_",region,"loose.csv"))
	 
	 sort(reg_r2)
		 
	 #this will cycle through all signification associations (row_col_signif) and plot that regression, pulling from row_col_signif
		if (p==1) { pdf(width=9, height=5,file=paste0("regr_gen_clim_loose_",region,".pdf")); par(mfrow=c(1,2))}
		if (p==2) { pdf(width=5, height=5,file=paste0("regr_gen_clim_loose_",region,".pdf")); par(mfrow=c(1,1))}
		if (p==3) { pdf(width=12, height=9,file=paste0("regr_gen_clim_loose_",region,".pdf")); par(mfrow=c(5,4),mar=c(4,4,2,2),oma=c(2,3,1,1))}
		for (i in 1:nrow(row_col_signif)){
			eq1<-gen_sum_stats[pop_set[[p]],stat_list[row_col_signif[i,2]]]~clim_reg[pop_set[[p]],row_col_signif[i,1]]
			 plot(eq1,col=point_colors[pop_set[[p]]],pch=point_shapes[pop_set[[p]]], 
			 xlab=list_bioc[clim_to_keep[row_col_signif[i,1]],3], ylab=colnames(gen_sum_stats)[stat_list[row_col_signif[i,2]]])
			 abline(lm(eq1))
			 mtext(paste(paste("R2= ",round(unlist(summary(lm(eq1))[8]),3),sep=""),paste("p= ",round(summary(lm(eq1))[4][[1]][8],3),sep=""),sep="  "),side=3)
			  #text(eq1, labels = gen_sum_stats[pop_set[[p]],2]) #optional to add labels for pop names
		 }
		dev.off() 
		
	}

	
	#identify rows and columns with significant p values 
		row_col_signif<-which( matrix(p.adjust(reg_pval,"BY"),nrow=length(clim_to_keep),ncol=3)<=0.054,arr.ind=T)
	 if (length(row_col_signif)!=0){
	 #store the list of which variables are significant and what genetic statistic
		variables_id<-matrix(nrow=length(row_col_signif[,1]),ncol=2)
		for (l in 1:length(row_col_signif[,1])) variables_id[l,1]<-as.character(list_bioc[row_col_signif[l,1],3])
		for (l in 1:length(row_col_signif[,1])) variables_id[l,2]<-colnames(gen_sum_stats)[stat_list[row_col_signif[l,2]]]
		write.csv(variables_id[order(variables_id[,1]),],file=paste0("clim_variables_signif_",region,"strict.csv"))
	 
		if (p==2) { pdf(width=5, height=5,file=paste0("regr_gen_clim_strict_",region,".pdf")); par(mfrow=c(1,1))}
		if (p==3) { pdf(width=11, height=9,file=paste0("regr_gen_clim_strict_",region,".pdf")); par(mfrow=c(3,3),mar=c(4,4,2,2),oma=c(2,3,1,1))}
		for (i in 1:nrow(row_col_signif)){
			 eq1<-gen_sum_stats[pop_set[[p]],stat_list[row_col_signif[i,2]]]~clim_reg[pop_set[[p]],row_col_signif[i,1]]
			 plot(eq1,col=point_colors[pop_set[[p]]],pch=point_shapes[pop_set[[p]]], 
			 xlab=list_bioc[clim_to_keep[row_col_signif[i,1]],3], ylab=colnames(gen_sum_stats)[stat_list[row_col_signif[i,2]]])
			 abline(lm(eq1))
			 mtext(paste(paste("R2= ",round(unlist(summary(lm(eq1))[8]),3),sep=""),paste("p= ",round(summary(lm(eq1))[4][[1]][8],3),sep=""),sep="  "),side=3)
			  #text(eq1, labels = gen_sum_stats[pop_set[[p]],2]) #optional to add labels for pop names
		 }
	dev.off() 
	}
}


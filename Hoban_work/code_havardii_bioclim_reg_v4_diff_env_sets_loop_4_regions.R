##################
# Get climate variables 
##################
library(stats); library(dismo); library(maptools); library(rgdal)

setwd("C:/Users/shoban/Downloads/wc2.0_2.5m_bio")
 files <- list.files(pattern='tif', full.names=TRUE) #Load climate files
 bioclim2.5 <- stack(files) #Create a raster stack

setwd("C:/Users/shoban/Documents/GitHub/QH_EnvironmentalAnalyses/Hoban_work/")
 QHloc <- read.csv("Hob_QHOccur_gen.csv", sep=",", header=T, stringsAsFactors = FALSE) #Loading locations

clim <- extract(bioclim2.5, QHloc[,1:2]) 

#add in lat/long- optional... I would argue that lat long are not biologically relevant here
#clim_and_loc<-cbind(clim, QHloc[,1:2])
 
setwd("C:/Users/shoban/Documents/GitHub/QH_EnvironmentalAnalyses/Hoban_work/")
 #read in the names of the bioclim variables for the axis names of the plots
 #the x is for those variables to be removed due to correlation
 list_bioc<-read.csv("list_bioclim_var.csv")
 list_bioc<-list_bioc[-which(list_bioc[,1]=="x"),]

 
#################
# Get genetic data
##################
setwd("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Qhavardii/7indiv_11loci/")
gen_sum_stats<-read.csv("Qha_summ_stats.csv") 
#these are the genetic statistics: number clones (5). heterozygosity (6), allelic richness (8), pw FST (9), relatedness (10)
stat_list<-c(5,6,8,9,10)

setwd("C:/Users/shoban/Documents/GitHub/QH_EnvironmentalAnalyses/Hoban_work/")


##################
# Do Regressions
################

pop_set<-list(1:13,14:26,1:26,15:26)
point_colors<-c(rep("red",13),rep("blue",13))
point_shapes<-c(rep(21,13),rep(24,13))
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
	reg_pval<-matrix(nrow=dim(clim_reg)[2],ncol=5)
	reg_r2<-matrix(nrow=dim(clim_reg)[2],ncol=5)

	#linear models
	for (ss in 1:length(stat_list)){
		for (c in 1:dim(clim_reg)[2]){
			reg_pval[c,ss]<-summary(lm(gen_sum_stats[pop_set[[p]],stat_list[ss]]~clim_reg[pop_set[[p]],c]))[4][[1]][8]
			reg_r2[c,ss]<-unlist(summary(lm(gen_sum_stats[pop_set[[p]],stat_list[ss]]~clim_reg[pop_set[[p]],c]))[8]	)
		}
	}
	
	#identify rows and columns with significant p values 
	 	row_col_signif<-which( matrix(p.adjust(reg_pval,"BH"),nrow=length(clim_to_keep),ncol=5)<0.05,arr.ind=T)
		if (p==4) row_col_signif<-which( matrix(p.adjust(reg_pval,"BH"),nrow=length(clim_to_keep),ncol=5)<0.075,arr.ind=T)
	#store the list of which variables are significant and what genetic statistic
		variables_id<-matrix(nrow=length(row_col_signif[,1]),ncol=2)
		for (l in 1:length(row_col_signif[,1])) variables_id[l,1]<-as.character(list_bioc[row_col_signif[l,1],3])
		for (l in 1:length(row_col_signif[,1])) variables_id[l,2]<-colnames(gen_sum_stats)[stat_list[row_col_signif[l,2]]]
		write.csv(variables_id[order(variables_id[,1]),],file=paste0("clim_variables_signif_",region,"loose.csv"))
	 
	 sort(reg_r2)
		 
	 #this will cycle through all signification associations (row_col_signif) and plot that regression, pulling from row_col_signif
	if (p!=2) {
		pdf(width=10, height=10,file=paste0("regr_gen_clim_loose_",region,".pdf"))
		par(mfrow=c(ceiling(nrow(row_col_signif)/4),plot_col),mar=c(4,4,2,2),oma=c(2,3,1,1))
		for (i in 1:nrow(row_col_signif)){
			eq1<-gen_sum_stats[pop_set[[p]],stat_list[row_col_signif[i,2]]]~clim_reg[pop_set[[p]],row_col_signif[i,1]]
			 plot(eq1,col=point_colors[pop_set[[p]]],pch=point_shapes[pop_set[[p]]], xlab=list_bioc[row_col_signif[i,1],3], ylab=colnames(gen_sum_stats)[stat_list[row_col_signif[i,2]]])
			 abline(lm(eq1))
			 mtext(paste("R2= ",round(unlist(summary(lm(eq1))[8]),3),sep=""),side=3)
			  #text(eq1, labels = gen_sum_stats[pop_set[[p]],2]) #optional to add labels for pop names
		 }
		dev.off() 
	}

	
	#identify rows and columns with significant p values 
		row_col_signif<-which( matrix(p.adjust(reg_pval,"BY"),nrow=length(clim_to_keep),ncol=5)<0.05,arr.ind=T)
	 if (length(row_col_signif)!=0){
	 #store the list of which variables are significant and what genetic statistic
		variables_id<-matrix(nrow=length(row_col_signif[,1]),ncol=2)
		for (l in 1:length(row_col_signif[,1])) variables_id[l,1]<-as.character(list_bioc[row_col_signif[l,1],3])
		for (l in 1:length(row_col_signif[,1])) variables_id[l,2]<-colnames(gen_sum_stats)[stat_list[row_col_signif[l,2]]]
		write.csv(variables_id[order(variables_id[,1]),],file=paste0("clim_variables_signif_",region,"strict.csv"))
	 
		pdf(width=10, height=4,file=paste0("regr_gen_clim_strict_",region,".pdf"))
		par(mfrow=c(ceiling(nrow(row_col_signif)/4),plot_col),mar=c(4,4,2,2),oma=c(2,3,1,1))
		for (i in 1:nrow(row_col_signif)){
			 eq1<-gen_sum_stats[pop_set[[p]],stat_list[row_col_signif[i,2]]]~clim_reg[pop_set[[p]],row_col_signif[i,1]]
			 plot(eq1,col=point_colors[pop_set[[p]]],pch=point_shapes[pop_set[[p]]], xlab=list_bioc[row_col_signif[i,1],3], ylab=colnames(gen_sum_stats)[stat_list[row_col_signif[i,2]]])
			 abline(lm(eq1))
			 mtext(paste("R2= ",round(unlist(summary(lm(eq1))[8]),3),sep=""),side=3)
			  #text(eq1, labels = gen_sum_stats[pop_set[[p]],2]) #optional to add labels for pop names
		 }
	dev.off() 
	}
}


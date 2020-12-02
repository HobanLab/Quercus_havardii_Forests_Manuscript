#Downloaded data from Binned_MP1_All_Samples.csv from Drive
#manual replacement of all cells with text in them (i.e. no alleles) with -1, saved as binned_dat.csv


##############
# SPIN UP	#
##############

setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Qhavardii")
setwd("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Qhavardii")
setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Documents/git/Qhavardii/")
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Qhavardii")
library("adegenet"); library("poppr"); library(hierfstat)

repool_new<- function(genind_obj,vect_pops){
	genind_obj_sep<-seppop(genind_obj)
	genind_obj_merge<-genind_obj_sep[[vect_pops[1]]]
	for (i in 1:(length(vect_pops)-1)) genind_obj_merge<-repool(genind_obj_merge,genind_obj_sep[[vect_pops[i+1]]])
	genind_obj_merge
}
remove3and2<- function(x) {substr(x,4,(nchar(x)-2))}
 
##############ADD COMMENT HERE EXPLAINING FILE TYPES FOR THE CODE##############

#to do for all subsets of the data

subset_names<-c("7indiv_7loci","7indiv_11loci_RedefinedPopsMay2019","10indiv_7loci","10indiv_11loci") 
min_indivs<-c(7,7,10,10)
for (subn in 2:2) {
#setwd(paste("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Qhavardii/",subset_names[subn],sep=""))

#------------------------------------------------------------------------------#
#				SECTION ONE
#------------------------------------------------------------------------------#
	
##############
# CLONES	#
##############
	
#import the data- note individuals with no data at all are dropped
ade_test<-read.genepop(paste("QH_GenePop_",subset_names[subn],".gen",sep=""),ncode=3)
pop(ade_test)<-remove3and2(as.character(ade_test@pop))
pop_names<-levels(ade_test@pop)

#First put into poppr format
popr_test <- as.genclone(ade_test)
strata(popr_test) <- other(popr_test)$population_hierarchy[-1]
list_a<-mlg.id(popr_test)
#Function to pull out individual indices where clone length greater than 1
clone_index<-which(sapply(list_a,function(x) length(x)>1))
list_a[clone_index]

#TO DO SOME KIND OF MAPPING OF CLONES WITHIN POPULATIONS

#This removes clones and then saves as new file for Genealex
popr_nocl<-clonecorrect(popr_test,strata=~Pop)
genind2genalex(genclone2genind(popr_nocl),file=paste("QH_clone_free_",subset_names[subn],".csv",sep=""))

#Create genpop and genind objects that now have no clones- GI_nocl, GP_nocl
GI_nocl<-genclone2genind(popr_nocl); 	GP_nocl<-genind2genpop(GI_nocl)

						  
###################
#  BASIC STATS	  #
###################
	
#narrow down to populations with >=10 individuals- call these GP_sub, GI_sub for subset
#pop_keep<- which(as.vector(table(GI_nocl@pop)>=10))
pop_keep<- which(as.vector(table(GI_nocl@pop)>=min_indivs[subn]))
GI_sub<-repool_new(GI_nocl,pop_keep);	GP_sub<-GP_nocl[pop_keep,]
samp_size<-table(GI_nocl@pop)[pop_keep]

#now we are down to no clones and decently sampled populations, with minimal missing data!
#Summary statistics
par(mfrow=c(2,1))
Qha_sumstats<-summary(GI_sub)
Qha_alleles<-Qha_sumstats$pop.n.all
Qha_all_rich<-colSums(allelic.richness(GI_sub)$Ar)
Qha_N_v_MLG<-poppr(popr_test)[pop_keep,2:3] #note this is not popr_nocl because we want the sample size vs. MLG
num_clones<-Qha_N_v_MLG[,1]-Qha_N_v_MLG[,2]
Qha_Hexp<-poppr(popr_nocl)[pop_keep,10]

#Calculating Hobs and Hexp per population by averaging loci
n.pop <- seppop(popr_nocl) 
mean.hobs <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs)))
mean.hobs[is.nan(mean.hobs)] <- NA 
barplot(mean.hobs)

mean.hexp <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hexp)))
mean.hexp[is.nan(mean.hexp)] <- NA 
barplot(mean.hexp)

East_pops<-which(substr(pop_names[pop_keep],1,1)=="E")
West_pops<-which(substr(pop_names[pop_keep],1,1)=="W")

###################
#  PAIRWISE FST	  #
###################

sm_fst_mat<-as.matrix(pairwise.fst(GI_sub))
rownames(sm_fst_mat)<-pop_names[pop_keep];	colnames(sm_fst_mat)<-pop_names[pop_keep]
sm_fst_mat[sm_fst_mat==0]<-NA
Qha_pwfst<-apply(sm_fst_mat,2,mean,na.rm=T)
#write.csv(sm_fst_mat,file="fst_mat.csv")

sm_fst_mat[sm_fst_mat==0]<-NA
pdf(file="Qha_EvW_fst.pdf")
boxplot(c(sm_fst_mat[East_pops,East_pops]),c(sm_fst_mat[West_pops,East_pops]),c(sm_fst_mat[West_pops,West_pops]),names=c("w/i E", "b/w E.W", "w/i W"),ylab="pairwise FST")
dev.off()
pdf(file="Qha_pwfst.pdf")
boxplot(sm_fst_mat,na.rm=T,names=pop_names[pop_keep],ylab="pairwise FST by population",xlab="",las=2)
mtext(side=1,line=4,"populations names")
dev.off()
write.csv(sm_fst_mat,file="fst_mat.csv")

########################
# Mantel test 
########################
#setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Documents/git/Qhavardii/")
Qha_gd<-read.csv("GeographicDist_SquareMatrix.csv")
Qha_gd<-as.matrix(Qha_gd[,-24])
Qha_gd_E<-Qha_gd[1:10,1:10]; Qha_gd_W<-Qha_gd[11:23,11:23]
Qha_fst_mat<-read.csv("fst_mat.csv")
Qha_fst_mat<-as.matrix(Qha_fst_mat[,-1])
Qha_fst_mat_E<-Qha_fst_mat[1:10,1:10]; Qha_fst_mat_W<-Qha_fst_mat[11:23,11:23]

pdf(file="qh_overall_mantel.pdf")
plot(Qha_fst_mat,Qha_gd, xlab="geographic distance (km)", ylab="pairwise FST")
lm0<-lm(c(Qha_fst_mat)~c(Qha_gd))
abline(lm0)
text(50,.12,paste("R2=", round(unlist(summary(lm0)[8]),3)))
dev.off()
summary(lm(c(Qha_fst_mat)~c(Qha_gd)))

pdf(file="qhe_mantel.pdf")
plot(Qha_fst_mat[1:10,1:10],Qha_gd[1:10,1:10], xlab="geographic distance (km)", ylab="pairwise FST")
dev.off()
summary(lm(c(Qha_fst_mat[1:10,1:10])~c(Qha_gd[1:10,1:10])))

pdf(file="qhw_mantel.pdf")
plot(Qha_fst_mat[11:23,11:23]~Qha_gd[11:23,11:23], xlab="geographic distance (km)", ylab="pairwise FST")
lm1<-lm(c(Qha_fst_mat[11:23,11:23])~c(Qha_gd[11:23,11:23]))
abline(lm1)
text(50,.12,paste("R2=", round(unlist(summary(lm1)[8]),3)))
dev.off()
summary(lm(c(Qha_fst_mat[11:23,11:23])~c(Qha_gd[11:23,11:23]),na.omit=T))


Qha_gd<-Qha_gd[lower.tri(Qha_gd)]
Qha_fst_mat<-Qha_fst_mat[lower.tri(Qha_fst_mat)]
plot(Qha_fst_mat~Qha_gd)
Qha_gd_E<-Qha_gd_E[upper.tri(Qha_gd_E)]
Qha_fst_mat_E<-Qha_fst_mat_E[upper.tri(Qha_fst_mat_E)]
Qha_gd_W<-Qha_gd_W[upper.tri(Qha_gd_W)]
Qha_fst_mat_W<-Qha_fst_mat_W[upper.tri(Qha_fst_mat_W)]
summary(lm(Qha_fst_mat_E~Qha_gd_E)); plot(Qha_fst_mat_E~Qha_gd_E)
summary(lm(Qha_fst_mat_W~Qha_gd_W)); plot(Qha_fst_mat_W~Qha_gd_W)


######################
# Pop size playing with
######################

summ_stats<-read.csv("7indiv_11loci/Qha_summ_stats.csv")
psizes<-read.csv("C:/Users/shoban.DESKTOP-DLPV5IJ/Documents/git/Qhavardii/pop_size_estim.csv")
boxplot(summ_stats[,8]~psizes[,3])
anova(lm(summ_stats[,8]~psizes[,3]))
anova(lm(summ_stats[14:26,9]~psizes[14:26,3]))
boxplot((summ_stats[14:26,9]~psizes[14:26,3]))
anova(lm(summ_stats[1:13,9]~psizes[1:13,3]))
boxplot((summ_stats[1:13,9]~psizes[1:13,3]))
t.test(summ_stats[,8]~psizes[,4])


###################
#  RELATEDNESS	  #
###################
#Note the relatedness file has a different format- essentially just delete the top two lines of a genalex file, add in column names for loci (every column) and that is all.  They are saved as "_rel.csv" to distinguish
library("Demerelate");		Qha_rel<-read.csv("QH_GenePop_REL_7Indiv_11loci_RedefinedPopsMay2019.csv",header=T)
wang.results<-Demerelate(Qha_rel,object=T,value="wang")
#means of pairwise relatedness within each population
wang.means<-unlist(lapply(wang.results$Empirical_Relatedness,mean))
						  
#calculate proportion of population that are half sibs or greater, and plot this
half_sib<-function(x) (sum(x>0.25)/length(x))
mean(unlist(lapply(wang.results$Empirical_Relatedness,half_sib))[1:10])
mean(unlist(lapply(wang.results$Empirical_Relatedness,half_sib))[11:22])
# 0.381 vs 0.1866
pdf(file="prop_half_sib.pdf"); barplot(unlist(lapply(wang.results$Empirical_Relatedness,half_sib)),las=2,ylab="average pairwise relatedness value",col=c(rep("firebrick1",10),rep("cornflowerblue",13))); abline(h=0.25, lty=3, col="red"); 
mtext(side=1,line=4,"populations names"); text(20,.4,"red line is expectation for full siblings",col="red",cex=.8)
dev.off()

			
#####################################						  
#  OUTPUT SUMMARY TABLE, DO T TESTS #
#####################################						  
write.csv(file="Qha_summ_stats.csv", cbind(pop_names[pop_keep],Qha_N_v_MLG,num_clones,Qha_Hexp,Qha_alleles,Qha_all_rich,Qha_pwfst,wang.means))

write.csv(file="normality_pval.csv", cbind(as.numeric(shapiro.test(Qha_alleles)[2]), as.numeric(shapiro.test(Qha_all_rich)[2]), as.numeric(shapiro.test(Qha_Hexp)[2]), as.numeric(shapiro.test(Qha_pwfst)[2]), as.numeric(shapiro.test(wang.means)[2])))

pdf("normality_plots.pdf")
par(mfrow=c(3,2))
qqnorm(Qha_alleles); qqline(Qha_alleles); qqnorm(Qha_all_rich); qqline(Qha_all_rich)
qqnorm(Qha_Hexp); qqline(Qha_Hexp); qqnorm(Qha_pwfst); qqline(Qha_pwfst); qqnorm(wang.means); qqline(wang.means)
dev.off()

t_test_res<-matrix(ncol=4,nrow=4)
t_test_res[1,]<-c(unlist(t.test(Qha_all_rich[East_pops],Qha_all_rich[West_pops])[c(3,5)]),
	as.numeric(wilcox.test(Qha_all_rich[East_pops],Qha_all_rich[West_pops])[3]))
t_test_res[2,]<-c(unlist(t.test(Qha_Hexp[East_pops],Qha_Hexp[West_pops])[c(3,5)]),
	as.numeric(wilcox.test(Qha_Hexp[East_pops],Qha_Hexp[West_pops])[3]))
t_test_res[3,]<-c(unlist(t.test(Qha_pwfst[East_pops],Qha_pwfst[West_pops])[c(3,5)]),
	as.numeric(wilcox.test(Qha_pwfst[East_pops],Qha_pwfst[West_pops])[3]))
t_test_res[4,]<-c(unlist(t.test(wang.means[East_pops],wang.means[West_pops])[c(3,5)]),
	as.numeric(wilcox.test(wang.means[East_pops],wang.means[West_pops])[3]))
colnames(t_test_res)<-c("p val T", "East mean", "West mean", "p val W"); rownames(t_test_res)<-c("Ar","Hexp","pw fst","relatedness")
write.csv(t_test_res,file="t_test_results.csv")
						  

#------------------------------------------------------------------------------#
#				SECTION TWO
#------------------------------------------------------------------------------#
							  
#################
#	CLUSTERING	#
#################

#TO DO CONSIDER ARRANGING POPULATIONS IN ORDER OF EAST TO WEST OR SOMETHING LIKE THAT- CLOSE POPULATIONS NEAR
#pdf(file="cluster_number_K.pdf"); grp<-find.clusters(GI_sub,max.n.clust=40,n.pca=100); dev.off()

#############
#	DAPC	#
#############
#pdf(file="Qha_dapc1.pdf")
#for (cl in 2:10){
#grp<-find.clusters(GI_sub,max.n.clust=12,n.pca=50,n.clust=cl)
#dapc1 <- dapc(GI_sub, grp$grp,n.pca=50,n.da=50);	scatter(dapc1)
#}
#dev.off()
pdf(file="Qha_str_like.pdf",width=40,height=9)
for (cl in 2:10){
grp<-find.clusters(GI_sub,max.n.clust=12,n.pca=50,n.clust=cl)
dapc1 <- dapc(GI_sub, grp$grp,n.pca=50,n.da=50);	compoplot(dapc1,show.lab=T)
}
dev.off()
pdf(file="Qha_dapc2.pdf")
dapc2 <- dapc(GI_sub,n.pca=50,n.da=50);	scatter(dapc2)
dev.off()
pdf("Qha_snapclust.pdf", height=5, width=18)
for (cl in 2:10){
grp<-find.clusters(GI_sub,max.n.clust=12,n.pca=50,n.clust=cl)
res<-snapclust(GI_sub,cl,pop.ini=grp$grp);		compoplot(res)
}
dev.off()


#############
#	PCA/CA	#
#############
#pdf(file="Qha_pca_bw.pdf")
#pca1<-dudi.pca(df = tab(GI_sub, NA.method = "zero"), scannf = FALSE, nf = 50)
#s.class(pca1$li,pop(GI_sub),xax=1,yax=2,sub="PCA 1-2",csub=2)
#dev.off()
#pdf(file="Qha_pca_col.pdf")
#col <- funky(length(pop(GI_sub)))
#s.class(pca1$li, pop(GI_sub),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,cstar=0, cpoint=3, grid=FALSE)
#dev.off()
pdf(file="Qha_ca.pdf")
ca1 <- dudi.coa(tab(GP_sub),scannf=FALSE,nf=3)
s.label(ca1$li, sub="CA 1-2",csub=2)
dev.off()


#####################
#	EAST & WEST SEP	#
#####################

GI_sub_E<-repool_new(GI_sub,East_pops);		GP_sub_E<-genind2genpop(GI_sub_E)
levels(GI_sub_E@pop)<-levels(GI_sub@pop)[East_pops]
pdf(file="QhaE_dapc2.pdf")
dapc2 <- dapc(GI_sub_E, n.pca=50,n.da=50);	scatter(dapc2)
dev.off()
pdf("QhaE_snapclust.pdf", height=5, width=18)
for (cl in 2:10){
grp<-find.clusters(GI_sub_E,max.n.clust=12,n.pca=50,n.clust=cl)
res<-snapclust(GI_sub_E,cl,pop.ini=grp$grp);		compoplot(res)
}
dev.off()
pdf(file="QhaE_ca.pdf")	
ca1 <- dudi.coa(tab(GP_sub_E),scannf=FALSE,nf=3); 	s.label(ca1$li, sub="CA 1-2",csub=2)
dev.off()

GI_sub_W<-repool_new(GI_sub,West_pops); 	GP_sub_W<-genind2genpop(GI_sub_W)
levels(GI_sub_W@pop)<-levels(GI_sub@pop)[West_pops]
pdf(file="QhaW_dapc2.pdf")
dapc2 <- dapc(GI_sub_W, n.pca=50,n.da=50);	scatter(dapc2)
dev.off()
pdf("QhaW_snapclust.pdf", height=5, width=18)
for (cl in 2:10){
grp<-find.clusters(GI_sub_W,max.n.clust=12,n.pca=50,n.clust=cl)
res<-snapclust(GI_sub_W,cl,pop.ini=grp$grp);		compoplot(res)
}
dev.off()
pdf(file="QhaW_ca.pdf")	
ca1 <- dudi.coa(tab(GP_sub_W),scannf=FALSE,nf=3); 	s.label(ca1$li, sub="CA 1-2",csub=2)
dev.off()

#can zoom in
compoplot(dapc1,show.lab=T,subset=1:50)
		  

#------------------------------------------------------------------------------#
#				SECTION THREE
#------------------------------------------------------------------------------#
							  		  
#################
#	RELATEDNESS	#
#################

library("Demerelate")
Qha_rel<-read.csv("Qha_binned_8loc_REL_min12.csv")
all_rel_res<-matrix(nrow=22,ncol=7)				

#calculations of each of measures
#tried ritlands and morans but these were inconsistent/ uncorrelated with the others and each other so removed
lois.results<-Demerelate(Qha_rel,object=T,value="loiselle")
all_rel_res[,1]<-unlist(lapply(lois.results$Empirical_Relatedness,mean))
wang.results<-Demerelate(Qha_rel,object=T,value="wang")
all_rel_res[,2]<-unlist(lapply(wang.results$Empirical_Relatedness,mean))
lxy.results<-Demerelate(Qha_rel,object=T,value="lxy") 	#Lynch Ritland
all_rel_res[,3]<-unlist(lapply(lxy.results$Empirical_Relatedness,mean))
Sxy.results<-Demerelate(Qha_rel,object=T,value="Sxy") 	#Lynch alleles sharing
all_rel_res[,4]<-unlist(lapply(Sxy.results$Empirical_Relatedness,mean))
rxy.results<-Demerelate(Qha_rel,object=T,value="rxy") 	#Queller, Goodnight
all_rel_res[,5]<-unlist(lapply(rxy.results$Empirical_Relatedness,mean))			  

par(mfrow=c(2,2))
barplot(unlist(lapply(wang.results$Empirical_Relatedness,half_sib))); abline(h=0.25, lty=3, col="red")
barplot(unlist(lapply(lois.results$Empirical_Relatedness,half_sib))); abline(h=0.25, lty=3, col="red")
barplot(unlist(lapply(rxy.results$Empirical_Relatedness,half_sib))); abline(h=0.25, lty=3, col="red")
barplot(unlist(lapply(lxy.results$Empirical_Relatedness,half_sib))); abline(h=0.25, lty=3, col="red")

}	
#TO DO MAKE SURE NOTE RELATED TO SAMPLE SIZE

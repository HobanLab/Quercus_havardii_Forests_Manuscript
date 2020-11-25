setwd("C:/Users/shoban/Dropbox/havardii/new")

test_dat<-read.csv("binned_dat.csv")


#Only to use the following if need to make single alleles into homozygotes
for (i in seq(4,14,by=2)) { test_NA<-is.na(test_dat[,i]); test_dat[test_NA,i]<-test_dat[test_NA,i-1] }

#Binding columns into a new object, new_dat, in order to prep for genepop file creation
	new_dat<-matrix(nrow=nrow(test_dat),ncol=(ncol(test_dat)-1)/2)
	for (i in 1:5)	new_dat[,i]<-paste(test_dat[,i*2],test_dat[,i*2+1],sep="")
	new_dat<-cbind(as.character(test_dat[,1]),",",new_dat)
#slick way to order on column 1
	new_dat<-new_dat[order(new_dat[,1]),]

#remove individuals with more than one locus missing data
	new_dat<-new_dat[-which(rowSums(new_dat[,]=="-1-1")>1),]
#subsitute in zeroes for -1 for missing data
	new_dat[,2:7]<-gsub("-1-1","000000",new_dat[,2:7])

#Saved as new file which I will then modify by hand
	write.table(new_dat,file="binned_dat.txt",sep="\t",quote=F,row.names=F)

#After this is saved, I opened this file and manually inserted a header, inserted POP before every population 
#This is all that is needed for genepop format, so I saved as binned_dat_form.gen (form means formatted)

#IGNORE FOR NOW looking at pop names to try to remove small pops and insert POP
#pops<-str_sub(new_dat[,1],1,-3)
#table(gsub("-","",pops))
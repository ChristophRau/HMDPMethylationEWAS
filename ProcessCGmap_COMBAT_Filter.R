#This script takes the raw outputs from BSSeeker2, combines the values into a single
#(large) table, removes SNPs, performs COMBAT correction for batch effects
#and then calculates hypervariable (25% shift in at least 5% of strains) CpGs 

library(sva) #for COMBAT


setwd("/media/user/Methylation/Raw_called/")
min_coverage=5  #minimum coverage to be considered


#read it all 'CGmap' files (the output from BSSeeker2)
files=dir()
files=grep(".CGmap",files,fixed=T,value = TRUE)



#Now we start working.  We're generating three files.
#First, a summary file that contains the information on how many 'Good' CpGs exist per sample.
#Second, a growing table of all actual CpGs.  We save these as Meth/Total for ease of use later

#We use the following format to designate
#specific CpGs -   Chromosome_C/G_BasePairLocation
summary_info=c()
splat <- read.table(gzfile(files[1]), header=F)   #Read in the FIRST file
splat=splat[splat$V4=="CG",]   #Filter for only CpG methlyation 
splat=splat[splat$V8>=min_coverage,]  #Filter for samples that have good coverage
summary_info=rbind(summary_info,c(files[1],nrow(splat))) #Add information on # of CpGs per Sample
info_col=paste(splat$V1,splat$V2,splat$V3,sep="_")  #Create the CpG identifier

tomerge=cbind(info_col,paste0(as.numeric(splat$V7),"/",as.numeric(splat$V8)))  #Create the first values for the big table
colnames(tomerge)=c("site",strsplit(files[1],".",fixed = TRUE)[[1]][1]) #and name it.

all_data=tomerge   #Prepare for iteration...

for(i in 2:length(files)){  #for every other file, we will do the same as above, except where noted
  print(i*100/length(files))
  print(files[i])
  splat <- read.table(gzfile(files[i]), header=F)
  splat=splat[splat$V4=="CG",]
  splat=splat[splat$V8>=min_coverage,]
  summary_info=rbind(summary_info,c(files[i],nrow(splat)))
  
  info_col=paste(splat$V1,splat$V2,splat$V3,sep="_")
  tomerge=cbind(info_col,paste0(as.numeric(splat$V7),"/",as.numeric(splat$V8)))
  colnames(tomerge)=c("site",strsplit(files[i],".",fixed = TRUE)[[1]][1])
  all_data=merge(all_data,tomerge,by.x=1,by.y=1,all=T) #Right here, we use the merge
            #function to combine the old data with the new data, inserting NAs for any missing values
}



#We now filter for any SNPs present in the data.
CpG_Names=all_data[,1]
setwd("/media/user/Methylation/Raw_called/SNPs")
All_SNPs=read.csv("All_SNPs.csv")
SNP_Names=All_SNPs[,1]
temp=match(CpG_Names,SNP_Names)
temp=temp[!is.na(temp)]
Relevant_SNPs=All_SNPs[temp,]

for(i in 1:nrow(Relevant_SNPs)){  #goes through each SNP, if it finds it in the data, replaces it with NAs
  print(i/nrow(Relevant_SNPs)*100)
  all_data_row=match(Relevant_SNPs[i,1],all_data[,1])
  all_data_cols=match(colnames(Relevant_SNPs)[which(!is.na(Relevant_SNPs[i,]))[-1]],colnames(all_data))
  all_data[all_data_row,all_data_cols]=NA
}


#and write to a file
setwd("/media/user/Methylation/Raw_called/")
write.csv(all_data,file="results/Merged_CpG_Counts_Cov5.csv",row.names=F)


#this function turns our Meth/Total values into percentages
calc_percent=function(x){  #runs across columns
  temp=c()
  for(i in 1:length(x)){
    vals=as.numeric(strsplit(x[i],"/")[[1]])
    temp=c(temp,vals[1]/vals[2])
  }
  return(temp)
}
#this function gets coverage values for each entry
get_coverage=function(x){ 
  temp=c()
  for(i in 1:length(x)){
    vals=as.numeric(strsplit(x[i],"/")[[1]])
    temp=c(temp,vals[2])
}
return(temp)
}
#this function gets methylated CpG values for each entry
get_Meth=function(x){ 
  temp=c()
  for(i in 1:length(x)){
    vals=as.numeric(strsplit(x[i],"/")[[1]])
    temp=c(temp,vals[1])
  }
  return(temp)
}



#Now we run COMBAT on our data.  Combat only runs on continuous data, so we use percentiles
rownames(all_data)=all_data[,1]
all_data=all_data[,-1]   # get the names into the rownames

all_percents=apply(all_data,2,calc_percent)
all_covs=apply(all_data,2,get_coverage)

Sample_Info=read.csv("../../SampleInfo.csv")   #The sample info for the RRBS libraries

#We found a batch effect based on the date of sequencing

Batches=Sample_Info$Batches

out_percent=ComBat(all_percents,Batches)   #This will take a LONG time.
#We create a version that is in our original Meth/Total version
out_Meth=round(out_percent*all_covs,0)
out_Meth=paste(out_Meth,all_covs,sep="/")

write.csv(out_Meth,"results/Merged_CpG_Counts_COMBAT_Cov5.csv")
write.csv(out_percent,"results/Merged_CpG_Counts_COMBAT_Cov5_Percents.csv")

#Now we filter for sites expressed in at least 70% of samples

present=apply(!is.na(out_percent),1,sum,na.rm=T)/ncol(out_percent)
tokeep=present>.7
out_percent=out_percent[tokeep,]
out_Meth=out_Meth[tokeep,]

write.csv(out_Meth,"results/Merged_CpG_Counts_COMBAT_Cov5_70.csv")
write.csv(out_percent,"results/Merged_CpG_Counts_COMBAT_Cov5_70_Percents.csv")

#Now we filter for hypervariable sites

means=apply(AllPercents,1,mean,na.rm=TRUE)
Differences=AllPercents-means
temp=apply(abs(Differences)>.25,1,sum,na.rm=TRUE)  # number of samples with at least 25% shift from the mean value
temp=temp>8  #Which sites are actually hypervariable
HypVar_percent=out_percent[temp,]
HypVar_Counts=out_Meth[temp,]

#Write to disk
write.csv(HypVar_Counts,"results/HypVar_COMBAT_Cov5_70.csv")
write.csv(HypVar_percent,"results/HypVar_COMBAT_Cov5_70_Percents.csv")

#separate to Meth and Total for MACAU
HypVarMeth=apply(HypVar_Counts,2,get_Meth)
HypVarTotal=apply(HypVar_Counts,2,get_coverage)

write.csv(HypVarMeth,"results/HypVar_Meth.csv")
write.csv(HypVarTotal,"results/HypVar_Total.csv")


#



              
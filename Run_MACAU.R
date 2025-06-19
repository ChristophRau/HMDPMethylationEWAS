#Code to take hypervariable sites and run MACAU
#Hyp_Meth -  m x n matrix of methylated Cs where m are hypervariable CpGs and n are strains 
#Hyp_Total - m x n matrix of all methylated and unmethylated Cs where m are hypervariable CpGs and n are strains

#Note:  Our Cpg names are in the format of "chr#_C_XXXXX" where XXXXX is the base-pair location

#pheno_raw - n x q matrix of phenotype data where n are strains and q are phenotypes



library(doParallel)  #for running MACAU in parallel (saves significant time)
library(qqman)  #for qq and manhattan plots
#Also requires installation of the pyLMM python module for kinship matrix calculation


#Read in Hypervariable Counts Data
setwd("/media/user/Methylation/MACAU/MACAU/EWAS/")
Hyp_Meth=read.csv("HypVar_Meth.csv",row.names=1,check.names=F)
Hyp_Total=read.csv("HypVar_Total.csv",row.names=1,check.names=F)

Hyp_Meth=Hyp_Meth[,order(colnames(Hyp_Meth))]
Hyp_Total=Hyp_Total[,order(colnames(Hyp_Total))]

Ctrls=grep("Ctrl",colnames(Hyp_Meth))  #extract Control/ISO values
Isos=grep("ISO",colnames(Hyp_Meth))
Hyp_Meth_Ctrl=Hyp_Meth[,Ctrls]
Hyp_Total_Ctrl=Hyp_Total[,Ctrls]
Hyp_Total_ISO=Hyp_Total[,Isos]
Hyp_Meth_ISO=Hyp_Meth[,Isos]

#Clean strain identifiers 

colnames(Hyp_Meth_Ctrl)=gsub("_Ctrl","",colnames(Hyp_Meth_Ctrl)) 
colnames(Hyp_Meth_ISO)=gsub("_ISO","",colnames(Hyp_Meth_ISO))
colnames(Hyp_Total_Ctrl)=gsub("_Ctrl","",colnames(Hyp_Total_Ctrl))
colnames(Hyp_Total_ISO)=gsub("_ISO","",colnames(Hyp_Total_ISO))


#Write for future use
write.table(cbind(Hypvar_names,Hyp_Meth_Ctrl),"Hyp_Ctrl_Meth.txt",row.names=F,sep="\t",quote=F)
write.table(cbind(Hypvar_names,Hyp_Meth_ISO),"Hyp_ISO_Meth.txt",row.names=F,sep="\t",quote=F)
write.table(cbind(Hypvar_names,Hyp_Total_Ctrl),"Hyp_Ctrl_Total.txt",row.names=F,sep="\t",quote=F)
write.table(cbind(Hypvar_names,Hyp_Total_ISO),"Hyp_ISO_Total.txt",row.names=F,sep="\t",quote=F)

#Prepare to Run MACAU


#Create the Kinship Matrix

#######################################INSERT CODE HERE, make it CTRL/ISO specific
write.table(Hyp_Meth_Ctrl/Hyp_Meth_Total,file="HypVar_Percentage_Ctrl.txt",row.names=F,sep="\t",col.names = F,quote=F)
write.table(Hyp_Meth_ISO/Hyp_Meth_ISO,file="HypVar_Percentage_ISO.txt",row.names=F,sep="\t",col.names = F,quote=F)

command=paste0("python pylmmKinship.py --SNPemma HypVar_Percentage_Ctrl.txt --NumSNPsemma ",
               nrow(Hyp_Meth_Ctrl)," Ctrl_Meth.k")
system(command)

command=paste0("python pylmmKinship.py --SNPemma HypVar_Percentage_ISO.txt --NumSNPsemma ",
               nrow(Hyp_Meth_ISO)," ISO_Meth.k")
system(command)


#The example code is for the control data, but the ISO data is the same with some name changes

Meth=read.delim("Hyp_Ctrl_Meth.txt",check.names=F,row.names=1)
Total=read.delim("Hyp_Ctrl_Total.txt",check.names=F,row.names=1)

#Read in Phenotypes
pheno_raw=read.csv("Final ISO Project Data by Strain.csv",row.names=1)

pheno=pheno_raw[match(colnames(Meth),rownames(pheno_raw)),] #ensure order matches EWAS, remove any that don't line up

write.csv(colnames(pheno),"Phenotype_Names.csv")

#MACAU has a bug where it doesn't like values that are smaller than 1.
#So to compensate, we scale all phenotypes so the minimum value is 1.

#Run either of the two scripts below:

#This code writes the phenotypes straight (no scaling)
for(i in 1:ncol(pheno)){ 
  this_col=matrix(pheno[,i],ncol=1)
  write.table(this_col,paste0(i,".in"),row.names=F,col.names=F,quote=F,sep=" ")
}   

#This version scales 
for(i in 1:ncol(pheno)){
  pheno[which(pheno[,i]==0),i]=NA
this_col=matrix(pheno[,i],ncol=1)/min(pheno[,i],na.rm=T) #Unfortunate that this is necessary!
write.table(this_col,paste0(i,".in"),row.names=F,col.names=F,quote=F,sep=" ")
}



#And now we run MACAU

cl <- makeCluster(15)
registerDoParallel(cl)

run_script    <-   function(val){  #The very simple script that will run MACAU
  command=paste0("../../macau -g Hyp_Ctrl_Meth_NoCorrect.txt -t Hyp_Ctrl_Total_NoCorrect.txt -p ",val,
                 ".in -k Ctrl_Meth.k -bmm -o ",val)
  system(command)
}


list=c(1:ncol((pheno)))
foreach(i=list) %dopar% run_script(i)



#After this is run, we combine everything back into a single file for analysis 

#getting file names
setwd("/media/user/Methylation/MACAU/MACAU/EWAS/Control_Covar/output")
toread=grep("assoc",dir(),value = T)

#getting phenotype names and ensuring they are in the right order after output
pheno_names=colnames(pheno)
map_names=pheno_names[as.numeric(gsub(".assoc.txt","",toread,fixed=T))]  #the actual phenotype names in the ORDER of the other ones.


CpGNames=rownames(Meth)   #get the names of the CpGs ready 

#This grabs the pvalues for each phenotype and puts it all in one table
outdata=matrix(NA,nrow=length(CpGNames),ncol=length(map_names))
for(i in 1:length(toread)){
  print(i/length(toread)*100)
  this_data=read.delim(toread[i])
  outdata[,i]=this_data$pvalue[match(CpGNames,this_data$id)]
}
#add phenotype names
colnames(outdata)=map_names

#This extracts the chromosome and basepair location for each CpG
CpGInfo=sapply(CpGNames,function(x){strsplit(x,"_")[[1]][c(1,3)]})
CpGInfo=matrix(CpGInfo,ncol=2,byrow = T)
CpGInfo[,1]=gsub("chr","",CpGInfo[,1])
CpGInfo=cbind(CpGInfo,CpGNames)
colnames(CpGInfo)=c("CHR","BP","NAME")

#And then we add it to the pvalues
outdata=cbind(CpGInfo,outdata)

#And save it.
write.csv(outdata,file="Combined_Pvalues_Ctrl_MACAU.csv",row.names=F)



#Make QQs and Manhattans
for(i in 4:ncol(outdata)){  #for each phenotype
  print(i*100/ncol(outdata))   #progress bar
  
  cur_data=as.data.frame(outdata[,c(1:3,i)])  #subset data, remove missing values
  cur_data=cur_data[!is.na(cur_data[,4]),]
  cur_data=cur_data[cur_data[,4]!="NaN",]
  
  #convert X and Y chromosomes to numeric
  cur_data$CHR[cur_data$CHR=="X"]=20
  cur_data$CHR[cur_data$CHR=="Y"]=21
  cur_data$CHR=as.numeric(cur_data$CHR)
  cur_data$BP=as.numeric(cur_data$BP)
  
  #Grab the phenotype name, do final checks
  name=colnames(cur_data)[4]
  colnames(cur_data)[4]="P"
  cur_data$P=as.numeric(cur_data$P)
  
  #create Manhattans and QQ plots
  png(paste0(name,".Manhattan.png"))
  manhattan(cur_data,genomewideline = -log10(4.2e-6),chrlabs = c(1:19,"X","Y"),main=name)
  dev.off()
  png(paste0(name,".qq.png"))
  qq(cur_data$P,main=name)
  dev.off()
}



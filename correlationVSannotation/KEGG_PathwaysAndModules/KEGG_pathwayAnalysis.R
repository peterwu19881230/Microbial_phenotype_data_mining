# Use KEGG pathway to do the correlationVSannotation experiment that I did for EcoCyc pathways
#(!) The definition for pathways in KEGG is too broad and over-encompassing, so this appraoch doesn't work very well

KEGGpwy=read.table("http://rest.kegg.jp/link/eco/pathway",stringsAsFactors = F)
colnames(KEGGpwy)=c("Pathway","ACCESSION-1 - ")
KEGGpwy$"ACCESSION-1 - "=sub("eco:","",KEGGpwy$"ACCESSION-1 - ",fixed = T) #remove "eco:" at the front


##I am trying to remove some pathways that are too broad (ex. path:eco01100: Metabolic pathways):
KEGGpwy=KEGGpwy[!KEGGpwy$"Pathway" %in% c("path:eco01100",
                                         "path:eco01110",
                                         "path:eco01120",
                                         "path:eco02010",
                                         "path:eco02020",
                                         "path:eco01230",
                                         "path:eco01200"),]



##Get bxxxx accession No. from genes.dat
name_ACC1_ACC2=sapply(seq_along(genes.dat),FUN=function(i){
  c(names(genes.dat[i]),genes.dat[[i]]$'ACCESSION-1 - ',genes.dat[[i]]$'ACCESSION-2 - ')
})
name_ACC1_ACC2=Reduce(rbind,name_ACC1_ACC2)
colnames(name_ACC1_ACC2)=c("EcoCycID","ACCESSION-1 - ","ACCESSION-2 - ")


##use bxxxx in name_ACC1_ACC2 to associate ECK_1st_table to KEGGpwy => To get Nichols' ID
#(!)Not sure these 2 merging steps are what I want. Also, even if this is not good enough, I shouldn't have gotton the result I have
temp_table=merge(ECK_1st_table,name_ACC1_ACC2,by="EcoCycID") 
temp_table=merge(temp_table,KEGGpwy,by="ACCESSION-1 - ")

NicholsID_bxxxx=temp_table[,c("ids","ACCESSION-1 - ")]
KEGGpwy=merge(KEGGpwy,NicholsID_bxxxx,by="ACCESSION-1 - ")
KEGGpwy$ids=as.character(KEGGpwy$ids)  

##Fix the duplicated rows
library(dplyr)
KEGGpwy=distinct(KEGGpwy)

attribute_list=attr_list(KEGGpwy$ids,KEGGpwy$Pathway)


listOfResults=dist_TF_cumsum_old(All_Data_NAimputed[unique(KEGGpwy$ids),],attribute_list)


#Have to finish the MHD part...



##Plot them
slopeForControl=sum(listOfResults[[1]]$sameORnot)/length(listOfResults[[1]]$sameORnot)
samples=seq(from=0,to=681528,by=73); samples[1]=1
library(ggplot2)
library(scales) #This is for label=comma (not using scientific notation and put commas) to work
ggplot()+geom_line(data=data.frame(no=samples,CumSum=listOfResults[[1]]$cumsum[samples]),aes(no,CumSum,color="PCC"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[2]]$cumsum[samples]),aes(no,CumSum,color="PCC_squared"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[3]]$cumsum[samples]),aes(no,CumSum,color="euclidean"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[4]]$cumsum[samples]),aes(no,CumSum,color="maximum"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[5]]$cumsum[samples]),aes(no,CumSum,color="manhattan"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[6]]$cumsum[samples]),aes(no,CumSum,color="canberra"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[7]]$cumsum[samples]),aes(no,CumSum,color="binary"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[8]]$cumsum[samples]),aes(no,CumSum,color="minkowski"))+
  geom_line(data=data.frame(no=samples,CumSum=listOfResults[[9]]$cumsum[samples]),aes(no,CumSum,color="spearman"))+  
  geom_abline(intercept=0,slope=slopeForControl,linetype='dotted',color="black")+ #negative control
  geom_abline(intercept=0, slope=1,linetype='dotted',color="black")+ #positive control
  scale_y_continuous(expand = c(0,0),label=comma)+scale_x_continuous(expand = c(0,0),label=comma)+
  labs(x = "Low--Distance--High",y="No. of being in the same pwy/ptcomplex",aesthetic='custom text')+
  scale_color_discrete(name="Distance")








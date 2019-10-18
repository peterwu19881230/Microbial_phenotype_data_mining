#CGSC data
currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)

CGSC=read.table(paste(currentDir,"CGSC/JWstrainInfoFromCGSC.txt",sep="/"),stringsAsFactors = F,quote="",sep="\t")[,c(1,2,4)]
names(CGSC)=c("strain","JWsyn","url")



##======8/6/2019 gsheet package doesn't work anymore======

##listOfPackages=c("googlesheets", "tidyverse")
##new_pack=listOfPackages[!(listOfPackages %in% installed.packages()[,"Package"])]
##if(length(new_pack)) install.packages(new_pack)

##library(googlesheets); library(tidyverse)

##gs_ls() #This will open a browser and need user's authorization. .httr_oauth is stored in the same directory so the next time this question won't be asked again 
##sheet=gs_title("Nichols_FinalListOfStrainsForMakingPagesInOMPwiki") #Read the sheet identified by the name
##gs_ws_ls(sheet)

  ##Try to add some stuff (Just for practice. Not necessary because there are only 2 different strings for the new column):
  #Strain_bkgd=gs_read_cellfeed(sheet, ws = "2018-07-14_ReadyToMakeStrainPages", range = "E1:E3956") %>% gs_simplify_cellfeed()
  #new_Strain_bkgd=paste("parent: ",Strain_bkgd,sep="") 
  ##If I don't use fixed(), the content will be interpreted as regular expressions. Can read fixed() from R Documentation
  #dataToPut=gs_edit_cells(sheet, ws = "2018-07-14_ReadyToMakeStrainPages",input = c("Strain bkgd",new_Strain_bkgd), anchor = "E1")

##fromDebby=gs_read(sheet,ws="2018-07-14_ReadyToMakeStrainPages")
##========================================================


filePath=paste(currentDir,"2018-07-14_ReadyToMakeStrainPages.tsv",sep="/")
fromDebby=read.table(filePath,sep="\t",fill=T,quote="",header=T)
names(fromDebby)=c("ECK","strain","gene","allele","background","category","reference") #Note: some of them don't have JW number 



##Double check if anything is duplicated
sum(duplicated(CGSC$strain)) #0

fromDebby$strain[duplicated(fromDebby$strain)] #Except this, nothing is duplicated: "BW25113 + allele in column D"

mergedData=merge(x=fromDebby,y=CGSC,all.x=T,all.y=F,by="strain")

##replace "" with NA so the file looks cleaner. ("" will cause php to complain "Undefined offset" when parsing the resulting text file)
mergedData$category[mergedData$category==""]=NA


##Check what strains that don't have JW numbers look like
mergedData$strain[!grepl("JW[0-9]{4}",mergedData$strain)]



#write.table(mergedData,"strainInfo.txt",row.names = F,col.names = F,quote=F,sep="\t")

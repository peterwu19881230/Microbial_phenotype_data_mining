#Examine the fitness scores under different treatment more closely
#uniqueChemIndex is used (It is defined in obj_prep.R)


##I am curious in determining how many treatments have conflicting fitness scores (eg. one is -4 and the other is 4). 
#Or maybe those shouldn't be called conflicting scores because different concentrations can have different effect
system.time(
  
  inconsistencyNo<-sapply(uniqueChemIndex,FUN=function(indices){
    
    if(length(indices)==1){return(rep(0,3979))} #treatment with 1 concentration won't have inconsistency
    
    treatmentData=Ternary_Data[,indices]
    
    inconsistentTr=apply(treatmentData,1,FUN=function(row){
      ifelse((unique(row) %>% length)==1,0,1)
    })
    
    return(inconsistentTr)
    
  })
  
  
)

##Total number of inconsistency for the whole ternary data
sapply(inconsistencyNo,sum) %>% sum #7422

##Percentage across all scores(doesn't seem to be too many inconsistent phenotypes)
(sapply(inconsistencyNo,sum) %>% sum)/(114*3979)*100 #1.636222



##Categorization of scores under 114 unique treatments: For each treatment, highlight the knockouts that have at least 1 significant phenotype => upload to a google sheet

if(!exists("tr_sigPheno")){
  
  dat=Ternary_Data_324cutff  #Note: before I used this as the input: Ternary_Data_NAnotimputed 
  geneNames=read.csv("Data/allData.csv",stringsAsFactors = F)[,1] #Get the original names of knockouts from Nichols' original file 
  colnames(dat)=colnames(read.csv("Data/allData.csv",stringsAsFactors = F,check.names = F))[-1] # A temporary fix for the colnames
  
  tr_sigPheno=list()
  count=1
  for(i in 1:length(uniqueChemIndex)){
    
    complete_mat=dat[,uniqueChemIndex[[i]]]
    
    if(is.null(dim(complete_mat)[2])) next  #if there is only 1 concentration for that treatment, skip and run the next i 
    
    index=apply(complete_mat,1,FUN=function(row){
      ifelse( (1 %in% row) | (-1 %in% row) ,T,F)
    })
    
    out=complete_mat[index,]
    
    tr_sigPheno[[count]]=out; names(tr_sigPheno)[count]=names(uniqueChemIndex)[count]
    count=count+1
  }
  save(tr_sigPheno,file="Data/sourced/tr_sigPheno.RData")
  
}

##Verify that the above for() works correctly:
##temp=tr_sigPheno[[1]]
##temp2=tr_sigPheno[[2]]
##temp3=tr_sigPheno[[3]]




##Upload tr_sigPheno. Each element of the list is uploaded as a sheet. The sheet name is the treatment
gs_ls() #This will open a browser and needs user's authorization 
file_name="significant_phenotypes_under_114_unique_treatments"
gObj_tr_sigPheno=gs_new(file_name) #Create a new spread sheet on google drive:


for(i in 1:20){  #Should be 1:length(tr_sigPheno), but it exceeds the maximum No. See explainations below this for()
  gs_ws_new(gObj_tr_sigPheno, ws_title = names(tr_sigPheno[i])) #Add a worksheet inside the google spread sheet
  gObj_tr_sigPheno=gs_title(file_name) #Have to re-register the title after the worksheet is created, otherwise the worksheet won't be recognized (Ref: https://stackoverflow.com/questions/47178305/googlesheets-r-not-recognising-worksheet-number-2-to-add-rows?rq=1)
  strainNameForTr=geneNames[rownames(tr_sigPheno[[i]]) %>% as.numeric]
  gs_edit_cells(gObj_tr_sigPheno,ws=names(tr_sigPheno[i]), input = cbind("Original name in Nichols'"=strainNameForTr,tr_sigPheno[[i]]),col_names = T)
}

gs_ws_delete(gObj_tr_sigPheno,ws="Sheet1") #Delete the default worksheet after the above is done

##If I create different treatment in separate tabs the maximum tabs is 76. After that this error will appear: Error in function_list[[k]](value) : Bad Request (HTTP 400).
##How did I confirm this? : 1. I created a sheet of 75 tabs => got the error  2. Deleted 1 tab => Can make 1 additional tab, then the error re-occured

##Move the new sheet inside this folder on google drive: OMP/NicholsDataMining_Data/Data
##This function is from "googledrive" package which was created by the same author that created the "googlesheet" package -> Jennifer Bryan

drive_mv(file=file_name,path="OMP/NicholsDataMining_Data/Data/")





##Categorization of scores under 114 unique treatments: For each treatment, highlight the knockouts that have at least one ternary fitness score=1 and one ternary fitness score=-1 => upload to a google sheet
##(updated 7/24/2018) quantitative scores are included underneath ternary scores

if(!exists("tr_conflictPheno")){
  
dat=Ternary_Data_324cutff  #Note: before I used this as the input: Ternary_Data_NAnotimputed 
colnames(dat)=colnames(read.csv("Data/allData.csv",stringsAsFactors = F,check.names = F))[-1] # A temporary fix for the colnames

dat_quant=All_Data
colnames(dat_quant)=colnames(dat) # A temporary fix for the colnames

geneNames=read.csv("Data/allData.csv",stringsAsFactors = F)[,1] #Get the original names of knockouts from Nichols' original file 

  

tr_conflictPheno=list()
count=1
for(i in 1:length(uniqueChemIndex)){
  
  if(length(uniqueChemIndex[[i]])<3) next 
  #If No. of concentrations < 3 it's easy to explain any score combination
  
  complete_mat=dat[,uniqueChemIndex[[i]]] %>% as.data.frame # as.data.frame is to prevent an 1-row complete_mat from being converted to vecotr and loss the rowname
  complete_mat_quant=dat_quant[,uniqueChemIndex[[i]]] %>% as.data.frame
  
  index=apply(complete_mat,1,FUN=function(row){
    ifelse( (1 %in% row) & (-1 %in% row) ,T,F)
  })
  
  if(sum(index)==0) next #If there is no row found, skip to the next iteration
  
  out=complete_mat[index,]
  out_quant=complete_mat_quant[index,]
  
  out_combined=rbind(out,out_quant)
  
  
  tr_conflictPheno[[count]]=out_combined; names(tr_conflictPheno)[count]=names(uniqueChemIndex)[i]
  count=count+1
}
save(tr_conflictPheno,file="Data/sourced/tr_conflictPheno.RData")

}


##Upload tr_conflictPheno
gs_ls() #This will open a browser and needs user's authorization 
file_name="possible_conflicts_phenotypes_under_114_unique_treatments"
gObj_tr_conflictPheno=gs_new(file_name) #Create a new spread sheet on google drive:
gs_ws_rename(gObj_tr_conflictPheno,"Sheet1","possibleConflicts") #Rename the 1st tab
gObj_tr_conflictPheno=gs_title(file_name) #Have to re-register the title after the worksheet is created (or changed name), otherwise the worksheet won't be recognized (Ref: https://stackoverflow.com/questions/47178305/googlesheets-r-not-recognising-worksheet-number-2-to-add-rows?rq=1)

bound=1 #Initialize value for bound cell
for(i in 1:length(tr_conflictPheno)){  
  strainNameForTr=geneNames[rownames(tr_conflictPheno[[i]])[1: (dim(tr_conflictPheno[[i]])[1]/2) ] %>% as.numeric]
  
  gs_edit_cells(gObj_tr_conflictPheno,ws="possibleConflicts", input = cbind("Original name in Nichols'"=strainNameForTr,tr_conflictPheno[[i]]),col_names = T,
                anchor =paste("R",bound,"C1",sep=""))
  
  bound=bound+dim(tr_conflictPheno[[i]])[1]+1 # +1 because of colnames
}


##Move the new sheet inside this folder on google drive: OMP/NicholsDataMining_Data/Data
##This function is from "googledrive" package which was created by the same author that created the "googlesheet" package -> Jennifer Bryan

drive_mv(file=file_name,path="OMP/NicholsDataMining_Data/Data/")


#How many total rows were extracted in tr_conflictPheno (not including the corresponding quantitative scores)?
(sapply(tr_conflictPheno,dim)[1,] %>% sum)/2 #52





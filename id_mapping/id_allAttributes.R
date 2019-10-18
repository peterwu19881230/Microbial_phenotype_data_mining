#Create a master object for Nichols' ids


#The following objs are used
##ECK_1st_table #name and strain availability
apply(ECK_1st_table,2,class)
class(ECK_1st_table)


##GO term
load("Data/id.GO") #defined in 2nd_Table_in_SQL.R
used_id.GO=id.GO[,-3]
names(used_id.GO)[2]="GO"
apply(used_id.GO,2,class)
class(used_id.GO)


##UniProt ID 
load("Data/id.UniProt.RData") #defined in 2nd_Table_in_SQL.R
used_id.UniProt=id.UniProt[,-3]
names(used_id.UniProt)[2]="UniProtID"
apply(used_id.UniProt,2,class)
class(used_id.UniProt)


##Pathway 
load("Data/ECK.Pathway.RData")
used_ECK.Pathway=ECK.Pathway[,-3] #defined in 2nd_Table_in_SQL.R
names(used_ECK.Pathway)[2]="Pwy"
apply(used_ECK.Pathway,2,class)
class(used_ECK.Pathway)


#protein complex (from proteinComplex.R) 
load("Data/ids.pcomplex.RData")
apply(ids.pcomplex,2,class)
class(ids.pcomplex)


#Regulon ECK12 id (from parse_regulonDB.R)
load("Data/used_id.regulatorID.RData")


  
#Operon ECK12 id (from parse_regulonDB.R)
load("Data/used_id.operonID.RData")




#kegg modules (from scrapeKEGG.R)
used_id.KEGGmodules=data.frame()
for(i in 1:length(KEGGmodulesForNichols)){
  id=names(KEGGmodulesForNichols)[i]
  
  for(j in 1:length(KEGGmodulesForNichols[[i]])){
    if(is.na(KEGGmodulesForNichols[[i]][j])) break
    used_id.KEGGmodules=rbind(used_id.KEGGmodules,c(id,KEGGmodulesForNichols[[i]][j]))
  }
}
colnames(used_id.KEGGmodules)=c("ids","kegg_modules")
str(used_id.KEGGmodules)


#Map the ids from this paper: Genomewide landscape of gene–metabolome associations in Escherichia coli

# A column that tells whether the strain has any significant phenotypes (added 3/21/19)
TF=apply(Ternary_Data_324cutff_NAremoved,1,FUN = function(row){
  sum(row!=0)>0  
})
id.TF=data.frame(ids=as.character(1:3979), AnySigFitness=TF)


#Merge them altogether (!)Note: for the protein complex, the unwanted annotations for 
id_allAttributes=list(ECK_1st_table,used_id.GO,used_id.UniProt,used_ECK.Pathway,ids.pcomplex,used_id.regulatorID,used_id.operonID,used_id.KEGGmodules,id.TF) %>%
  Reduce(function(df1,df2) left_join(df1,df2,by="ids"),.)

#Check if it is correct before saving 
apply(id_allAttributes,2,FUN=function(attr){
  return(c(class(attr),attr %>% unique %>% length)) #(!) class(attr) doesn't necessary return the right type and I don't know why
})

str(id_allAttributes)


#Add common names that are searchable on regulonDB website (ECK12 identifiers don't find stuff when searching on regulonDB)

#Regulators' name
generegulation_tmp=read.table("Data/regulon/generegulation_tmp.txt",sep="\t",quote="",stringsAsFactors = F) #no. of rows match the file: 4986-37+1=4950
ECK12_regulator=generegulation_tmp[,c(1,4)] %>% unique; names(ECK12_regulator)=c("regulator","regulator_name")
##Note: ECK12_regulator is a many-to-many-relationship table

ECK12_regulator_duplicated=ECK12_regulator[duplicated(ECK12_regulator$regulator),]



#Operon name
ECK12_operon=read.table("Data/regulon/operon.txt",comment.char = "#",fill=T,sep="\t",quote="",header=F,check.names = F)[,1:2]
names(ECK12_operon)=c("operon","operon_name")


#join the 2 columns: regulator_name, operon_name
temp=left_join(id_allAttributes,ECK12_regulator,by="regulator")
temp=left_join(temp,ECK12_operon,by="operon")


#drag these 2 colums to near their identifier columns: regulator_name, operon_name
temp=temp[,c(1:16,17,21,18,22,19:20)] 

id_allAttributes=temp

#save(id_allAttributes,file="Data/sourced/id_allAttributes.RData")
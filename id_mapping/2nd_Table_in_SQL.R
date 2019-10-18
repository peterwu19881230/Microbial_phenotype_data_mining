#Goal: Create the 2nd table for building the relational database of Nichols' strains - GO - UniProt
#1. Use clean columns from 1st_tabe.R for gene names in Nichols' data 
#2. Map the GO IDs to Nichols' strains by the .GAF file Suzi made (Note that the GO annotations here are different than if retrieved from UniProt)
#3. Map the UniProt IDs to Nichols' strains
#4. Map the pathways based on pathwayscol (a file from EcoCyc) to Nichols' strains
#5. transform the format to id-data-datatype and upload to SQL


#1. Use ECK_1st_table created by id_mapping/1st_table.R 
load("Data/ECK_1st_table.RData")



#Step 2: Map the GO IDs to Nichols' strains by the .GAF file Suzi made

##read.table doesn't work for .txt files that have NA values in some rows. I have to convert the file "2017_05_ECgene_association.ecocyc.txt" into .csv and then use read.csv

##Read the file
ECgene_association.ecocyc=read.csv("Data/2017_05_ECgene_association.ecocyc.csv",stringsAsFactors=F) 
##There are newer versions. Probably doesn't differ much



##Give the column names according to the teamwiki page: https://hexamer.tamu.edu/team/wiki/index.php/Making_GAF
colnames(ECgene_association.ecocyc)<-c("DB","DB Object ID","DB Object Symbol","Qualifier","GO ID","DB:Reference (IDB:Reference)","Evidence Code","With (or) From","Aspect","DB Object Name","DB Object Synonym (ISynonym)","DB Object Type","Taxon(ITaxon)","Date","Assigned By","Annotation Extension","Gene Product Form ID")


##map based on DB Object ID (EcoCyc symbol)
tab1=ECK_1st_table[,c("ids","EcoCycID")]
tab2=ECgene_association.ecocyc[,c("DB Object ID","GO ID")]
temp=inner_join(tab1,tab2,by=c("EcoCycID"="DB Object ID"))

names(temp)[3]="Data"
temp$'Data Type'="GO ID"
temp=temp[,c("ids","Data","Data Type")]
temp=temp[order(temp$ids %>% as.numeric),]
id.GO=temp


save(id.GO,file="Data/id.GO")




#Step 3. Map the UniProt IDs to Nichols' strains


##Get product from genes.dat obj. 
###Some of them doesn't have product (makes sense because not all genes produce proteins)
###Some of them have more than 2 products
EG.product.table=data.frame()
for(i in 1:dim(ECK_1st_table)[1]){ #here i is not Nichols' ID
  
  EcoCycID=ECK_1st_table$EcoCycID[i]
  
  product=paste("ECOCYC:",genes.dat[[EcoCycID]]$`PRODUCT - `,sep="")
  
  if(length(product)>=2) print(product)
  
  
  if(is.null(product)==T){  ##If there is no product
    #print(i); print("=> Doesn't have product") 
    next
  }
  
  
  EG.product.table=rbind(EG.product.table,
  cbind(ECK_1st_table$ids[i],EcoCycID,product), 
  stringsAsFactors=F
  )
  
}

#clean the propagated rows caused by propagated rows in ECK_1st_table (ECK_1st_table has 3986 rows, not 3979)
EG.product.table=unique(EG.product.table)
names(EG.product.table)[1]="ids"

##This shows ECK with more than 2 products
###EG.product.table$ids[duplicated(EG.product.table$ids)]
###ECK.id=EG.product.table$ids[duplicated(EG.product.table$ids)] %>% unique %>% as.numeric
###EG.MoreThan2=EG.product.table[EG.product.table$ids %in% ECK.id,]


write.table(EG.product.table$product,"Data/EcoCycProduct.txt",quote=F,col.names=F,row.names=F)
##Use this file to get UniProt IDs: on http://www.uniprot.org/uploadlists/ (select: from "BioCyc" to "UniProtKB")=> download the text file (I used tab separated format ) and renamed to: uniprot.txt


###Why using quote="" in read.csv?
###https://biowize.wordpress.com/2013/10/08/quotation-marks-and-rs-read-table-function/
UniProt_Mapping=read.table("Data/uniprot.txt",sep="\t",header=T,stringsAsFactors = F,quote="")
colnames(UniProt_Mapping)[1:2]=c("yourlist","Entry")
##Not all the product ID (eg. EcoCyc:G7002-MONOMER) returns an UniProt ID. (Updated 9/11/2019) EcoCyc:G7002-MONOMER now gets P64512

##Not all the input queries find a UniProtID:
is.found=list()
for(i in 1:length(EG.product.table$product)){
  product=EG.product.table$product[i]
  is.found[[i]]=which(UniProt_Mapping$yourlist %in% product)
  if(identical(is.found[[i]],integer(0))==T){ print(paste("id= ",EG.product.table$i[i],",",product," is not found"))}
  }
##The above shows the UniProt IDs not found from the EcoCyc:XXXX. If there is only EcoCyc:(blank), it means the ECK doesn't find the EcoCyc product identifier




##Create the SQL table: id (From ECK)- UniProt ID - Data Type (UniprotID)
id.UniProt=data.frame()
for(i in 1:dim(EG.product.table)[1]){
  
  product=EG.product.table$product[i]  
  index=which(UniProt_Mapping$yourlist %in% product) ##This index is the index for UniProt_Mapping
  
  if(identical(index,integer(0))==F){ #Will there be situs that index has more than 2?
    UniProtID=UniProt_Mapping$Entry[index] ##$Entry is the UniProt ID
    id.UniProt=rbind(id.UniProt,c(EG.product.table[i,1],UniProtID,"UniProt ID"),stringsAsFactors=F)
  }
}
colnames(id.UniProt)=c("ids","Data","Data Type")

###How many genes have 2 protein products? => None (by looking at the ECK - Uniprot ID)
morethan2=list()
count=1
for(i in 1:3979){
  if(length(which(id.UniProt$ids %in% i))>=2){
    morethan2[[count]]=i
    count=count+1
  }
}
morethan2
###Or simply do the following to verify:
###duplicated(id.UniProt$ids) %>% sum  

save(id.UniProt,file="Data/id.UniProt.RData")





#4. Map the pathways based on pathwayscol (a file from EcoCyc) to Nichols' strains

pathwayscol=read.csv("Data/pathwayscol.csv",comment.char="#",stringsAsFactors=F)[-1,c(1:2,112:220)]
names(pathwayscol)=pathwayscol[1,]
pathwayscol=pathwayscol[-1,]
##1074 genes are in at least 1 pathway, others are not


##EcoCyc ID - Pathway ID - Pathway Name
EcoCycID.Pwys=data.frame()
for(i in 1:dim(pathwayscol)[1]){
  
  j=3 #Gene IDs start from the 3rd column
  while(is.na(pathwayscol[i,j])==F){ ##empty cells were read as NAs from read.xlsx()
    EcoCycID=pathwayscol[i,j]
    
    
    EcoCycID.Pwys=rbind(EcoCycID.Pwys,
                        c(EcoCycID,pathwayscol$`UNIQUE-ID`[i],pathwayscol$NAME[i]),
                        stringsAsFactors=F)
    
    j=j+1
    if(j>dim(pathwayscol)[2]) break
  }
  
}
colnames(EcoCycID.Pwys)=c("EcoCyc ID","Pathway ID","Pathway Name")

##Combine EcoCycID.Pwys with ids-ECK

###genes to pathways is a many to many relationship:
###length(unique(EcoCycID.Pwys$`EcoCyc ID`))  ###1074 unique genes
###length(EcoCycID.Pwys$`EcoCyc ID`) ###2681
###length(unique(EcoCycID.Pwys$`Pathway ID`)) ###424 unique pathways
###length(EcoCycID.Pwys$`Pathway ID`) ###2681



ECK.Pathway_table=data.frame()
for(i in 1:3979){
  
  index=which(EcoCycID.Pwys$`EcoCyc ID` %in% ECK_1st_table$EcoCycID[i])
  
  pathwayIDs=EcoCycID.Pwys$`Pathway ID`[index]  ##1 gene might be involved in multiple pathways. 
  ##Some genes might not match a pathway?
  
  if(identical(pathwayIDs,character(0))==T){
    pathwayIDs="Not in any pathway"
  }
  
  
  for(j in 1:length(pathwayIDs)){ 
    ECK.Pathway_table=rbind(ECK.Pathway_table,
                            c(ECK_1st_table$ids[i],
                              ECK_1st_table$`Original Name`[i],
                              ECK_1st_table$sorted_ECK_missing_gene_names_added[i],
                              pathwayIDs[j]),
                            stringsAsFactors=F)
  }
  
}
colnames(ECK.Pathway_table)=c("ids","Original Name","sorted_ECK_missing_gene_names_added","Pathway")

##Create the SQL table: id (From ECK)- Pathway - Data Type (Pathway)
ECK.Pathway=cbind(ECK.Pathway_table[,c(1,4)],Data="Pathway",stringsAsFactors=F)
ECK.Pathway=ECK.Pathway[ECK.Pathway$Data!="Not in any pathway",] #remove things not involved in any pathway


colnames(ECK.Pathway)=c("ids","Data","Data Type")
save(ECK.Pathway,file="Data/ECK.Pathway.RData")

second.table=rbind(id.GO,id.UniProt,ECK.Pathway)




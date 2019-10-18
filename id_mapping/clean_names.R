#Goal: Correct the original ECKs from Nichols'
#Step 1: Sort the name column by ECK. Label duplicated strains with -1, -2

###Note: allData.csv is directly from from Nichols' supplemental
ECKs_rows_fixed=read.csv(file="Data/allData.csv",colClasses=c("character", rep("NULL", 324))) %>% t %>% as.vector ##colClasses is used to skip all the other columns
ECKs_rows_fixed=sort(ECKs_rows_fixed)

dup1=ECKs_rows_fixed[duplicated(ECKs_rows_fixed,fromLast=T)] %>% paste0("-1")
dup2=ECKs_rows_fixed[duplicated(ECKs_rows_fixed,fromLast=F)] %>% paste0("-2")

index_dup_1=which(duplicated(ECKs_rows_fixed,fromLast=T))
index_dup_2=which(duplicated(ECKs_rows_fixed,fromLast=F))

ECKs_rows_fixed[index_dup_1]=dup1
ECKs_rows_fixed[index_dup_2]=dup2




#Step 2: Correct the strain names


#With SPA tag
SPAtag<-ECKs_rows_fixed[grep(pattern="- SPA",x=ECKs_rows_fixed)]
#With DAS tag
DAStag<-ECKs_rows_fixed[grep(pattern="- DAS",x=ECKs_rows_fixed)]
#9 point-mutant alleles of essential genes (plus corresponding ‘‘linked’’ strains for each of the 9 alleles, in which the antibiotic-resistance cassette was linked to the wild-type allele as a control)
mt_linked<-ECKs_rows_fixed[grep(pattern="Linked",x=ECKs_rows_fixed)]
#2 truncation of essential genes
truncations<-ECKs_rows_fixed[grep(pattern="Truncation",x=ECKs_rows_fixed)]


#others
isnot<-which(
  ECKs_rows_fixed %in%
    
    c(
      SPAtag,
      DAStag,
      mt_linked,
      truncations
    ) 
)

other_genes<-ECKs_rows_fixed[-isnot]
rm(isnot)

##In other_genes, I want to cut them into sub-groups

##Most common expression in the whole data
most_common_genes<-other_genes[
  c(
    grep(pattern="^ECK[0-9]{4}-[a-zA-Z]{3}[a-zA-Z0-9]$",x=other_genes)
  )
  ]




##Only ECKXXXX <= This is where we have to add gene names. I don't know why there weren't gene names.
onlyECK<-other_genes[
  grep(pattern="^ECK[0-9][0-9][0-9][0-9]$",x=other_genes)
  ]

##Add gene names for ECK_only manully: from OMP shared/Nichols_phenotypic-landscape/Nichols-Strains-Sorted.xlsx (search for all tabs)
##Some gene names might not be in Nichols-Strains-Sorted.xlsx?
##format: ECKXXXX->ECKXXXX-YYYY

onlyECK_gene_name_added<-c(
  "ECK0012-HTGA", "ECK0017-HOKC", "ECK0266-YKGN", "ECK0320-YAHH", "ECK0359-YAIF", "ECK0367-YKIB", "ECK0369-YAIU", "ECK0503-YBBV",
  "ECK0619-YBEM", "ECK0679-YBFH", "ECK1128-YMFH", "ECK1132-CROE", "ECK1159-YMGG", "ECK1160-YMGH", "ECK1453-YNCM", "ECK1933-YEDM",
  "ECK1990-YOEE", "ECK2132-YOHH", "ECK2331-YCFT", "ECK2636-YPJM", "ECK2637-YPJM", "ECK2647-YPJC", "ECK2650-YGAR", "ECK2651-YGAC",
  "ECK2652-YGAD", "ECK2675-YGAY", "ECK2854-YGEL", "ECK2856-YGEN", "ECK2859-YGEQ", "ECK2994-YGHY", "ECK3474-YHIK", "ECK3672-YSDC",
  "ECK3675-GLVC", "ECK3769-YIFN", "ECK3802-YZCX", "ECK4097-PHNE", "ECK4219-YZFA", "ECK4265-YJGW", "ECK4330-YJIQ", "ECK4334-YJIV",
  "ECK4426-TISA"
)

'
Note:
*No gene name is found for ECK0012 within Nichols-Strains-Sorted.xlsx. However I did find it inside OMP shared/Nichols_phenotypic-landscape/Sorted_ECK_by_Peter/Keio_strains_with_verified_JW. But I fotgot where I got the JW numbers and the matched genes from 
*ECK1132 has gene name synonyms: croE, ymfT (http://ecoliwiki.net/colipedia/index.php/croE:Gene)
*Ecoliwiki indicates that ECK2636 and ECK2637 are the same 
'

##remaining
isnot<-which(
  
  other_genes %in%
    c(
      most_common_genes,
      onlyECK
    ) 
)

remaining_genes<-other_genes[-isnot]

rm(isnot)

##check the total length (should be 3979):
##length(SPAtag)+length(DAStag)+length(mt_linked)+length(truncations)+length(most_common_genes)+length(onlyECK)+length(remaining_genes)


##Start making columns

##Sorted by category
sorted_ECK<-c(SPAtag,DAStag,mt_linked,truncations,most_common_genes,onlyECK,remaining_genes)
sorted_ECK_missing_gene_names_added<-c(SPAtag,DAStag,mt_linked,truncations,most_common_genes,onlyECK_gene_name_added,remaining_genes)

##ECK only
library(stringr)

ECK_only<-str_extract(sorted_ECK_missing_gene_names_added,"^ECK[0-9][0-9][0-9][0-9]")


##gene name only (This is based on the assumption that all the full names start with ECKXXXX-)
ECKXXXXremoved<-str_extract(sorted_ECK_missing_gene_names_added,"-(.*)") 


associated_gene_names<-str_extract(ECKXXXXremoved,"[A-Z0-9]{3,}") 
###These are just gene names, not genotypes (Genes that are associated with the genotypes. Some genotypes are actually more like wild-type )
###This is a little dangerous since whether the names are actual gene names cannot be verified. Will verify after mapping to gene names in the .GAF file

#manualy correct istR -> istR-1 (this is a phantom gene)
associated_gene_names[3977]="ISTR-1"


rm(ECKXXXXremoved)

##bind them 
ECK_name_columns<-cbind(sorted_ECK,sorted_ECK_missing_gene_names_added,ECK_only,associated_gene_names)


##Replace names for exceptions (The rules above don't apply)
###ECK0005-TP2
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK5005-TP2","associated_gene_names"]<-"TP2"

###ECK0086-A-Linked: This is just the wild-type E.coli that has the antibiotic resistance gene linked to MURE
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK0086-A - Linked","associated_gene_names"]<-"MURE"

###ECK0086-B-Allele: This is the actual mutant (the mutant is called MurE1)
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK0086-B - Allele","associated_gene_names"]<-"MURE1"

###ECK0086-C-SPA: This is also a kind of wt that has SPA+KAN after the orf of interest
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK0086-C - SPA","associated_gene_names"]<-"MURE"


###ECK0055-D-IMP4213 - Allele
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK0055-D-IMP4213 - Allele","associated_gene_names"]<-"imp"

###ECK4142-4143-ECNAB
ECK_name_columns[ECK_name_columns[,"sorted_ECK"]=="ECK4142-4143-ECNAB","associated_gene_names"]<-"ecnA,ecnB"


#Step2: Corrected strain names (Labeled duplicated strains with -1, -2. Added missing ECK names)


#Original strain names (created with an id column in which ids serve as identifiers)
ids=1:3979 ##identifiers
originalECKs=read.csv(file="Data/allData.csv",colClasses=c("character", rep("NULL", 324))) %>% t %>% as.vector ##colClasses is used to skip all the other columns
id_ECKs<-cbind(ids,originalECKs) #Bind 2 columns




##Get the indices from ECK_name_column and bind -- this sapply() take about 10 sec
indices<-sapply(id_ECKs[,2],FUN=function(x){
  
  ###Escape all metas in regex
  x<-gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", x)
  
  length<-length(grep(x,ECK_name_columns[,"sorted_ECK_missing_gene_names_added"]))
  if(length<2&identical(grep(x,ECK_name_columns[,"sorted_ECK_missing_gene_names_added"]), integer(0))==FALSE){  
    ###if there are >=2 (duplicated ECK names), I want to deal with it later
    ###integer(0) is what gets returned when grep() finds nothing
    return(grep(x,ECK_name_columns[,"sorted_ECK_missing_gene_names_added"]))   
  }else(return(paste("Have to fix->",x))) 
}
)
indices<-unlist(indices) #simplify the result so all the integers are all in 1 list instead of many
indices<-unname(indices) #Get rid of names 

###manaully clean the duplicated strains (Total 12. But ECK2019 - HISA was found twice because there are: ECK2019 - HISA, ECK2019-HISA - SPA)
###
indices[indices=="Have to fix-> ECK1323-YMJC'"]<-c(3739,3740) ###The first index is from -1, the second from -2. 
indices[indices=="Have to fix-> ECK1824-MGRB"]<-c(3789,3790)
indices[indices=="Have to fix-> ECK3357-YHFL"]<-c(3901,3902)
indices[indices=="Have to fix-> ECK3531-DPPA"]<-c(3916,3917)
indices[indices=="Have to fix-> ECK2613-SMPA"]<-c(3855,3856)
indices[indices=="Have to fix-> ECK2019-HISA"]<-1708
indices[indices=="Have to fix-> ECK0295-YKGO"]<-c(3660,3661)
indices[indices=="Have to fix-> ECK1544-GNSB"]<-c(3763,3764)
indices[indices=="Have to fix-> ECK4410-YDGU"]<-c(3968,3969)
indices[indices=="Have to fix-> ECK4415-YPFM"]<-c(3970,3971)
indices[indices=="Have to fix-> ECK1556-HOKD"]<-c(3766,3767)
indices[indices=="Have to fix-> ECK4416-RYFB"]<-c(3972,3973)
indices[indices=="Have to fix-> ECK2593-A-YFIO\\* - Truncation"]<-c(132,133)

###convert indices to numberic
indices<-as.numeric(indices)


###order + cbind() + remove unnecessary columns from ECK_name_columns + correct the gene names (All capital -> the first 3 small)
ordered_ECK_name_columns<-ECK_name_columns[indices,]
id_ECKs_ordered_ECK_name_columns<-cbind(id_ECKs, ordered_ECK_name_columns)
id_ECKs_CorrectedECKs_AssociatedGeneNames<-id_ECKs_ordered_ECK_name_columns[,c("ids","originalECKs","sorted_ECK_missing_gene_names_added","associated_gene_names")]


id_ECKs_CorrectedECKs_AssociatedGeneNames[,"associated_gene_names"]<-sapply(id_ECKs_CorrectedECKs_AssociatedGeneNames[,"associated_gene_names"],function(x) ###(All capital -> the first 3 small)
{
  sub("[A-Za-z]{3}",str_extract(x,"[A-Za-z]{3}") %>% tolower,x)
}
)

id_ECKs_CorrectedECKs_AssociatedGeneNames[,"associated_gene_names"][id_ECKs_CorrectedECKs_AssociatedGeneNames[,"associated_gene_names"]=="TP2"]="tp2"  ###clean the exception
id_ECKs_CorrectedECKs_AssociatedGeneNames=as.data.frame(id_ECKs_CorrectedECKs_AssociatedGeneNames)



###(Updated 9/3/2019) manually correct the fake ECKs made by Nichols (search for gene names in EcoCyc and put back the right ECK)

###====================================================================================================================================================
###The list of 22 fake ECKs (I pulled them out by manually looking at the sorted ECK of the original data file that Nichols' provides) -> correct ECK
###ECK4466-MOKC ECK0018
###ECK4472-YOAI ECK1786
###ECK5000-SROH ECK4505
###ECK5001-SGRT ECK4477
###ECK5002-ISTR-1 G0-10202 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5003-RYBD ECK4621
###ECK5004-RYEF ECK4574
###ECK5005-TP2 G0-8894 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5006-TPKE70 G0-8906 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5007-YKGR ECK4486
###ECK5008-YMIB ECK4487
###ECK5009-YMJD ECK4488
###ECK5010-YNBG ECK4489
###ECK5011-YOAJ ECK4490
###ECK5012-YOAK ECK4491
###ECK5013-YOBI ECK4492
###ECK5014-YOEI ECK4493
###ECK5015-YOHP ECK4494
###ECK5016-YPDK ECK4495
###ECK5017-YQCG ECK4497
###ECK5018-YQEL ECK4498
###ECK5019-YQFG ECK4499


###correct the rows of ECK_1st_table for the above strains

##correct the 2nd, 3rd columns of "id_ECKs_CorrectedECKs_AssociatedGeneNames" generated at the end of clean_names.R (associated_gene_names shouldn't have to change)
fake_ECK_genes=c("ECK4466-MOKC","ECK4472-YOAI","ECK5000-SROH","ECK5001-SGRT","ECK5002-ISTR-1","ECK5003-RYBD","ECK5004-RYEF","ECK5005-TP2","ECK5006-TPKE70",
                 "ECK5007-YKGR","ECK5008-YMIB","ECK5009-YMJD","ECK5010-YNBG","ECK5011-YOAJ","ECK5012-YOAK","ECK5013-YOBI","ECK5014-YOEI","ECK5015-YOHP",
                 "ECK5016-YPDK","ECK5017-YQCG","ECK5018-YQEL","ECK5019-YQFG")

corrected_ECK=c("ECK0018","ECK1786","ECK4505","ECK4477","","ECK4621","ECK4574","","","ECK4486","ECK4487","ECK4488","ECK4489","ECK4490","ECK4491","ECK4492","ECK4493","ECK4494","ECK4495",
                "ECK4497","ECK4498","ECK4499")

temp=id_ECKs_CorrectedECKs_AssociatedGeneNames
for(i in 1:length(fake_ECK_genes)){
  index_=grep(fake_ECK_genes[[i]],id_ECKs_CorrectedECKs_AssociatedGeneNames$originalECKs)
  replacement=str_replace(fake_ECK_genes[[i]],"^ECK[0-9]{4}",corrected_ECK[[i]])
  if(corrected_ECK[[i]]=="") replacement=NA #if no ECK,give NA
  temp[index_,2]=replacement
  temp[index_,3]=temp[index_,2]

}



id_ECKs_CorrectedECKs_AssociatedGeneNames=temp
###====================================================================================================================================================


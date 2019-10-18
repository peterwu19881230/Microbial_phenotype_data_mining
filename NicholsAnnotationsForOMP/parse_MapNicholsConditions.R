#Create the table ready to make Nichols' annotations
#1. table for the virtual parent annotations
#2. table for the child strain annotations

#The google sheet: https://docs.google.com/spreadsheets/d/16ZYwtAqf0Yl8oN4erbLnehFpXGOPbmGmNFG8vhI7F2k/edit?pli=1#gid=404903003

currentDir=dirname(rstudioapi::getActiveDocumentContext()$path)

#I made the following file (Map Nichols conditions to OMP_for_parent_strain_Peter_9.19.2019.xlsx) was created by manually cleaning the above google doc that Debby made
###table_=read.xlsx(paste(currentDir,"Nichols_parent_term_annot.xlsx",sep="/"),sheetIndex = 1,startRow=1,encoding="UTF-8",quote="",header = TRUE)
table_=read_xlsx(paste(currentDir,"Nichols_parent_term_annot.xlsx",sep="/"),sheet=1)
##(on windows I have to specify the encoding. Otherwise I get weird characters: https://dss.iq.harvard.edu/blog/escaping-character-encoding-hell-r-windows)


table_[is.na(table_)]=""  #replace NA with "" 

class_=table_$childTerm

child_terms=unique(class_)
child_terms=child_terms[!is.na(child_terms)]

parent_terms=c("OMP:0007876 ! resistant to an antimicrobial agent",
               "OMP:0007247 ! resistant to a chemical",
               "OMP:0007373 ! resistant to alkylating agent",
               "OMP:0007909 ! resistant to a chelator",
               "OMP:0006039 ! resistant to macrolide",
               "OMP:0006065 ! resistant to tetracycline",
               "OMP:0006060 ! resistant to aminoglycoside",
               "OMP:0006047 ! resistant to beta-lactam",
               "OMP:0006043 ! resistant to quinolone",
               "OMP:0007171 ! resistant to antimicrobial peptide", #10
               "OMP:0007829 ! resistant to sulfonamide",
               "OMP:0007884 ! resistant to a metal",
               "OMP:0007892 ! resistant to a surfactant",
               "OMP:0007899 ! resistant to bile acids", 
               "OMP:0007517 ! presence of population growth in high osmolarity",
               "OMP:0007018 ! does not require growth factor(s)",
               "OMP:0007884 ! resistant to a metal",
               "OMP:0007962 ! presence of cell population growth during nutrient limitation",
               "OMP:0007306 ! presence of inorganic nitrogen source utilization",
               "OMP:0006024 ! presence of organic carbon source utilization", #20
               "OMP:0005190 ! presence of anaerobic growth", 
               "OMP:0007282 ! presence of UV radiation resistance",
               "OMP:0007956 ! presence of population growth at low temperature", 
               "OMP:0007950 ! presence of population growth at high temperature", 
               "OMP:0007945 ! presence of population growth at acidic pH",
               "OMP:0007771 ! presence of population growth at alkaline pH"
               )





#map back the parent terms to the corresponding row
row_index_for_parent_term=match(class_,child_terms)
parent_terms_mapped=parent_terms[row_index_for_parent_term]






#create virtual parent annotation table for .php script 
##page name to annotate to: OMP ST:74667 ! Virtual reference strain for Nichols RJ (2011) Cell 144, 143-156.
page_title=rep("OMP_ST:74667_!_Virtual_reference_strain_for_Nichols_RJ_(2011)_Cell_144,_143-156.",length(parent_terms_mapped))

annotation_table_for_parent=data.frame(page_title=page_title)
library(stringr)
annotation_table_for_parent$OMP_ID=str_extract(parent_terms_mapped,"OMP:[0-9]{7}")
annotation_table_for_parent$OMP_term=str_replace(parent_terms_mapped,"OMP:[0-9]{7} ! ","")
annotation_table_for_parent$relative_phenotype=""

chemical=table_$condition
annotation_table_for_parent$experimental_condition=paste0("\\n*temperature:37C\\n*medium:",table_$medium,"\\n*medium:",chemical) #temperature + medium + chemical name
##Note 1: I have to have \\n at the very beginning to get the proper bullet point for "temperature:37C"
##Note 2: I have to do \\n to escape otherwise the printed .txt will have the wrong format


annotation_table_for_parent$ECO_ID="ECO:0007005"
annotation_table_for_parent$ECO_term="high throughput direct assay evidence used in manual assertion"
annotation_table_for_parent$ref_PMID="PMID:21185072"

CHEBI_ID=table_$chemicalID
annotation_table_for_parent$Annot_Extension=paste0("towards some ","(",CHEBI_ID,")"," ",chemical) #(1st line): towards(CHEBI) (2nd line): SDS
annotation_table_for_parent$Annot_Extension[grep("()",annotation_table_for_parent$Annot_Extension,fixed=T)]="" #clean conditions that don't have CHEBI ID
annotation_table_for_parent$Notes=table_$concentration #concentrations or other details
annotation_table_for_parent$Notes=sapply(annotation_table_for_parent$Notes,FUN=function(note_){
  note_=stringr::str_replace_all(note_,"µ","u") #change "µ" to "u" ("µ" doesn't display on wiki)
  note_=stringr::str_replace_all(note_,",",", ") #add extra space after ","
  note_=stringr::str_replace_all(note_,"µ","") #remove "µ" before C
})
annotation_table_for_parent$Notes[grep("percent",annotation_table_for_parent$Notes)]=annotation_table_for_parent$Notes[grep("percent",annotation_table_for_parent$Notes)] %>% str_replace("percent","%")



write.table(annotation_table_for_parent,file=paste(currentDir,"annotation_table_for_parent.txt",sep="/"),sep="\t",quote=F,fileEncoding="UTF-8",row.names = F,col.names = F)
#write.csv(annotation_table_for_parent,file=paste(currentDir,"annotation_table_for_parent.csv",sep="/"),fileEncoding="UTF-8",row.names = F) #even if I change the file encoding to UTF-8 it still outputs weird characters




#This part creates the table for child strain annotations


##ref: https://stackoverflow.com/questions/15347282/split-delimited-strings-in-a-column-and-insert-as-new-rows
temp=table_
temp[,"childTerm"]=parent_terms_mapped; names(temp)[5]="parentTerm"  #replace childTerm column with parent term column, correct the column name


temp[temp$name_in_the_R_object_I_defined__uniqueChemIndex=="Anaerobic","concentration"]=" " #give " " for anaerobic condition because "" doesn't split when using strsplit() 
row_id=seq(dim(temp)[1])
s=strsplit(temp$concentration,split=",")

##if there are >= 2 concentrations, add the unit to all the concentrations (except for the last one I don't have to)
clean_s=lapply(s,FUN=function(concentrations){
  if(length(concentrations)>=2 & 
     !identical(concentrations,c("0.5%/0.1mM","0.5%/0.5mM","1%/0.5mM")) & 
     !identical(concentrations,c("pH 4","4.5","5","6")) &
     !identical(concentrations,c("pH 8","9","9.5","10"))
     ){ #I excluded 3 exceptions in this "if" because the following regex only works on the rest of the cases
    unit_=tail(concentrations,1) %>% str_replace(pattern="[0-9.]{1,10}","") %>% trimws 
    return(c(paste(head(concentrations,-1),unit_),
             tail(concentrations,1))
           )
  }else{
    if(identical(concentrations,c("0.5%/0.1mM","0.5%/0.5mM","1%/0.5mM"))) return(concentrations)
    if(identical(concentrations,c("pH 4","4.5","5","6"))) return(c("pH 4","pH 4.5","pH 5","pH 6"))
    if(identical(concentrations,c("pH 8","9","9.5","10"))) return(c("pH 8","pH 9","pH 9.5","pH 10"))
    return(concentrations)
  }
})
clean_s=unlist(clean_s)
clean_s=clean_s %>% str_replace("µ","u") %>% str_replace("°","") %>% str_replace("percent","%")


temp=base::subset(temp,select=-c(concentration))
temp=temp[rep(row_id,sapply(s,length)),]
temp$splitted_concentration=clean_s

parentTerm_childTermIncreased_childTermDecreased=temp[,c("parentTerm","childTerm_increased","childTerm_decreased")]

#no. of concentration=334 because of 10 extra annotations made to minimal media


##For decreased resistance pick the lowest concentration, for increased resistance pick the highest

new_uniqueChemIndex=uniqueChemIndex
names(new_uniqueChemIndex)[grep("Calcofluor",names(new_uniqueChemIndex))]="Calcofluor (F3543 Fluorescent Brightener 28)" #remove a EOL character that is in uniqueChemIndex just for this script
new_uniqueChemIndex=new_uniqueChemIndex[match(table_$name_in_the_R_object_I_defined__uniqueChemIndex,names(new_uniqueChemIndex))] %>% unlist



temp=cbind(chem_index=new_uniqueChemIndex,temp) 
 
  
  

##now the above "temp" table contains these columns:
###[1] "chem_index"                                      "name_in_the_R_object_I_defined__uniqueChemIndex" "condition"                                      
###[4] "chemicalID"                                      "childTerm"                                       "childTerm_increased"                            
###[7] "childTerm_decreased"                             "medium"                                          "splitted_concentration"        



#Check:
#concentration column => checked
#chem index column => checked


#create this table that only contains significant phenotypes: strain ID - condition ID  => left_join "temp" table by chem_index but truncate "increased", "decreased" term columns into 1

annotation_table_for_children=data.frame()
for(strain_id in 1:3979){
  for(condition_id in 1:324){
    sig_fitness=Ternary_Data_324cutff_NAremoved[strain_id,condition_id]
    fitness=All_Data_NAimputed[strain_id,condition_id]
    
    if(sig_fitness==1){
      annot="childTerm_increased"
      
    }else if(sig_fitness==-1){
      annot="childTerm_decreased"
    
      }else{
      next   
      
        }
    
    rows=cbind(strain_id,condition_id,temp[temp$chem_index==condition_id,c("name_in_the_R_object_I_defined__uniqueChemIndex","condition","chemicalID",annot,"medium","splitted_concentration")])
    rows$fitness=round(fitness,2)
    colnames(rows)=c("strain_id","condition_id","name_in_the_R_object_I_defined__uniqueChemIndex","condition","chemicalID","child_term","medium","concentration","fitness") 
    
    annotation_table_for_children=rbind(annotation_table_for_children,
                                        rows) 
  }
  
}


#Filter out annotations based on this: For "increased" phenotypes only annotate the highest concentration, for "decreased" only annotate the lowest concentration
#=> For the same strain + same condition (but different concentration) + same OMP term, only annotate once


start.time=Sys.time()

temp=data.frame()
for(strain_id in unique(annotation_table_for_children$strain_id)){
  
  subset1=annotation_table_for_children[annotation_table_for_children$strain_id==strain_id,]
  
  for(condition in unique(subset1$condition)){
    
    subset2=subset1[subset1$condition==condition,]
    
    for(child_term in unique(subset2$child_term)){
      
      subset3=subset2[subset2$child_term==child_term,]
      
      
      if( sum(grep("! increased",subset3$child_term,fixed=T))==length(subset3$child_term) ){ #this logical setup can ignore strains that have both "increased" and "decreased" conflict conditions
        subset3=subset3[subset3$condition_id==max(subset3$condition_id),] #this is based on assuming the latter ids contain higher concentrations
      }else if(sum(grep("! decreased",subset3$child_term,fixed=T))==length(subset3$child_term)){
        subset3=subset3[subset3$condition_id==min(subset3$condition_id),]
      }else{
        next
      }
      
      temp=rbind(temp,subset3)
      
      
    }
  }
  
}


annotation_table_for_children=temp

end.time=Sys.time()
end.time-start.time #Time difference of 16.0336 secs








##map in the strain page names on OMP

##from  2018-07-14_ReadyToMakeStrainPages.tsv, id_allAttributes: map strain page names to annotation_table_for_children using ECK + manual corrections


filePath=paste(currentDir,"2018-07-14_ReadyToMakeStrainPages.tsv",sep="/")
fromDebby=read.table(filePath,sep="\t",fill=T,quote="",header=T)


id_ECK=id_allAttributes[,c("ids","ECK")] %>% unique
names(id_ECK)[2]="corrected.ECK.ID"

temp=left_join(id_ECK,fromDebby,by="corrected.ECK.ID") 
dim(temp)



## see what's going on with the ids that map to more than 1 ECK (or no match so the result is NA)
index_for_problematic_ids=which(duplicated(temp$ids,fromLast=F) | duplicated(temp$ids,fromLast=T))
duplicated_rows=temp[index_for_problematic_ids,] #some duplications resulted from: Same ECK, same id, but different alleles (e.g. look at ECK0055)

## manually clean this table should be enough: duplicated_rows
id_=duplicated_rows$ids %>% unique #all unique ids that need cleaning


correct_allele=c("lptD::DAS-kan","lptD(imp4213)-kan","lptD+-kan","lptD::DAS+4-kan", 
                 "murE+-kan","murE::SPA-kan","murE-kan",
                 "ftsA::SPA","ftsA+-kan","ftsA(R286W)-kan",
                 "lpxC+-kan","lpxC(G210S)-kan","lpxC::SPA",
                 "bamA66-kan","bamA+-kan","bamA6-kan",
                 "fabZ+-kan","fabZ(F101Y)-kan",
                 "lolA::SPA","lolA::DAS-kan","lolA::DAS+4-kan",
                 "msbA+-kan","msbA(P18S)-kan",
                 "hisA::SPA-kan","hisA722(del)::FRT-kan-FRT",
                 "yhjD+-kan","yhjD777(del)::FRT-kan-FRT",
                 "kdtA::SPA","kdtA+-kan","kdtA(E234V)-kan"
                 )


corrected_sub_table=data.frame(ids=id_,Allele.at.CGSC=correct_allele)
corrected_sub_table=left_join(corrected_sub_table,temp,by=c("ids","Allele.at.CGSC"))

rows_that_has_no_problem=temp[-index_for_problematic_ids,]

id_strainPageTable=rbind(rows_that_has_no_problem,corrected_sub_table)


#get a table that has 2 columns: 1. Nichols' ids 2. strain page names (not including the OMP identifiers because they are difficult to get) <= this will be used to search unique strain pages in OMP



#2 kinds of strain pages:

#if there is "allele in column D" -> Escherichia coli K-12 + (search and replace "allele in column D" by "Allele.at.CGSC")
#if not just use  Escherichia coli K-12 + "Name.for.Strain.Page"


strain_page_name=c()
for(i in 1:dim(id_strainPageTable)[1]){
  
  if(grepl("allele in column D",id_strainPageTable$Name.for.Strain.Page[i])){
    strain_page_name[i]=paste("Escherichia coli K-12",str_replace(id_strainPageTable$Name.for.Strain.Page[i],"allele in column D",id_strainPageTable$Allele.at.CGSC[i]))
  }else if(is.na(id_strainPageTable$Allele.at.CGSC[i])){
    strain_page_name[i]=NA #assign NA and remove later
  }else{
    strain_page_name[i]=paste("Escherichia coli K-12",id_strainPageTable$Allele.at.CGSC[i])
  }
  
  strain_page_name[i]=str_replace_all(strain_page_name[i]," ","_") #have to replace " " to "_" so SQL search would work
  
}


id_nameForStrainPage=data.frame(id=id_strainPageTable$ids,strain_page_name=strain_page_name )
id_nameForStrainPage=id_nameForStrainPage[!is.na(id_nameForStrainPage$strain_page_name),]
dim(id_nameForStrainPage)



#why 3955 strain pages propagating to -> 3966: First reason: 12 Duplicated strains in Nichols'
temp2=id_nameForStrainPage[duplicated(id_nameForStrainPage$strain_page_name,fromLast=F) | duplicated(id_nameForStrainPage$strain_page_name,fromLast=T),]
propogated_ids=temp2$id

##check the 12 duplicated strains in Nichols'
id_originalName=id_allAttributes[,c("ids","Original Name")] %>% unique
dim(id_originalName)
duplicated_ids=id_originalName$ids[duplicated(id_originalName$`Original Name`,fromLast=T) | duplicated(id_originalName$`Original Name`,fromLast=F)]


### check what are ids that weren't propagated because of those 12 duplicated strains
temp2[!(propogated_ids %in% duplicated_ids),] 

id_strainPageTable[id_strainPageTable$ids %in% as.character(c(3514,367,2923,3862,3198,3315,3333)),c("ids","corrected.ECK.ID","Allele.at.CGSC")]

#(id=367, ECK2652), (id=2923, ECK2651) not found on EcoCyc

#(id=3862, original name: ECK4142-4143-ECNAB) should be annotated differently (Now it's been linked to this allele: ecnA773(del)::FRT-kan-FRT)
#(id=3315, original name: ECK4418/2590-RYFD(WITHCLPB)) should be annotated differently (Now it's been linked to this allele: ryfD(del)(clpB intact)::FRT-kan-FRT)


#temporarily ignore id=367, 3862, 3315
id_nameForStrainPage=id_nameForStrainPage[!(id_nameForStrainPage$id %in% c("367","3862","3315")),]


#create the final table
names(id_nameForStrainPage)[1]="strain_id"
id_nameForStrainPage$strain_id=as.integer(id_nameForStrainPage$strain_id)
temp3=inner_join(id_nameForStrainPage,annotation_table_for_children)
dim(temp3)
dim(annotation_table_for_children)
annotation_table_for_children=temp3


##Get this mapping: parent_ANNO_id - parent OMP ID - child OMP ID

#I ran this on BCBP to get the parent annotation table:
##Run get_parent_annotations_PW.php by: php php get_parent_annotations_PW.php -w /var/www/html/omp/wiki -t Virtual_reference_strain_for_Nichols_RJ -f vitural_parent.csv
parent_annot_from_OMP=read.csv(paste(currentDir,"vitural_parent.csv",sep="/"),header = F)
annotID_parentOMPID=parent_annot_from_OMP[,c(1,3)]

annotID_parentOMPID_childOMPIDandTerm=cbind(annotID_parentOMPID,table_[,c("childTerm_increased","childTerm_decreased")])



##Final form for child strain annotations:  
#-- page title
#-- OMP ID 
#-- OMP term name (can leave blank)
#-- relative phenotype info (parent term annotation ID)
#-- experimental condition (medium, temperature...etc)
#-- ECO ID
#-- ECO term name (can leave blank)
#-- PMID
#-- annotation extension (leave blank)
#-- notes (fitness score: (number))

#annotation_table_for_children[,c("strain_page_name","child_term","","")]...


#start from here and finish...
annotID_parentOMPID_childOMPIDandTerm
parentTerm_childTermIncreased_childTermDecreased
annotation_table_for_parent
annotation_table_for_children



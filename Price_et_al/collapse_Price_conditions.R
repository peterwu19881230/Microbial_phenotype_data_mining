
#Place to download the data: http://genomics.lbl.gov/supplemental/bigfit/
#An earlier version (with only 92 successful conditions) from the previous paper is here: http://genomics.lbl.gov/supplemental/rbarseq/html/Keio/

dat=read.table("Data/fit_logratios_good.tab.txt",quote="",fill=T,sep="\t",header = T,check.names=F) 
dat_scoresOnly=dat[,-(1:4)]

conditions=names(dat_scoresOnly)

#Supplementary Table 18: Supplementary Table 18 – Media formulations used in this study
##download all sub tables in one excel file: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0124-0/MediaObjects/41586_2018_124_MOESM3_ESM.xlsx
conditions_tag_removed=str_replace(conditions,"^[a-zA-Z0-9]{9} ","")


###(!)Not sure why they have multiple D-Glucose, D-Fructose...etc without labeling concentrations 
##(are they repeats? Here I temporarily treat them as different concentrations)

##for the identical names I will use a loop to map them
collapsed_condition_index=list("D-Glucose (C)"=1:2,"D-Fructose (C)"=3:4,"Sucrose (C)"=5:6,"D-Maltose (C)"=7:8,
                               "D-Xylose (C)"=9:10,"D-Galactose (C)"=11:12,"D-Ribose (C)"=13:14,"L-Fucose (C)"=15:16,
                               "a-Ketoglutaric (C)"=17:18,"D-Glucuronic Acid (C)"=19:20,"D-Gluconic Acid (C)"=21:22,
                               "D-Glucose-6-Phosphate (C)"=23:24,"acetate (C)"=25:26,"D-Lactate (C)"=27:28,"D,L-Lactate (C)"=29:30,
                               "pyruvate (C)"=31:32,"succinate (C)"=33:34,"Glycolic Acid (C)"=35:36,"D-Mannitol (C)"=37:38,
                               "D-Galacturonic Acid (C)"=39:40,"Glycerol (C)"=41:42,"D-Sorbitol (C)" =43:44,"L-Malic (C)"=45:46,
                               "D-Trehalose (C)"=47:48,"D-Serine (C)"=49:50,"casaminos (C)"=51:52,"NAG (C)"=53:54,
                               "D-Glucosamine Hydrochloride (C)"=55:56,"D-Mannose (C)"=57:58,"D-Galacturonic Acid (C) MOPS"=59:60,
                               "D-Glucose (C) MOPS"=61:62,"L-Arginine (N)"=63:64,"L-Aspartic Acid (N)"=65:66,"L-Serine (N)"=67:68,
                               "L-Asparagine (N)"=69:70,"L-Glutamine (N)"=71:72,"Glycine (N)"=73:74,"L-Alanine (N)"=75:76,
                               "D-Serine (N)"=77:78,"D-Alanine (N)"=79:80,"Gly-DL-Asp (N)"=81:82,"Gly-Glu (N)"=83:84,
                               "casaminos (N)"=85:86,"Putrescine (N)"=87:88,"Adenosine (N)"=89:90,"Cytidine (N)"=91:92,
                               "Ammonium chloride (N)"=93:94,"nitrite"=95:96,"Chlorite"=97:98,
                               "Chloride"=99:100,
                               
                               "A22"=101:102,"Doxycycline hyclate"=103,"Nickel (II) chloride"=104,"Cobalt chloride"=105,
                               "copper (II) chloride"=106,"sodium fluoride"=107,"Thallium(I) acetate"=108,"Cisplatin"=109:111,
                               "LB"=c(112,131),"Benzalkonium Chloride"=113,"Aluminum chloride"=114,"benzoic"=115,"Bacitracin"=116:117,
                               "Fusidic"=118:119,"Dimethyl Sulfoxide"=120,"Dimethyl Sulfoxide"=121,
                               "Chloramphenicol"=122:124,"Tetracycline"=125:126,"Spectinomycin"=127:128,"Carbenicillin"=129:130,
                               "methylglyoxal"=132,"2-Furfuraldehyde"=133:134,"5-Hydroxymethylfurfural"=135,"Vanillin"=136,
                               "syringaldehyde"=137:138,"Dimethyl Sulfoxide"=139,"Cholin acetate"=140:141,"acetate"=142:143,
                               "1-ethyl-3-methylimidazolium acetate"=144:145,"1-ethyl-3-methylimidazolium chloride"=146:147,
                               "L-Lysine"=148,"D-Galacturonic Acid (C) MOPS"=149,"D-Glucose (C) MOPS"=150,
                               "D-Galacturonic Acid (C) M9"=151,"D-Glucose (C) M9"=152,"Nalidixic"=153,"Phosphomycin"=154,
                               "D-Cycloserine"=155,"outer cut, LB soft agar motility assay" =156:158,
                               "inner cut, LB soft agar motility assay"=159:162
)


#double check that I have included every index
#=============================================
all_index_in_list=c()
for(i in 1:length(collapsed_condition_index)){
  all_index_in_list=c(all_index_in_list,collapsed_condition_index[[i]])
}

sum(base::sort(all_index_in_list)==1:162)==162 #TRUE
#=============================================


#To collapse conditions: For nutrients: average the scores. For stresses: pick the severest phenotype
##double check the pcc of the duplicated nutrients (up to index=94)
i=1
pcc=c()
while(names(collapsed_condition_index)[i]!="nitrite"){
  
  cond1=dat_scoresOnly[collapsed_condition_index[[i]]][1]
  cond2=dat_scoresOnly[collapsed_condition_index[[i]]][2]
  
  pcc=c(pcc,cor(cond1,cond2))
  i=i+1
}

summary(pcc) 
hist(pcc)
##=>all of them have high pcc (then why did they keep the replicates instead of averaging them?)


##create a quantitative data with averaged fitness for nutrients
nutrientOnly_odd=dat_scoresOnly[,1:94][,as.logical(1:94%%2==1)]
nutrientOnly_even=dat_scoresOnly[,1:94][,as.logical(1:94%%2==0)]
nutrientOnly_averaged=(nutrientOnly_odd+nutrientOnly_even)/2  
names(nutrientOnly_averaged)=paste(names(collapsed_condition_index)[1:(94/2)],"-averaged",sep="")


price_nutrientsAveraged=cbind(nutrientOnly_averaged,dat_scoresOnly[,95:dim(dat_scoresOnly)[2]])

#save(price_nutrientsAveraged,file="Data/sourced/price_nutrientsAveraged.RData")


#Next: map conditions to Nichols (find overlapped conditions)
#==The mapping is in: Nichols_Price_overlappedConds.docx ==

##test whether the fitness under the same condition for the same strain are similar or not



##strain mapping (modified from correlationVSannotation_test.R)
#====================================================================================
dat$sysName=as.character(dat$sysName)
names(dat)
names(dat)[2]="bNumber"
duplicated(dat$bNumber) %>% sum #No duplicates in Deutschbauer's data (Awesome!)

bNumber_inDeutschbauer=cbind(dat$bNumber,T) %>% as.data.frame(stringsAsFactors=F)
names(bNumber_inDeutschbauer)=c("bNumber","in_Deutschbauer_s")

temp=left_join(id_allAttributes,bNumber_inDeutschbauer,by = "bNumber")
temp$in_Deutschbauer_s[is.na(temp$in_Deutschbauer_s)]=F

temp2=unique(temp[,c("ids","in_Deutschbauer_s")])

sum(temp2$in_Deutschbauer_s %>% as.logical) #3525 strains mapped to Nichols



##The id for Nichols' / bNumber for Deutschbauer's
strain_to_merge=temp[as.logical(temp$in_Deutschbauer_s)==T,c("ids","bNumber")] %>% unique
##Nichols' has duplicates => some of the Deutschbauer's will be duplicated

##get the indices for subtracting Deutschbauer's data
index=c()
count=0
for(bNumber in strain_to_merge$bNumber){
  count=count+1
  index[count]=which(dat$bNumber==bNumber)
}

Nich_Deut=cbind(All_Data_NAimputed[strain_to_merge$ids %>% as.numeric,],
                dat_scoresOnly[index,])

Deut=dat_scoresOnly[index,] #This only contains strains that are also in Nichols'
rownames(Deut)=rownames(Nich_Deut)
str(Deut)


#=====Data that have the same strains mached on rows=======
Nich=All_Data_NAimputed[strain_to_merge$ids %>% as.numeric,]
Price=dat_scoresOnly[index,]; rownames(Price)=rownames(Nich)  
#==========================================================


#====================================================================================





##condition mapping

###collapse the quantitative data using the most significant fitness scores -> compare

collapse_condition=function(cols_){ #function to select the most severe fitness across the same condition with different concentrations
  
  if(!is.null(dim(cols_))){
    most_severe=apply(cols_,1,FUN=function(x){
      
      result=x[which(abs(x)==max(abs(x)))]
      
      if(length(result)>1){ #if there are 2 same fitness scores (eg. when there are same inputed fitness score), only pick 1
        result=result[1]
      }
      
      return(result)
      
      
    } 
    )
  }else{
    most_severe=cols_
  }
  
  return(most_severe)
}


#acetate (carbon source)
cor(collapse_condition(Nich[,19]),collapse_condition(Price[,25:26])) #0.58




#function to get tables that show most different fitness between the 2 papers
#============================================================================
id_ECKandGene_all=unique(id_allAttributes[,c("ids","sorted_ECK_missing_gene_names_added")])
#dim(id_ECKandGene_all) #3979    2
compute_difference_table=function(collapsedCond1,collapsedCond2){
  diff=collapsedCond1-collapsedCond2
  rank_=base::rank( -1*abs(diff) ) #-1 is to let rank work in reverse
  ids=rownames(Nich)
  ECKandGene=id_ECKandGene_all$sorted_ECK_missing_gene_names_added[match(ids,id_ECKandGene_all$ids)]
  diff_sorted=cbind(ECKandGene,collapsedCond1,collapsedCond2,diff)[order(rank_),]
  
  return(diff_sorted)
}
#============================================================================

#function to compute pcc for all combinations of conditions
##============================================================================
pairwise_pcc_of_conds=function(conds1,conds2){
  cors=c()
  
  conds1=as.matrix(conds1)
  conds2=as.matrix(conds2)
  
  for(i in 1:dim(conds1)[2]){
    for(j in 1:dim(conds2)[2]){
      cors=c(cors,cor(conds1[,i],conds2[,j]))
    }
  }
  
  return(cors)
}
#============================================================================


#functions to calculate similarity based on direction (positive or negative) of the phenotypes
##============================================================================

###similarity is defined as: no. of same element/total no. of elements

to_ternary=function(conds){
  conds[conds>0]=1
  conds[conds<0]=-1
  
  return(conds)
}


cal_sim_ternary=function(conds1,conds2){
  sims=c()
  
  conds1_ternary=as.matrix(conds1) %>% to_ternary
  conds2_ternary=as.matrix(conds2) %>% to_ternary
  
  for(i in 1:dim(conds1_ternary)[2]){
    for(j in 1:dim(conds2_ternary)[2]){
      sims=c(sims,sum(conds1_ternary[,i]==conds2_ternary[,j])/length(conds1_ternary))
    }
  }
  
  return(sims)
}
#============================================================================


#benzalkonium Chloride
cor(collapse_condition(Nich[,40:42]),collapse_condition(Price[,113])) #0.1967572
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,40:42]),collapse_condition(Price[,113]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,40:42],Price[,113]) # -0.02778135  0.16882590  0.24250316
cal_sim_ternary(Nich[,40:42],Price[,113]) #0.1639716 0.1694563 0.1739007



#cisplatin
cor(collapse_condition(Nich[,84:86]),collapse_condition(Price[,109:111])) #0.004155488
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,84:86]),collapse_condition(Price[,109:111]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,84:86],Price[,109:111]) #-0.0120499079 -0.0089946784  0.0074761643 -0.0216977432 -0.0224078298 -0.0067128991 -0.0172739081 -0.0009027592  0.0029812689
cal_sim_ternary(Nich[,84:86],Price[,109:111]) #0.1653901 0.1667139 0.1636879 0.1667139 0.1655792 0.1658629 0.1641608 0.1669976 0.1648227


#cobalt chloride
cor(collapse_condition(Nich[,91:92]),collapse_condition(Price[,105])) #0.03665826
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,91:92]),collapse_condition(Price[,105]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,91:92],Price[,105])  #0.02568027 0.03118506
cal_sim_ternary(Nich[,91:92],Price[,105]) #0.2537589 0.2486525


#copper (II) chloride
cor(collapse_condition(Nich[,93:95]),collapse_condition(Price[,106])) #0.2324331
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,93:95]),collapse_condition(Price[,106]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,93:95],Price[,106]) #0.07325352 0.12302751 0.29762860 
cal_sim_ternary(Nich[,93:95],Price[,106]) #0.1685106 0.1689835 0.1753191


#glucose (carbon source)
cor(collapse_condition(Nich[,129]),collapse_condition(Price[,1:2])) #0.5915173

#glycerol (carbon source)
cor(collapse_condition(Nich[,130]),collapse_condition(Price[,41:42])) #0.6187663

#maltose (carbon source)
cor(collapse_condition(Nich[,140]),collapse_condition(Price[,7:8])) #0.6057956

#nickel (II) chloride
cor(collapse_condition(Nich[,152:153]),collapse_condition(Price[,104])) #0.1685807
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,152:153]),collapse_condition(Price[,104]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,152:153],Price[,104]) #-0.0153105  0.2233411
cal_sim_ternary(Nich[,152:153],Price[,104]) #0.2487943 0.2609929




#succinate (carbon source)
cor(collapse_condition(Nich[,204]),collapse_condition(Price[,33:34])) #0.5410673

#bacitracin
cor(unlist(collapse_condition(Nich[,234:236])),collapse_condition(Price[,116:117])) #0.1027618
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,234:236]),collapse_condition(Price[,116:117]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,234:236],Price[,116:117]) #0.05356377 0.01007168 0.10488942 0.05940903 0.20135804 0.11219441
cal_sim_ternary(Nich[,234:236],Price[,116:117]) #0.1689835 0.1683215 0.1713475 0.1674704 0.1738061 0.1673759


#carbenicillin
cor(collapse_condition(Nich[,237:239]),collapse_condition(Price[,129:130])) #0.1294743
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,237:239]),collapse_condition(Price[,129:130]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,237:239],Price[,129:130]) #0.09399372 0.07726313 0.19142419 0.14161970 0.15057777 0.11488754
cal_sim_ternary(Nich[,237:239],Price[,129:130]) #0.1683215 0.1683215 0.1734279 0.1707801 0.1734279 0.1717258



#chloramphenicol
cor(collapse_condition(Nich[,247:250]),collapse_condition(Price[,122:124])) #0.1958572
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,247:250]),collapse_condition(Price[,122:124]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,247:250],Price[,122:124]) #0.10440896 0.11186432 0.06726894 0.17137299 0.16902621 0.09584131 0.23285460 0.21940823 0.12522696 0.27594246 0.24809851 0.13136647
cal_sim_ternary(Nich[,247:250],Price[,122:124]) #0.1288652 0.1268085 0.1263830 0.1285106 0.1290071 0.1282270 0.1313475 0.1314184 0.1295035 0.1328369 0.1329078 0.1295745



#D-Cycloserine
cor(collapse_condition(Nich[,254]),collapse_condition(Price[,155])) #0.08853055
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,254]),collapse_condition(Price[,155]))[1:30,]
cal_sim_ternary(Nich[,254],Price[,155]) #0.5228369


#Doxycycline 
cor(collapse_condition(Nich[,255:258]),collapse_condition(Price[,103])) #0.2411233
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,255:258]),collapse_condition(Price[,103]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,255:258],Price[,103]) #0.1803125 0.2301841 0.2377179 0.2659150
cal_sim_ternary(Nich[,255:258],Price[,103]) #0.1272340 0.1310638 0.1305674 0.1290071



#Nalidixic acid
cor(collapse_condition(Nich[,274:277]),collapse_condition(Price[,153])) #0.2459178
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,274:277]),collapse_condition(Price[,153]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,274:277],Price[,153]) #-0.02306796  0.15038612  0.26897900  0.34682549
cal_sim_ternary(Nich[,274:277],Price[,153]) #0.1224823 0.1278014 0.1265248 0.1312057



#Spectinomycin
cor(collapse_condition(Nich[,304:305]),collapse_condition(Price[,127:128])) #0.2341706
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,304:305]),collapse_condition(Price[,127:128]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,304:305],Price[,127:128]) #0.2529728 0.2607911 0.2233260 0.2335108
cal_sim_ternary(Nich[,304:305],Price[,127:128]) #0.2605674 0.2544681 0.2643972 0.2587234


#Tetracycline
cor(collapse_condition(Nich[,312:315]),collapse_condition(Price[,125:126])) #0.1773482
##Genes that have the most different scores under the above overlapped conditinos
result=compute_difference_table(collapse_condition(Nich[,312:315]),collapse_condition(Price[,125:126]))[1:30,]
##pcc of individual condition
pairwise_pcc_of_conds(Nich[,312:315],Price[,125:126]) #0.06133567 0.09343241 0.15233124 0.19006481 0.21767751 0.23340651 0.18580539 0.20673399
cal_sim_ternary(Nich[,312:315],Price[,125:126]) #0.1243972 0.1257447 0.1256028 0.1275177 0.1295745 0.1293617 0.1293617 0.1293617


#summary google sheet:
##https://docs.google.com/spreadsheets/d/14HhhuU6UmCxsYdaeEnkPdABk0rAfJ5DqGEFJwc4VlCw/edit#gid=0









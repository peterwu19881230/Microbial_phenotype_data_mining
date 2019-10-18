pos=read.table("Data/genomwide_metabolomic/zscore_pos.tsv")
neg=read.table("Data/genomwide_metabolomic/zscore_neg.tsv")


all=cbind(t(pos),t(neg)) #id from top to bottom matches those in the excel?
##=>comfirmation using several examples:
##By Fig.4 in the paper (although they don't seem to have used pcc, the concept should be similar)
##The paper that describe the method used in Fig.4 (CLR algorithm <= it uses mutual information):
##https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050008
##The corresponding R packages: 
###https://rdrr.io/cran/parmigene/man/clr.html
###https://rdrr.io/bioc/minet/man/clr.html

##, the following genes have high pairwise pcc:
##1. fepB, fepD, fepG, fepC: 566,568,570,567
##2. sdhA, sdhB, sdhC, sdhD: 1761:1764
##3. cyoA, cyoB, cyoC, cyoD: 355:358
##These 2 pairs should have low pcc: (frdA, dppF): 650,457     (fepE, fepG): 569, 570     

#Calculate similarity based on PCC:

##if id (1~3807) is from top to bottom:
cor(t(all[c(566,568,570,567),]))
cor(t(all[1761:1764,]))
cor(t(all[355:358,]))
cor(t(all[c(650,457),]))
cor(t(all[c(569,570),]))
##Almost all pccs are not as expected!

##if id (1~3807) is from bottom to top:
cor(t(all[3807+1-c(566,568,570,567),]))
cor(t(all[3807+1-1761:1764,]))
cor(t(all[3807+1-355:358,]))
cor(t(all[3807+1-c(650,457),]))
cor(t(all[3807+1-c(569,570),]))
##Almost all pccs are not as expected!


#Calculate similarity based on MI:
if(!require(mpmi)){
  install.packages("mpmi")
  library(mpmi)
}

##if id (1~3807) is from top to bottom:
cminjk(t(all[c(566,568,570,567),]))
cminjk(t(all[1761:1764,]))
cminjk(t(all[355:358,]))
cminjk(t(all[c(650,457),]))
cminjk(t(all[c(569,570),]))
##Diagnals have MI ~ 2, other MIs are generally low (e.g. 0.03)

##if id (1~3807) is from bottom to top:
cminjk(t(all[3807+1-c(566,568,570,567),]))
cminjk(t(all[3807+1-1761:1764,]))
cminjk(t(all[3807+1-355:358,]))
cminjk(t(all[3807+1-c(650,457),]))
cminjk(t(all[3807+1-c(569,570),]))
##Diagnals have MI ~ 2, other MIs are generally low (e.g. 0.03)


#try minet package that calculates MI (https://rdrr.io/bioc/minet/man/build.mim.html)
if(!require(minet)){
  BiocManager::install("minet")
}


start.time = Sys.time()
mi_all=build.mim(t(all)) #the default estimator is spearman, but I am not sure it's correct (or acceptable)
end.time = Sys.time()
end.time - start.time #Time difference of 1.16055 mins


##Calculate similarity based on CLR:
net = clr(mi_all)

##if id (1~3807) is from top to bottom:
net[c(566,568,570,567),c(566,568,570,567)]
net[1761:1764,1761:1764]
net[355:358,355:358]
net[c(650,457),c(650,457)]
net[c(569,570),c(569,570)]


##if id (1~3807) is from bottom to top:
net[3807+1-c(566,568,570,567),3807+1-c(566,568,570,567)]
net[3807+1-1761:1764,3807+1-1761:1764]
net[3807+1-355:358,3807+1-355:358]
net[3807+1-c(650,457),3807+1-c(650,457)]
net[3807+1-c(569,570),3807+1-c(569,570)]

##(!)the diagnal elemetns are all 0 doesn't make sense to me for now






gene_names=read_xls("Data/genomwide_metabolomic/sample_id_zscore.xls",sheet=1,col_names=F)




##compute pairwise MI matrix from Fuhrers
if(!require(mpmi)){
  install.packages("mpmi")
  library(mpmi)
}

start.time = Sys.time()
mi_all=cminjk(t(all))
end.time = Sys.time()
end.time - start.time #running only the 1:10 features takes 7 sec, 1:20 takes 9 sec, 1:100 takes 2 min. Takes forever to run the fullset on my PC
fuhrer_mi=mi_all
save(fuhrer_mi,file="Data/fuhrer_mi.RData")


#try to recover some associtaion of genes in the Fuhrer paper to decide the mapping of gene names to the z-score data
load("Data/fuhrer_mi.RData")

if(!require(parmigene)){
  install.packages("parmigene")
  library(parmigene)
}

temp=clr(fuhrer_mi)


#map the gene names back to Nichols' id

id_geneName=id_allAttributes[,c("ids","associated_gene_names")]









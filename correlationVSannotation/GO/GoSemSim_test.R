#Ref: https://bioconductor.org/packages/devel/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html

#Ref: http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
#=> http://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html

#(!)This once happened: before updating some of the packages, the distances of some really close terms are 0

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("GOSemSim", version = "3.8")
BiocManager::install("org.Hs.eg.db", version = "3.8")
BiocManager::install("org.EcK12.eg.db", version = "3.8")
BiocManager::install("org.Sc.sgd.db", version = "3.8")

library(GOSemSim)
library(org.Hs.eg.db)
library(org.EcK12.eg.db)
library(org.Sc.sgd.db)

#Build godata for both human, E. coli K-12, yeast

##human
hsGO_MF <- godata('org.Hs.eg.db', ont="MF",computeIC = T)
hsGO_BP <- godata('org.Hs.eg.db', ont="BP",computeIC = T)
hsGO_CC <- godata('org.Hs.eg.db', ont="CC",computeIC = T)

##E. coli K-12
ECK_GO_MF=godata('org.EcK12.eg.db',ont="MF",computeIC = T)
ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = T)
ECK_GO_CC=godata('org.EcK12.eg.db',ont="CC",computeIC = T)


##yeast
sc_GO_MF=godata('org.Sc.sgd.db',ont="MF",computeIC = T)
sc_GO_BP=godata('org.Sc.sgd.db',ont="BP",computeIC = T)
sc_GO_CC=godata('org.Sc.sgd.db',ont="CC",computeIC = T)


#Test on human data
##=========================================================================
goSim("GO:0004022", "GO:0005515", semData=hsGO_MF, measure="Jiang") #0.162
goSim("GO:0004022", "GO:0005515", semData=hsGO_MF) #0.158


## GO:0006635: fatty acid beta-oxidation, GO:0019395: fatty acid oxidation
goSim("GO:0006635", "GO:0019395", semData=hsGO_BP) #Got 0.753 and afterward it was 0. Got 0.753 on 2/13/2019 again
#2 terms are close to each other. The result doesn't make sense: https://www.ebi.ac.uk/QuickGO/term/GO:0006635

#geneSim: 
geneSim("241", "251", semData=hsGO_MF, measure="Wang", combine="BMA") #0.207
##=========================================================================



#Test on E. coli K-12 data
##=========================================================================
##GO:0004022: alcohol dehydrogenase (NAD) activity, GO:0005515: protein binding
goSim("GO:0004022", "GO:0005515", semData=ECK_GO_MF) #0.158

## GO:0006635: fatty acid beta-oxidation, GO:0019395: fatty acid oxidation
goSim("GO:0006635", "GO:0019395", semData=ECK_GO_MF) #Got 0.753 and afterward it was 0. Got 0.753 on 2/13/2019 again

goSim("GO:0006635", "GO:0019395", semData=ECK_GO_BP) #Got 0.753 and afterward it was 0. Got 0.753 on 2/13/2019 again



##GO:0006635: fatty acid beta-oxidation (a child term of GO:0055114), GO:0055114: oxidation-reduction process
goSim("GO:0006635","GO:0055114",semData=ECK_GO_MF) #0 => Got 0.091 on 2/13/2019
goSim("GO:0006635", "GO:0055114", semData=ECK_GO_BP) #0 => Got 0.091 on 2/13/2019
##For other methods c("Resnik", "Lin", "Rel", "Jiang") they are all non-0

## GO:0006635: fatty acid beta-oxidation, GO:0008150: biological_process
goSim("GO:0006635", "GO:0008150", semData=ECK_GO_BP) #0.076 => 0 => Got 0.076 on 2/13/2019 again


##Similarity of a term to itself
goSim("GO:0006635", "GO:0006635", semData=ECK_GO_MF) #1


##The strain pair in Nichols with the highest pcc (0.959): 3061 3107
##3061: "GO:0022627" "GO:0005829" "GO:0005515" "GO:0048027" "GO:0070181" "GO:0003735"
##3107: "GO:0005829"

go1=c("GO:0022627","GO:0005829","GO:0005515","GO:0048027","GO:0070181","GO:0003735")
go2="GO:0005829"
mgoSim(go1, go2, semData=ECK_GO_MF, combine=NULL)
mgoSim(go1, go2, semData=ECK_GO_MF, combine="BMA") #1

#1685 2886 pcc~0.9 
go1=id_allAttributes$GO[id_allAttributes$ids==1685] %>% unique()
go2=id_allAttributes$GO[id_allAttributes$ids==2886] %>% unique() 
mgoSim(go1, go2, semData=ECK_GO_MF, combine=NULL)
mgoSim(go1, go2, semData=ECK_GO_MF, combine="BMA") #0.093 on 2/13/2019

#2515 3002 pcc~0.8
go1=id_allAttributes$GO[id_allAttributes$ids==2515] %>% unique()
go2=id_allAttributes$GO[id_allAttributes$ids==3002] %>% unique() 
mgoSim(go1, go2, semData=ECK_GO_MF, combine=NULL)
mgoSim(go1, go2, semData=ECK_GO_MF, combine="BMA") # 0.222 on 2/13/2019

#13 305 pcc~0.7
go1=id_allAttributes$GO[id_allAttributes$ids==13] %>% unique()
go2=id_allAttributes$GO[id_allAttributes$ids==305] %>% unique() 
mgoSim(go1, go2, semData=ECK_GO_MF, combine=NULL)
mgoSim(go1, go2, semData=ECK_GO_MF, combine="BMA") #0.833 on 2/13/2019

#398  2571  pcc=1.458438e-08
go1=id_allAttributes$GO[id_allAttributes$ids==398] %>% unique()
go2=id_allAttributes$GO[id_allAttributes$ids==2571] %>% unique() #no GO ID for this guy
mgoSim(go1, go2, semData=ECK_GO_MF, combine=NULL)
mgoSim(go1, go2, semData=ECK_GO_MF, combine="BMA") #0.272 on 2/13/2019


#45 548 pcc=0.01215538
go1=id_allAttributes$GO[id_allAttributes$ids==45] %>% unique()
go2=id_allAttributes$GO[id_allAttributes$ids==548] %>% unique() 
mgoSim(go1, go2, semData=ECK_GO_MF, combine=NULL)
mgoSim(go1, go2, semData=ECK_GO_MF, combine="BMA") #NA on 2/13/2019


#For Nichols': ilvC: 3061: P02358, GeneID:948723(from NCBI)   argB: 3107: P39347 (putative pseudo gene. Didn't find GeneID from NCBI)
ECK_GO_MF_symbol <- godata('org.EcK12.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)
geneSim("ilvC", "ilvC", semData=ECK_GO_MF_symbol, measure="Wang", combine="BMA") #1
geneSim("ilvC", "argB", semData=ECK_GO_MF_symbol, measure="Wang", combine="BMA") #0.289 => 0.278 on 2/13/2019

#GO annotations from our GAF
go1=id_allAttributes$GO[id_allAttributes$ids=="3061"] %>% unique
go2=id_allAttributes$GO[id_allAttributes$ids=="3107"] %>% unique
mgoSim(go1,go2, semData=ECK_GO_MF_symbol, measure="Wang", combine="BMA") #0.176 on 2/13/2019


# Got 2 EntrezID from str(ECK_GO_BP): 944744, 944748 (from NCBI search I know they are sspA, talB, respectively)
geneSim("944744", "944748", semData=ECK_GO_MF_symbol, measure="Wang", combine="BMA") #NA #It used to give me 0.18 (got NA on 1/29/2019) (still got NA on 2/13/2019)

go1=unique(id_allAttributes$GO[id_allAttributes$associated_gene_names=="sspA"])
go2=unique(id_allAttributes$GO[id_allAttributes$associated_gene_names=="talB"])
mgoSim(go1,go2, semData=ECK_GO_MF) # 0.125 => 0.272 on 2/13/2019



goSim("GO:0043229", "GO:0043231", semData=ECK_GO_CC, measure="Wang") #0.826
goSim("GO:0043229", "GO:0005622", semData=ECK_GO_CC, measure="Wang") #0.661

goSim("GO:0004022", "GO:0004024", semData=ECK_GO_MF, measure="Wang") #0.869




##=========================================================================




#Test on E. yeast data
##=========================================================================

goSim("GO:0043229", "GO:0043231", semData=sc_GO_CC, measure="Wang") #0.826 #In the paper it is 0.7727
goSim("GO:0043229", "GO:0005622", semData=sc_GO_CC, measure="Wang") #0.661 

goSim("GO:0004022", "GO:0004024", semData=sc_GO_MF, measure="Wang") #0.869 

mgoSim(c("GO:0004023","GO:0004024","GO:0046872","GO:0016491","GO:0004022","GO:0008270","GO:0004174"),
       c("GO:0046872","GO:0008270","GO:0009055","GO:0020037","GO:0005515"),
       semData=sc_GO_MF, measure="Wang") #0.51 (In the paper it is 0.693)

mgoSim(c("GO:0004023","GO:0004024","GO:0046872","GO:0016491","GO:0004022","GO:0008270","GO:0004174"),
       c("GO:0046872","GO:0008270","GO:0009055","GO:0020037","GO:0005515"),
       semData=sc_GO_MF, measure="Wang",combine=NULL)

##=========================================================================





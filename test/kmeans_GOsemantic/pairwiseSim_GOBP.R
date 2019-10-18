#calculate pairwise similarity of all GO terms (Biological Process) in Nichols' strains

#retrieve all GO terms from strains
allGO=id_allAttributes$GO %>% unique
allGO=allGO[!is.na(allGO)]

allGO=allGO[AnnotationDbi::Ontology(allGO)=="BP" & !is.na(AnnotationDbi::Ontology(allGO)=="BP")] #some terms will not belong to either BP, MF, CC so I have to remove them


#get all the pairwise combinations of terms and compute the similarity
GO1_GO2_sim=data.frame(t(combn(allGO,2)),NA)
names(GO1_GO2_sim)=c("GO1","GO2","similarity")


ECK_GO_BP=godata('org.EcK12.eg.db',ont="BP",computeIC = F)

start.time=Sys.time()

sim=pbapply(GO1_GO2_sim,1,FUN=function(Row){
  goSim(Row[1], Row[2], semData=ECK_GO_BP,measure="Wang")
})

end.time=Sys.time()
end.time-start.time #About 3 min on my PC

GO1_GO2_sim$similarity=sim

#correct the obj name to make it more precise -> save
GOBP1_GOBP2_similarity=GO1_GO2_sim

## I am surprised that no NA is observed in the similarity column:
## anyNA(GOBP1_GOBP2_similarity$similarity)


save(GOBP1_GOBP2_similarity,file="Data/GOBP1_GOBP2_similarity.RData")







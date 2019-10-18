#Examples to understand Wang method to calculate distance between gene pairs (Using 2 sets of GO annotatinos)


#Cherry picking pairs that have high, medium, low, 0 distance (Biological process)
'#############################
strain1 strain2 distance_BP
1573 3872 0.908
1618 2228 0.908
1695 2670 0.908
----------------
965 3543 0.516
969 2221 0.516
972 2559 0.516
----------------
1632 2539 0.138
1633 1823 0.138
1635 1931 0.138
----------------
1 494 0
2 2756 0
3 2757 0
'#############################


extract_GO_BP=function(strain){
  GOid=id_allAttributes$GO[id_allAttributes$ids==strain] %>% unique
  GOid=GOid[!is.na(GOid)]
  GOid=GOid[AnnotationDbi::Ontology(GOid)=="BP" & !is.na(AnnotationDbi::Ontology(GOid))]
  
  term=AnnotationDbi::Term(GOid)
  
  result=rep(NA,length(GOid))
  for(i in seq(GOid)){
    result[i]=paste(GOid[i],term[i],sep=" ")
  }
  
  return(result)
}



strain_GO_BP=list(strains=as.character(c(1573,3872,
                            1618,2228,
                            1695,2670,
                            965,3543,
                            969,2221,
                            972,2559,
                            1632,2539,
                            1633,1823,
                            1635,1931
                            )))


strain_GO_BP$BP=sapply(strain_GO_BP$strains,FUN=extract_GO_BP)

#Go get the graphs: 
##Ref: https://www.biostars.org/p/75802/







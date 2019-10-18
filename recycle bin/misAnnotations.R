##test example
#class(id_allAttributes)
#id_pwyAnnot=id_allAttributes[,c("ids","Pwy")] %>% unique



#fake_pwyAnnot_list=shuffleAnnotation(id_pwyAnnot)
#map_dbl(fake_pwyAnnot_list,function(x){x[[2]] %>% unique %>% length}) ##Double check that they contain all the ids (1 ~ 3979)

#fake_pwyAnnot_list=removeAnnotation(id_pwyAnnot)
#map_dbl(fake_pwyAnnot_list,function(x){x[[2]] %>% unique %>% length}) ##Double check that they contain all the ids (1 ~ 3979)

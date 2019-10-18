#Ref: https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html

#install.packages("VennDiagram")
library(VennDiagram)

#Example from documentation: https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf
'
venn.diagram(
  x = list(
    A = c(1:80),
    B = c(41:150),
    C = c(71:100)
  ),filename="test.tiff",
  cex = 2.5,
  cat.cex = 2.5,
  cat.dist = c(0.07, 0.07, 0.02),
  cat.pos = c(-20, 20, 20))
'  

#This give areas that are labeled 0 when two sets don't intersect -> undesirable
'
venn.diagram(
  x = list(
    A = sample(letters,15),
    B = sample(letters,10),
    C = sample(letters,5),
    D = sample(letters,5),
    E = sample(letters,5)
  ),filename="test.tiff",
  cex = 2.5,
  cat.cex = 2.5
  )
'

names(id_allAttributes)


#Here I haven't accounted strains that have 2 or more annotations from any annotation set
strain=3979

getAnnotatedID=function(attr){
  id_attr=id_allAttributes[,c("ids",attr)] %>% unique
  annotatedID=id_attr$ids[!is.na(id_attr[[paste(attr)]])] %>% unique
}


venn.diagram(
  x = list(
    Pwy=getAnnotatedID("Pwy"),
    kegg_modules=getAnnotatedID("kegg_modules"),
    pcomplex=getAnnotatedID("pcomplex"),
    regulator=getAnnotatedID("regulator"),
    operon=getAnnotatedID("operon")
  ),filename="test.tiff",
  cex = 2,
  cat.cex = 2
)
#will give areas that are labeled 0 when two sets don't intersect -> undesirable













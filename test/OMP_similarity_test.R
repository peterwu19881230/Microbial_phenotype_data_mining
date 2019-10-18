#Try to calculate term similarity in OMP

#ref (the paper): https://academic.oup.com/bioinformatics/article/33/7/1104/2843897
#ref (the software): https://cran.r-project.org/web/packages/ontologyIndex/vignettes/intro-to-ontologyX.html
packages=c("ontologyIndex","ontologySimilarity")
install.packages(packages)

for( package in packages){
  library(package, character.only = T)
}

#Get obo version of OMP from out GitHub
OMP=ontologyIndex::get_ontology("https://raw.githubusercontent.com/microbialphenotypes/OMP-ontology/master/omp.obo")
class(OMP) #"ontology_index"

#calculate ontology similarity (ref: https://cran.r-project.org/web/packages/ontologySimilarity/ontologySimilarity.pdf)




#Examples to get properties from OMP
ontologyIndex::get_term_property(OMP,"name","OMP:0000040")
ancestors_ID=ontologyIndex::get_term_property(OMP,"ancestors","OMP:0000040")

for(ID in ancestors_ID){
  cat(ontologyIndex::get_term_property(OMP,"name",ID))
  cat("\n")
}



#Examples for similarity between 2 terms: 

##OMP:0000336: beta-lactam resistance phenotype
##OMP:0006047: resistant to beta-lactam
##Our webpage: https://microbialphenotypes.org/wiki/index.php?title=Category:OMP:0006047_!_resistant_to_beta-lactam
ontologySimilarity::get_sim_grid(ontology=OMP, term_sets=list(
  "case 1"=c("OMP:0000336"),
  "case 2"=c("OMP:0006047"))) #0.8554507


##OMP:0000274: antimicrobial agent resistance phenotype 
##OMP:0000455: alkaloid chemical resistance phenotype
##Our webpage: https://microbialphenotypes.org/wiki/index.php?title=Category:OMP:0007195_!_chemical_resistance_phenotype
ontologySimilarity::get_sim_grid(ontology=OMP, term_sets=list(
  "case 1"=c("OMP:0000274"),
  "case 2"=c("OMP:0000455"))) #0.5724949


##OMP:0006024: presence of organic carbon source utilization
##OMP:0007012: presence of inorganic carbon source utilization
##Our webpage: https://microbialphenotypes.org/wiki/index.php?title=Category:OMP:0006023_!_carbon_source_utilization_phenotype
ontologySimilarity::get_sim_grid(ontology=OMP, term_sets=list(
  "case 1"=c("OMP:0006024"),
  "case 2"=c("OMP:0007012"))) #0.6670588


#An example of calculating similarity between sets of terms (using the combination of the above terms)
ontologySimilarity::get_sim_grid(ontology=OMP, term_sets=list(
  "case 1"=c("OMP:0000336","OMP:0000274","OMP:0006024"),
  "case 2"=c("OMP:0006047","OMP:0000455","OMP:0007012"))) #0.6983348










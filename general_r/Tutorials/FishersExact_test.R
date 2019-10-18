#Hypergeometric test
#Refs: 
##https://www.youtube.com/watch?v=udyAvvaMjfM
##https://rcompanion.org/rcompanion/b_07.html
##http://www.biostathandbook.com/fishers.html
##https://en.wikipedia.org/wiki/Fisher%27s_exact_test



#Example from the ref
#==================================================================
data = matrix(data=c(1,0,0,0,0,7, 5,6,8,8,5,8), nrow=6)

colnames(data) = c("handful", "bag")
rownames(data) = c("red", "yellow", "orange", "green", "brown", "blue")

data

fisher.test(data)
#==================================================================



#Example from the gene ontology handbook: chapter 13

##                             from sample   from population except the subset
##   annnotated to term t            3              1
##
##  not annnotated to term t         2              12


fisher.test(matrix(c(3,2, 1,12),ncol=2,byrow=F),alternative = "greater") 

#manual calculation:
choose(4,3)*choose(14,2)/choose(18,5)+choose(4,4)*choose(14,1)/choose(18,5)



#Example of using Nichols' co-annotation data (a little different than the above because a subset is compared to its population)
##Ref: The gene ontology handbook: chapter 13

##                    from sample   from population except the subset
##   co-annotated           1              7787-1
##
## not co-annotated         3              7914231-7787-3


fisher.test(matrix(c(1,3, 7787-1,7914231-7787-3),ncol=2,byrow=F),alternative = "greater") 

#manual calculation (p of getting 1 or more than 1 annotations out of 4 genes):
choose(7787,1)*choose(7914231-7787,3)/choose(7914231,1+3)+
choose(7787,2)*choose(7914231-7787,2)/choose(7914231,1+3)+
choose(7787,3)*choose(7914231-7787,1)/choose(7914231,1+3)+
choose(7787,4)*choose(7914231-7787,0)/choose(7914231,1+3)
  





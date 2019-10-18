#Sort Nichols' pairwise pccs in a KNN manner

#Note: to do KNN, neither diagnal elements nor upper (or lower) duplicated pccs are removed from cor_strain
dt=apply(cor_strains,2,FUN=function(strain){
  names(strain)=as.character(1:3979) 
  ##even if I don't change the type to character, it will get changed automatically once assigned
  sort(strain,decreasing = T) %>% as.data.frame
})


#Try using a cutoff: pcc=0.6
dt_0.6=lapply(dt,FUN=function(strain){
  filter=ifelse(abs(strain)>=0.6,T,F)
  strain_filtered=strain[filter,] %>% as.data.frame #I don't know why I have to do this again
  rownames(strain_filtered)=rownames(strain)[filter] #No idea why the names aren't preserved if not doing this again
  return(strain_filtered)
})

#The first example I want to do GO overrepresentation (selected arbitrarily): dt_0.6$`13`

##13  1.0000000 
##17  0.7738179
##527 0.7657747
##530 0.7489893
##453 0.7098366
##88  0.7070961
##92  0.7045078
##305 0.7004078
##457 0.6722504
##167 0.6443931

#corresponding UniProtID: P13738 P0AEM0 P28916 P75742 P77806 P0AE14 P06959 P18200 P77202 P63224

#I went through every annotation data set selection in PANTHER (http://www.pantherdb.org) but didn't get any significant result
##Curtis' example: https://docs.google.com/document/d/1shxzrMSrGUjN_aBsDhuailirkGWVNctLSPfsSJVdB2g/edit

#Distribution of no. of nearest strains 
x=sapply(dt_0.6,dim)[1,]
table(x) #The majority of them don't have that many neighbors. Only 265 (out of 3979) have them have >= 3 neighbors
hist(x,breaks=as.numeric(names(table(x))))




# Conclusion: I don't think KNN will work based on the above result



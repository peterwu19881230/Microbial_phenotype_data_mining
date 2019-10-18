##Impute Data with mean(all scores) to remove NA 
##-> This allows me to get various distances without worring to deal with NA

##For quantitative data

library(factoextra)
dist.All_Data_NAimputed.euclidean=get_dist(All_Data_NAimputed,method="euclidean") ##get_dist is an enhanced version of dist()
dist.All_Data_NAimputed.maximum=get_dist(All_Data_NAimputed,method="maximum") #also called Chebyshev distance (application in King's movement on a chess board)
dist.All_Data_NAimputed.manhattan=get_dist(All_Data_NAimputed,method="manhattan")
dist.All_Data_NAimputed.canberra=get_dist(All_Data_NAimputed,method="canberra")
dist.All_Data_NAimputed.binary=get_dist(All_Data_NAimputed,method="binary") #This distance wouldn't be useful (can look at ?dist)
dist.All_Data_NAimputed.minkowski=get_dist(All_Data_NAimputed,method="minkowski",p=3) #Have to change p to 3. Default is 2, which makes this equal to euclidean 
dist.All_Data_NAimputed.spearman=get_dist(All_Data_NAimputed,method="spearman")
#dist.All_Data_NAimputed.kendall=get_dist(All_Data_NAimputed,method="kendall") => This works with smaller data frame. But it takes long here so I haven't generated it
#=> Also by seeing this ref (http://people.revoledu.com/kardi/tutorial/Similarity/KendallDistance.html), kendall is not something that would make sense for analysing Nichols'

save(dist.All_Data_NAimputed.euclidean,
     dist.All_Data_NAimputed.maximum,
     dist.All_Data_NAimputed.manhattan,
     dist.All_Data_NAimputed.canberra,
     dist.All_Data_NAimputed.binary,
     dist.All_Data_NAimputed.minkowski,
     dist.All_Data_NAimputed.spearman,
     file="Data/other_distances.RData"
)



##For ternary data


MHD_ternary=modified_hamming_distance(Ternary_Data)
#save(MHD_ternary,file="Data/MHD_ternary")

MHD2_ternary=modified_hamming_distance_2(Ternary_Data)
#save(MHD2_ternary,file="Data/MHD2_ternary")

MHD3_ternary=modified_hamming_distance_3(Ternary_Data)
#save(MHD3_ternary,file="Data/MHD3_ternary")



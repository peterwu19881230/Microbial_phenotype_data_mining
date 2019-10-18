#Goal: Convert all data into binary (-1,0,1) using R



#Transpose for cor() to work
t_Binary_Data=t(Binary_Data)


#Calculate the correlation coefficient
cor_binary<-cor(t_Binary_Data,use="pairwise.complete.obs")

#Warning message:
#In cor(t_Binary_Data, use = "pairwise.complete.obs") :
#  the standard deviation is zero

#-> I guess the Warning message makes sense because: If the strains are not significant (all |scores|<3), standard deviation would be 0, and Pearson Correlation Coefficient used variances(sigma) in the denominator


#Write the output to a csv
##write.csv(cor_binary, file = "correlation_binary.csv")


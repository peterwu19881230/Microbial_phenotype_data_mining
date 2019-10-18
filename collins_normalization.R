#Example to get normalized scores in Nichols (I followed the material and method in Collins et al., 2006)
#Ref


#suppose I have a plate of 20 strains with 20 conditions (what about in double mutant experiement where a gene of 2 makers are crossed?)
set.seed(102)
raw_plate_mat=matrix(rnorm(20*20,mean=500),nrow=20) 
raw_plate=raw_plate_mat %>% as.numeric()

raw_plate_excludedRowCol_mat=raw_plate_mat[-c(1,20),-c(1,20)]
raw_plate_excludedRowCol=raw_plate_excludedRowCol_mat %>% as.numeric


PMM=mean(raw_plate_excludedRowCol[raw_plate_excludedRowCol>=quantile(raw_plate_excludedRowCol,probs=c(0.4,0.6))[1] & 
                                    raw_plate_excludedRowCol<=quantile(raw_plate_excludedRowCol,probs=c(0.4,0.6))[2]])

##normalize the outermost rows and columns
plate_mat=raw_plate_mat
plate_mat[1,]=plate_mat[1,]*PMM/median(plate_mat[1,])
plate_mat[20,]=plate_mat[1,]*PMM/median(plate_mat[1,])
plate_mat[,1]=plate_mat[1,]*PMM/median(plate_mat[1,])
plate_mat[,20]=plate_mat[1,]*PMM/median(plate_mat[1,])

plate=plate_mat %>% as.numeric

normalized=plate*516.1/PMM #516.1: median size of all colonies measured in the study (does "study" mean 1 plate or all screens?)





#Suppose I have 4 biological replicates of a double-mutant:
normalized_mutant=c(400, 500, 600, 400) #4 replicates from different plates. The scores have been normalized



#Following is the s-score without modifications
#========================================================================================================
var_exp=var(normalized_mutant)


n_exp=4
var_cont=5000 #cont stands for the control KAN-marked single mutant. I just faked a value here
n_cont=4 #the number of measurements of colony sizes for this control KAN-marked strain

S_var=( var_exp*(n_exp-1) + var_cont*(n_cont-1) ) / ( n_exp + n_cont - 2)

mu_exp=mean(normalized_mutant)

mu_cont=500 #cont stands for the control KAN-marked single mutant. I just faked a value here

S= (mu_exp - mu_cont) / sqrt( S_var/n_exp + S_var/n_cont )
#========================================================================================================



#Following is the s-score with modifications
#========================================================================================================
var_exp=10000 #I just faked a value here
#the maximum of the variance of normalized colony sizes for the double mutant of interest (why are there multiple variances?)
#or a minimum bound described in Collins et al


n_exp=4
var_cont=5000 #I just faked a value here
#median of the variances in normalized colony sizes observed for all double mutants containing the KAN-marked mutant of interest 
#or a minimum bound described below

n_cont=6 # The median number of experimental replicates over all the experiments

S_var=( var_exp*(n_exp-1) + var_cont*(n_cont-1) ) / ( n_exp + n_cont - 2)

mu_exp=mean(normalized_mutant)

mu_cont=600 # median of normalized colony sizes for all double mutants containing the KAN-marked mutant of interest. I just faked a value here

S= (mu_exp - mu_cont) / sqrt( S_var/n_exp + S_var/n_cont )
#========================================================================================================


#Suppose I have a set of S-scores that need to be normalized again by IQR (1.35) of Z
set.seed(102)
S=runif(1000)
hist(S)
normalized_S=( S-mean(S) )*1.35/IQR(S) 
hist(normalized_S)











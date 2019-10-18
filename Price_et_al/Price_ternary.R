#Generate ternary version of the Price's data using their suggested cutoff: For E. coli: |fitness|>0.5 and |combined t|>4

fitness=read.table("Data/fit_logratios_good.tab.txt",quote="",fill=T,sep="\t",header = T,check.names=F) 
fitness_scoresOnly=fitness[,-(1:4)]

combined_t=read.table("Data/fit_t.tab.txt",quote="",fill=T,sep="\t",header = T,check.names=F) 
combined_t_scoresOnly=combined_t[,-(1:3)]


## All column names in the fitness scores are found in combined_t scores, but not vice versa
sum(names(combined_t_scoresOnly) %in% names(fitness_scoresOnly))


combined_t_for_fitness=combined_t_scoresOnly[,names(combined_t_scoresOnly) %in% names(fitness_scoresOnly)]


significant_or_not=abs(fitness_scoresOnly)>0.5 & abs(combined_t_for_fitness)>4
significant_quantitative_Price=fitness_scoresOnly*significant_or_not


ternary_Price=significant_quantitative_Price
ternary_Price[ternary_Price>0]=1
ternary_Price[ternary_Price<0]=-1

#save(ternary_Price,file="Data/sourced/ternary_Price.RData")




#Exploratory data analysis
#==================================================================================================

sum(as.numeric(ternary_Price==1)) #5039 (no. of sinificant positve phenotypes)

sum(as.numeric(ternary_Price==-1)) #22186 (no. of sinificant negative phenotypes)

22186/5039 #this ratio is similar as the Nichols' data


(5039+22186)/(3789*162) #~ 4.4% of the total data are significant phenotypes


sum(as.numeric(Ternary_Data_324cutff_NAremoved!=0))/(3979*324) # In Nichols ~ 1.2% are significant phenotypes
#==================================================================================================







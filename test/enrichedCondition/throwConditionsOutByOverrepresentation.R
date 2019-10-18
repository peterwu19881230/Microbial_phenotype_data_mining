#Remove genes that have no significant phenotypes -> Remove conditions that don't have enriched no. of significant phenotypes



#For data before conditions are collapsed
#============================================================================================================================================
class(Ternary_Data_324cutff_NAremoved)
noOfSigPhentypes=apply(Ternary_Data_324cutff_NAremoved,2,FUN=function(col) sum(col!=0))

names(noOfSigPhentypes)=names(All_Data)

noOfSigPhentypes=sort(noOfSigPhentypes,decreasing = T)

#barplot to show no. of significant phenotypes for ranked conditions
#(!)labels are too small
#dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
#pdf(file=paste(dir_of_workingScript,"noOfSigPhentypes.pdf",sep="/"))
#barplot(noOfSigPhentypes,las=2)
#dev.off()

#table to show no. of significant phenotypes for ranked conditions
noOfSigPhentypesCondition=data.frame(condition=names(noOfSigPhentypes),noOfSigPhentypes=noOfSigPhentypes)
#============================================================================================================================================
#Conclusion: each condition has at least 2 significant phenotypes





#For collapsed conditions
#============================================================================================================================================
class(Ternary_Data_324cutff_condCollapsed) #This data shouldn't have NA because they disappear after conditions are collapsed
noOfSigPhentypesCollapsedCond=apply(Ternary_Data_324cutff_condCollapsed,2,FUN=function(col) sum(col!=0))
names(noOfSigPhentypesCollapsedCond)=names(uniqueChemIndex)
noOfSigPhentypesCollapsedCond=sort(noOfSigPhentypesCollapsedCond,decreasing = T)

#barplot to show no. of significant phenotypes for ranked conditions
#(!)labels are too small
#dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
#pdf(file=paste(dir_of_workingScript,"noOfSigPhentypesCollapsedCond.pdf",sep="/"))
#barplot(noOfSigPhentypesCollapsedCond,las=2)
#dev.off()

#table to show no. of significant phenotypes for ranked conditions
noOfSigPhentypesCollapsedCondCondition=data.frame(condition=names(noOfSigPhentypesCollapsedCond),noOfSigPhentypesCollapsedCond=noOfSigPhentypesCollapsedCond)
#============================================================================================================================================
#Conclusion: each condition has at least 8 significant phenotypes



#Find enriched conditions in terms of fitness scores
#=> Use the enriched conditions and delete the rest

#What's the rationale of using overrepresentation here? What's the difference between term enrichment and over-representation?

#Refs:
##http://blog.nextgenetics.net/?e=94
##https://en.wikipedia.org/wiki/Hypergeometric_distribution
##Here is a tool https://systems.crump.ucla.edu/hypergeometric/


##Have to discuss the following to make sure my definition is right :

#Population: 3979*324 fitness scores, whether they are significant (sucess) or nonsignificant 


# For: Ternary_Data_324cutff_NAremoved
#no. that are significant:
class(Ternary_Data_324cutff_NAremoved)
significant=apply(Ternary_Data_324cutff_NAremoved,1,FUN=function(Row){
  ifelse(sum(Row!=0)>0,T,F)
})
ternary_noSigGeneRemoved=Ternary_Data_324cutff_NAremoved[significant,]
dim(ternary_noSigGeneRemoved) #2210  324
fitness=as.numeric(Ternary_Data_324cutff_NAremoved)
success=sum(fitness!=0) #15833
total=2210*324

#perform over- or under-representation (Ref: http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html)
#Should I use over- or under- ?

#over-representation
hitInSample=noOfSigPhentypes #noOfSigPhentypes is defined above
hitInPop=15833
failInPop=2210*324-15833
sampleSize=2210
  
pVals=phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
sigPVals=pVals[pVals<=0.05]
length(sigPVals) #91 conditions to be used


# For: Ternary_Data_324cutff_condCollapsed
#no. that are significant:
class(Ternary_Data_324cutff_condCollapsed)

significant=apply(Ternary_Data_324cutff_condCollapsed,1,FUN=function(Row){
  ifelse(sum(Row!=0)>0,T,F)    
})
ternaryCondCollapsed_noSigGeneRemoved=Ternary_Data_324cutff_condCollapsed[significant,]
dim(ternaryCondCollapsed_noSigGeneRemoved)
fitness=as.numeric(Ternary_Data_324cutff_condCollapsed)
sum(fitness!=0) #11317
total=2209*114

#perform over- or under-representation (Ref: http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html)
#Should I use over- or under- ?

#over-representation
hitInSample=noOfSigPhentypesCollapsedCond #noOfSigPhentypesCollapsedCond is defined above
hitInPop=11317
failInPop=3979*114-11317
sampleSize=3979

pVals=phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
sigPVals=pVals[pVals<=0.05]
length(sigPVals) #40 conditions to be used








#functions:
#This file is by default sourced

#Function Add Columns
addCols=function(matrix,ncol){
  ##This will make the added columns to be all 0s
  matrix<-cbind(matrix,matrix(0,nrow(matrix),ncol))
  return(matrix)
}

#Convert original data to ternary. Input is the original matrix and the threshold (cutoff) for scores. Input can be a matrix or a data frame
TernaryConvert=function(matrix,thresh){
  ternary<-matrix(,nrow(matrix),ncol(matrix)) #This creates an empty matrix (Precisely, matrix with all NAs)
  rownames(ternary)<-rownames(matrix)
  colnames(ternary)<-colnames(matrix)
  ternary[matrix>-1*thresh&matrix<thresh]=0
  ternary[matrix<=-1*thresh]=-1
  ternary[matrix>=thresh]=1
  return(ternary) 
}

#Convert original data to binary. Input is the original matrix and the threshold (cutoff) for scores. Input can be a matrix or a data frame
BinaryConvert=function(matrix,thresh){
  binary<-matrix(,nrow(matrix),ncol(matrix)) #This creates an empty matrix (Precisely, matrix with all NAs)
  rownames(binary)<-rownames(matrix)
  colnames(binary)<-colnames(matrix)
  binary[matrix>=thresh]=1
  binary[matrix<thresh]=0
  return(binary) 
}


#Convert a 2D matrix to SQL type (Or I can simply use melt() in library(reshape2). cast() is the opposite function )
matrix_to_SQL=function(matrix,y_attribute="y_attribute",x_attribute="x_attribute",data_attribute="data_attribute"){
  
  SQL_data.frame<-data.frame(x_attribute=character(),y_attribute=character(),data_attribute=double())
  for(row in 1:nrow(matrix)){
    
    for(col in 1:ncol(matrix)){
      SQL_data.frame<-rbind(SQL_data.frame,
                            data.frame(                  
                              x_attribute=rownames(matrix)[row],
                              y_attribute=colnames(matrix)[col],
                              data_attribute=matrix[row,col]
                            )
      )
    }
  }
  #I have to declare names here otherwise after rbind they are going to be changed to the names of the newly bound data frames
  colnames(SQL_data.frame)<-c(y_attribute,x_attribute,data_attribute)
  return(SQL_data.frame)
}

#Input is the original matrix and a SQL dataframe (Use my function: matrix_to_SQL to convert matrix to dataframe)
matADDdf=function(matrix,df){
  for(row in 1:nrow(df)){
    #These 2 ifs are for creating new cols and rows if needed
    if(!df[row,1] %in% rownames(matrix)){
      
      #create a new row with all NAs
      new_row<-c(NA)
      #Append the new row
      matrix<-rbind(matrix,new_row)
      #Name the new row
      rownames(matrix)[nrow(matrix)]<-paste(df[row,1])
    }
    
    if(!df[row,2] %in% colnames(matrix)){
      
      #create a new row with all NAs
      new_col<-c(NA)
      #Append the new row
      matrix<-cbind(matrix,new_col)
      #Name the new row
      colnames(matrix)[ncol(matrix)]<-paste(df[row,2])
    }
    
    #If the cell is empty (NA), put the value in it
    if(is.na(matrix[paste(df[row,1]),paste(df[row,2])])){
      matrix[paste(df[row,1]),paste(df[row,2])]<-df[row,3]
    }
    #Complain if there are 2 values that have the exact same x_attribute and y_attribute
    else{
      print("Duplicated values!")
      print(df[row,])  
    }
    
  }
  return(matrix)
}

#Some narrower functions are defined as below (Previous ones are for more general use):

##Distance matrix based on pearson correlation coefficient. Correlations for rows will be used. This is what I thought Nasos has used: use="pairwise.complete.obs"(He hasn't replied for my last question)
##I use abs() because anti-correlation is better than not having any relationships
pcc_dist=function(matrix){
  ###pairwise distance would cause bias, but this is the best I can do now
  dissimilarity_pearson_pairwise<- 1 - abs(cor(t(matrix),use="pairwise.complete.obs", method="pearson"))
  distance_pcc <- as.dist(dissimilarity_pearson_pairwise)
  return(distance_pcc)
}

##Distance matrix based on Spearman correlation coefficient. Correlations for rows will be used. 
##Note: I used absolute value while in factoextra:get_dist the spearman based method doesn't
##Note: The reason why I used absolute value is the same reason as when using pcc
spearman_dist=function(matrix){
  dissimilarity_spearman_pairwise<- 1 - abs(cor(t(matrix),use="pairwise.complete.obs", method="spearman"))
  distance_spearman <- as.dist(dissimilarity_spearman_pairwise)
  return(distance_spearman)
}

##Distance matrix based on Euclidean distance
euclidean_dist=function(matrix){
  get_dist(matrix,method="euclidean")
}


#Distance matrix based on Mutual Information
my_mutual_info_dist=function(matrix){
  mi=infotheo::mutinformation(t(matrix) %>% as.data.frame) %>% natstobits #I have to convert the matrix to a dataframe. Otherwise, infotheo::mutinformation() will fail
  mi_distObj=as.dist(mi)
  mi_distance=1-mi_distObj
  return(mi_distance)
}


#This functions return the corresponding indices to the input vector "list"(smaller) 
#from the input vector "target"(larger)
##Example: list=c("A","B","C"), target=c("A","B","C","D","E") 
GetIndices=function(list,target){ 
  indices=c()
  for(i in 1:length(list)){
    indices[i]=which(target==list[i])
  }
  return(indices)
}

#This function determines the elements in a sub-cluster in an hclust obj in the nth step of merging(NoOfMerging).
elements.hclust=function(hclustering,NoOfMerging){
  node=c()  
  for(j in 1:2){
    if(hclustering$merge[NoOfMerging,j]<0){
      node=c(node,hclustering$merge[NoOfMerging,j]) ##This is created to see if 2 nodes are both singletons
      
      #Stop if both are negative. It means the algorithm reaches the end
      if(sum(node<0)==2){
        break
      }
    }else{
      node=c(node,elements.hclust(hclustering,hclustering$merge[NoOfMerging,j]))
    }
  }
  return(abs(node)) #abs is to make those singletons positive
}


#This function determines the fraction by counting No. of found genes in a length(sub-cluster)>=n sub-cluster. 
##If length(sub-cluster)>n, additional elements are ignored (Not included in the calculation of fraction)
##nodes mean the positions of elements in the original order (the original order of hclust$labels)
##When there are only 1 elements in the nodes obj fraction is meaningless. The returned value is set to NA
fractionBottomUp=function(hclustering,nodes){
  n=length(nodes)
  found=c()
  for(j in 1:length(nodes)){
    startNode=nodes[j]
    listOfOtherNodes=nodes[nodes!=startNode]
    
    if(identical(listOfOtherNodes,numeric(0))==T){return(NA)} ###NA is returned when there is only 1 element in the node obj
    
    row.first.found=which(hclustering$merge==-1*startNode,arr.ind = T)[1]
    elements=elements.hclust(hclustering,row.first.found)
    
    ##If the cluster is joined with new elements, recalculate.
    sub.cluster=row.first.found
    for(i in sub.cluster:dim(hclustering$merge)[1]){  
      ##If the name of the subcluster is found in later merging, recalculate the total number of elements after merging
      if(sub.cluster %in% hclustering$merge[i,]){
        elements=elements.hclust(hclustering,i) ##elements.hclust is a self-defined function
        sub.cluster=i ##re-define sub.cluster when new cluster is formed
      }
      
      leng=length(elements)
      if(leng>=n){break}
    }
    found[j]=sum(elements %in% listOfOtherNodes)
  }
  
  return((found+1)/n)
}

#Similar as fractionBottomUp() except it calculates the iterations until every gene in the pwy is found
##(!)This is the older version that failed to be useful
iterationBottomUp.old=function(hclustering,nodes){
  n=length(nodes)
  found=c()
  iterations=c() ##The obj that stores iterations for all starting nodes
  for(j in 1:length(nodes)){
    startNode=nodes[j]
    listOfOtherNodes=nodes[nodes!=startNode]
    
    row.first.found=which(hclustering$merge==-1*startNode,arr.ind = T)[1]
    elements=elements.hclust(hclustering,row.first.found)
    
    ##If the cluster is joined with new elements, recalculate.
    sub.cluster=row.first.found
    iteration=1 ## initialize value for No. of iterations for each starting node
    for(i in sub.cluster:dim(hclustering$merge)[1]){  
      ##If the name of the subcluster is found in later merging, recalculate the total number of elements after merging
      if(sub.cluster %in% hclustering$merge[i,]){
        elements=elements.hclust(hclustering,i) ##elements.hclust is a self-defined function
        sub.cluster=i
        iteration=iteration+1
      }
      
      found[j]=sum(elements %in% listOfOtherNodes)
      if(found[j]==n-1){
        iterations[j]=iteration
        break
      }
    }
    
  }
  return(iterations)
}

#Similar as fractionBottomUp() except: 
##1. it calculates the no. of elements until every gene in the pwy is found
##2.If the merging incorporates k genes,no. of elements=no. of elements + k.
noOfElementsBottomUp=function(hclustering,nodes){
  n=length(nodes)
  
  if(n==1){return(1)} ##If in nodes there is only 1 element, return 1
  
  found=c()
  noOfElements.all=c() ##The obj that stores no. of elements for all starting nodes
  for(j in 1:length(nodes)){
    startNode=nodes[j]
    listOfOtherNodes=nodes[nodes!=startNode]
    
    row.first.found=which(hclustering$merge==-1*startNode,arr.ind = T)[1]
    elements=elements.hclust(hclustering,row.first.found) ##elements.hclust is a self-defined function
    
    ##If the cluster is joined with new elements, recalculate.
    sub.cluster=row.first.found
    
    noOfElements=length(elements) 
    ## initialize value for no. of elements for each starting node. 
    ##It means: No. of elements that are already there after the start node is first joined
    
    
    for(i in sub.cluster:dim(hclustering$merge)[1]){  
      ##If the name of the subcluster is found in later merging, recalculate the total number of elements after merging
      if(sub.cluster %in% hclustering$merge[i,]){
        noOfElements=noOfElements+length(elements.hclust(hclustering,i))-length(elements.hclust(hclustering,sub.cluster))
        elements=elements.hclust(hclustering,i) 
        sub.cluster=i
      }
      
      found[j]=sum(elements %in% listOfOtherNodes)
      if(found[j]==n-1){
        noOfElements.all[j]=noOfElements
        print(noOfElements)
        break
      }
    }
    
  }
  return(noOfElements.all)
}

##A general function for testing bottom up tree traversing experiment. Self-defined functions are reused
##fraction, iteration (= No. of elements in the sub-cluster - 1), NoOfGenes and iteration/NoOfGenes are calculated
##Works with many to many relationship as long as the format is a 2 column table (The first column can have partial ids out of the whole clustering tree).
##---The first column contains the labels for each node.
##---The second column contains the tags to the labels of the first column. If a particular tag only appears once, it is filtered out and the fraction and so on won't be calculated.
bottomUpExp=function(hclustering,dat.lab){  
  
  all.fraction=list()
  all.elementsInCluster=list()
  NoOfElements=c()
  labs=c()
  count=0
  
  for(i in 1:length(unique(dat.lab[,2]))){
    label=unique(dat.lab[,2])[i]
    indices=dat.lab[,1][dat.lab[,2]==label] #indices are the indices of the hclust tree for the correspding tags in dat.lab
    positions=which(hclustering$labels %in% indices) # positions are the positions for hclustering$labels
    
    if(length(positions)==1){next} ## If No. of elements=1, there is no need to calculate fraction or other things
    
    count=count+1
    labs[count]=as.character(label) ## Solve the level problem if label variable has the level property
    NoOfElements[count]=length(positions)
    print(count)
    print(positions)
    all.fraction[[count]]=fractionBottomUp(hclustering,positions)
    all.elementsInCluster[[count]]=noOfElementsBottomUp(hclustering,positions) ##This is the rate limiting step
  }
  names(all.fraction)=labs
  names(all.elementsInCluster)=labs
  
  ##Calculate the average:
  avg.fraction=sapply(all.fraction,mean)
  avg.fraction.table=data.frame(avg.fraction)
  
  avg.all.elementsInCluster=sapply(all.elementsInCluster,mean)
  avg.all.elementsInCluster.table=data.frame(avg.all.elementsInCluster)
  
  result.table=cbind(rownames(avg.fraction.table),avg.fraction.table,avg.all.elementsInCluster.table,NoOfElements,NoOfElements/avg.all.elementsInCluster.table)
  if(nrow(result.table)!=0){ ##If the table is not NULL, correct the column names
    colnames(result.table)=c("Tags","Avg Fraction","Average Elements in the Cluster","No of Elements for the tag","Elements Ratio")
  }
  ##Note: additional column containing row names are bound for easy search (rownames cannot be searched in Rstudio)
  return(result.table)
}

# Throw out conditions when max(pcc)>=cutof
#The first argument is the melted table that has 3 columns: Object1 Object2 Distance(or Correlation) 
throwPcc=function(pcc.table,threshold){
  max.pcc.row=pcc.table[pcc.table[,3]==max(pcc.table[,3]),] #the row that contains the max pcc
  max.pcc=max.pcc.row[1,3]
  while(max.pcc>=threshold){
    ## sub-data that contain condition1
    subData1=pcc.table[pcc.table[,1]==max.pcc.row[1,1]|pcc.table[,2]==max.pcc.row[1,1],] #determine which one to delete 
    
    ##average pcc of strain1
    avg.pcc1=mean(subData1[,3])
    
    ## sub-data that contain condition2
    subData2=pcc.table[pcc.table[,1]==max.pcc.row[1,2]|pcc.table[,2]==max.pcc.row[1,2],] #determine which one to delete 
    
    ##average pcc of condition2
    avg.pcc2=mean(subData2[,3]); print(avg.pcc1); print(avg.pcc2); 
    
    if(avg.pcc1>avg.pcc2){pcc.table=pcc.table[pcc.table[,1]!=max.pcc.row[1,1]&pcc.table[,2]!=max.pcc.row[1,1],]
    print(paste(max.pcc.row[1,1],"with average pcc=",avg.pcc1," has been thrown out"))
    }
    else if(avg.pcc1==avg.pcc2){ ## If 2 are the same randomly throw 1 out
      i=sample(c(1,2),1)
      pcc.table=pcc.table[pcc.table[,1]!=max.pcc.row[1,i]&pcc.table[,2]!=max.pcc.row[1,i],]
      print(paste(max.pcc.row[1,i]," has been randomly selected and thrown out"))
    } 
    else{pcc.table=pcc.table[pcc.table[,1]!=max.pcc.row[1,2]&pcc.table[,2]!=max.pcc.row[1,2],]
    print(paste(max.pcc.row[1,2],"with average pcc=",avg.pcc2," has been thrown out"))
    }
    
    max.pcc.row=pcc.table[pcc.table[,3]==max(pcc.table[,3]),] #the row that contains the max pcc
    max.pcc=max.pcc.row[1,3]
  }
  return(pcc.table)
}


#Based on throwPcc() but only gives the name of condition having been thrown out each time 
throwPcc_sequence=function(pcc.table){
  namesThrownOut=c(); count=1
  
  max.pcc.row=pcc.table[pcc.table[,3]==max(pcc.table[,3]),] #the row that contains the max pcc
  max.pcc=max.pcc.row[1,3]
  
  while(dim(pcc.table)[1]>=2){
    ## sub-data that contain condition1
    subData1=pcc.table[pcc.table[,1]==max.pcc.row[1,1]|pcc.table[,2]==max.pcc.row[1,1],] #determine which one to delete 
    
    ##average pcc of strain1
    avg.pcc1=mean(subData1[,3])
    
    ## sub-data that contain condition2
    subData2=pcc.table[pcc.table[,1]==max.pcc.row[1,2]|pcc.table[,2]==max.pcc.row[1,2],] #determine which one to delete 
    
    ##average pcc of condition2
    avg.pcc2=mean(subData2[,3]); print(avg.pcc1); print(avg.pcc2); 
    
    if(avg.pcc1>avg.pcc2){pcc.table=pcc.table[pcc.table[,1]!=max.pcc.row[1,1]&pcc.table[,2]!=max.pcc.row[1,1],]
    print(paste(max.pcc.row[1,1],"with average pcc=",avg.pcc1," has been thrown out"))
    namesThrownOut[count]=max.pcc.row[1,1]; count=count+1
    }
    else if(avg.pcc1==avg.pcc2){ ## If 2 are the same randomly throw 1 out
      i=sample(c(1,2),1)
      pcc.table=pcc.table[pcc.table[,1]!=max.pcc.row[1,i]&pcc.table[,2]!=max.pcc.row[1,i],]
      print(paste(max.pcc.row[1,i]," has been randomly selected and thrown out"))
      namesThrownOut[count]=max.pcc.row[1,i]; count=count+1
    } 
    else{pcc.table=pcc.table[pcc.table[,1]!=max.pcc.row[1,2]&pcc.table[,2]!=max.pcc.row[1,2],]
    print(paste(max.pcc.row[1,2],"with average pcc=",avg.pcc2," has been thrown out"))
    namesThrownOut[count]=max.pcc.row[1,2]; count=count+1
    }
    
    max.pcc.row=pcc.table[pcc.table[,3]==max(pcc.table[,3]),] #the row that contains the max pcc
    max.pcc=max.pcc.row[1,3]
  }
  return(namesThrownOut)
}

#This function takes a distance object as input and output a melted dataframe
#!(NA will be kicked out)
melt_dist=function(distance){
  m1=as.matrix(distance)
  m1[upper.tri(m1)]=NA #I suppose that no distance can be NA, so I can use this to do filtering
  diag(m1)=NA
  library(reshape2)
  m2=melt(m1) #I suppose that no distance can be NA, so I can use this to do filtering
  if(class(m2$Var1)!="character") m2$Var1=as.character(m2$Var1) #This is to prevent numeric names (the class should still be "character") being converted to "numeric" by melt()
  if(class(m2$Var2)!="character") m2$Var2=as.character(m2$Var2) #This is to prevent numeric names (the class should still be "character") being converted to "numeric" by melt()
  
  m2=m2[!is.na(m2[,3]) | is.nan(m2[,3]),] # "| is.nan(m2[,3])" is used because I want to keep the NaN values
  
  m2=m2[,c(2,1,3)] ##reorder the columns
  names(m2)=c("object_1","object_2","value") ##name the columns
  return(m2)
}


#This function takes a distance object as input and output a melted and sorted dataframe
#!(NA will be kicked out)
meltANDsort_dist=function(distance,decreasing=F){
  m2=melt_dist(distance) #melt_dist is a self-defined function
  m2=m2[order(m2[,3],decreasing=decreasing),] ##reorder by distance
  return(m2)
}


#give a molten matrix by inputting a pairwise similarity/distance matrix
#Note: this function is the same as melt_dist except the input can be a matrix
melt_pairwiseMatrix=function(matrix_){
  molten=melt_dist(as.dist(matrix_))
  return(molten)
}


#Can test the function using the following:
'
set.seed(101)
mat=matrix(sample(c(0,1,-1),20,replace = T),4,5)
'

#Performance test:
'
set.seed(1)
mat=matrix(sample(c(0,1,-1),3979*324,replace = T),3979,324)

start.time=Sys.time()
modified_hamming_distance(mat)
end.time=Sys.time()
end.time-start.time #Time difference of 1.081457 mins
'

#Calculate the ratio of T in a TF vector
ratio=function(TF){
  if(sum(TF %in% c(1,0))!=length(TF)) stop("This is not a vector that contains only T and F") 
  TF=as.numeric(TF)
  sum(TF)/length(TF)
}


#Calculate hamming distance
hamming_distance = function(X) {
  uniqs <- unique(as.vector(X))
  U <- X == uniqs[1]
  H <- U %*% t(U)
  for ( uniq in uniqs[-1] ) {
    U <- X == uniq
    H <- H + U %*% t(U)
  }
  ncol(X) - H
}



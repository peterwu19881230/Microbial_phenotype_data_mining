##====Package:tableSMY====


#Check if the matrix, data frame or data table contains any NA, NAN, NULL or ""
anyIncomplete=function(table){
  if(!(class(table) %in% c("matrix","data.frame","data.table"))) stop("Input is not a matrix, data frame or data table")
  
  out=list()
  out$dimension=paste("Dimension: ",dim(table)[1]," rows * ",dim(table)[2]," columns",sep="")
  
  #Check NA
  na=apply(table,2,FUN=function(column){
    any=sum(is.na(column))
    return(any)
  })
  out$na=na[na>=1]
  
  #Check NULL
  null=apply(table,2,FUN=function(column){
    any=sum(is.null(column))
    return(any)
  }
  )
  out$null=null[null>=1]
  
  #Check NaN
  nan=apply(table,2,FUN=function(column){
    any=sum(is.nan(column)) #is.nan works fine for this because the columns called are vectors (if I call rows they are gonna be data frames)
    return(any)
  }
  )
  out$nan=nan[nan>=1]
  
  #Check ""
  empty=apply(table,2,FUN=function(column){
    any=sum(column=="",na.rm=T)
    return(any)
  }
  )
  out$empty=empty[empty>=1]
  
  total=sum(c(na,null,nan,empty))
  out$completeness=paste(total,' of NA, NAN, NULL, or empty character is found from ',dim(table)[1]*dim(table)[2], ' data points. They constitute ',total/(dim(table)[1]*dim(table)[2])*100,'%',sep="")
  
  return(out)
}


#Filter the matrix, data frame or data table so that all rows/columns with NA/NAN/NULL/"" are removed
filterTable=function(table){
  if(!(class(table) %in% c("matrix","data.frame","data.table"))) stop("Input is not a matrix, data frame or data table")
  
  
  colToBeRemoved=apply(table,2,FUN=function(col){
    ifelse(sum(is.na(col)+is.nan(col)+is.null(col)+sapply(col,FUN=function(point){identical(as.character(point),"")}))>=1,T,F) 
  })
  #Note: I don't have to do the same thing for rows because all the undesired values are removed after doing columns
  
  new_table=table[,!colToBeRemoved]
  return(new_table)
}

#A heatmap to visualize where NA/NAN/NULL/"" are located (Note: I am using the pheatmap package, and the speed can be terrible when large dataset is applied)
##Ref: https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps/15506652


#add threshold so that it can be set to plot the data above that value
graphTable=function(table){
  if(!(class(table) %in% c("matrix","data.frame","data.table"))) stop("Input is not a matrix, data frame or data table")
  
  #summarize the table if it contains only numbers
  if(class(table)=="matrix" && typeof(table) %in% c("integer","single","double")){
    table[]=sapply(table,FUN=function(point)as.numeric(point)) # [] is to preserve the structure of the matrix/data frame/data table
    print("From all data in the table"); print(summary(as.numeric(table))) #Note: most numbers (including 1,2,3...) are represented as double by default in R
  }else table[]=sapply(table,FUN=function(point)as.numeric(point)) 
  
  #Remove the names so they don't get plotted
  colnames(table)=NULL
  rownames(table)=NULL
  
  
  #Check to see if required packages have to be installed
  #Ref: https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
  listOfPackages=c("grid", "pheatmap")
  new_pack=listOfPackages[!(listOfPackages %in% installed.packages()[,"Package"])]
  if(length(new_pack)){
    print(paste(new_pack," is required and being installed...",sep=""))
    install.packages(new_pack)
  } 
  
  
  #Reset the range for colors: middle color: median, color range: lower-upper 
  #(This way the breaks below median and above median are not equal, but it should be fine)
  q=quantile(table,na.rm=T); q1=q["25%"]; med=q["50%"]; q3=q["75%"]; iqr=q3-q1
  lower=q1-3*iqr
  upper=q3+3*iqr
  palette.breaks=c(0:50*(med-lower)/50+lower,1:50*(upper-med)/50+med) #Sometimes (eg. when inputting ternary data) there will be this issue: 'breaks' are not unique. Have to fix
  
  ##print it
  pheatmap::pheatmap(table,cluster_rows=F,cluster_cols=F,breaks=palette.breaks,show_rownames =F, show_colnames = F)
}

#Check if an inpute vector has duplicates. If so, return the frequency table
checkDuplicates_vect=function(vect){
  if(sum(duplicated(vect))>=1){
    print("Some duplicates are found:")
    table(vect)                
  }else{return("Everything in this vector is unique")}
}



#change rownames or colnames of a matrix/dataframe/datatable based on another matrix/dataframe/datatable

#nameForTable is a 2 column matrix/dataframe/datatable: 1st column/row: original names; 2nd column/row: new names
changeNames=function(rowOrCol,Table,nameForTable){
  
  if(rowOrCol=="row"){
    IndexForNewName=match(rownames(Table),nameForTable[,1])
    rownames(Table)=nameForTable[,2][IndexForNewName]
  }else if(rowOrCol=="col"){
    IndexForNewName=match(colnames(Table),nameForTable[,1])
    colnames(Table)=nameForTable[,2][IndexForNewName]
  }else{
    stop("Enter either \"row\" or \"col\"")  
  } 
  
  return(Table)
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


#Convert every element in a matrix or dataframe to another type (eg. integer -> character)
##The reason I wrote this is because things such as as.numeric() and as.character... 
##...convert a df or matrix to a vector and doesn't preserve the original data structure
convertMatOrDF=function(dat,function_new_type){
  if(class(dat)=="matrix"|class(dat)=="data.frame"){
    for(i in 1:dim(dat)[2]){
      dat[,i]=function_new_type(dat[,i])
    }
  }else print("The input data is not a matrix or dataframe")
  
  return(dat)
}



#Filter out NA from a vector
removeNA=function(vec){
  vec[!is.na(vec)]
}






#Test code

# mat=matrix(c(1,2,3,4,5,6),ncol=2)
# graphTable(mat)
# 
# 
# set.seed(101)
# random.matrix=matrix(runif(500, min = -1, max = 1), nrow = 50)
# graphTable(random.matrix)
# 
# set.seed(101)
# random.matrix[sample(1:50,10),sample(1:10,2)]=NA
# graphTable(random.matrix)
# 
# anyIncomplete(random.matrix)
# 
# filtered_random.matrix=filterTable(random.matrix)
# str(filtered_random.matrix)
# 
# checkDuplicates_vect(c(1,1,2,3,4,4,4,5,6,7,8,9,10))
# 
# 
# 
# Table=matrix(rnorm(2*3),ncol=2,nrow=3)
# rownames(Table)=c("one","two","three")
# colnames(Table)=c("col_one","col_two")
# Table
# 
# 
# rowNameForTable=matrix(c("two","one","three","TWO","ONE","THREE"),ncol=2,byrow=FALSE)
# colNameForTable=matrix(c("col_two","col_one","COL_TWO","COL_ONE"),ncol=2,byrow=FALSE)
# 
# #newTable=changeNames(rowOrCol="test",Table,nameForTable) #test the error message of the function
# newTable=changeNames(rowOrCol="row",Table,rowNameForTable) #test rownames
# newTable=changeNames(rowOrCol="col",Table,colNameForTable) #test colnames 


# Table=matrix(rnorm(2*3),ncol=3,nrow=4)
# Table
# 
# correlation_dist_Table=as.dist(cor(t(Table)))
# 
# melt_dist(correlation_dist_Table)
# 
# meltANDsort_dist(correlation_dist_Table)
# 
# melt_pairwiseMatrix(cor(t(Table)))



# x=data.frame(a=c(1,2,3),b=c(4,5,6),c=c(1.5,2.5,3.5))
# str(x)
# y=convertMatOrDF(x,as.character)
# str(y)




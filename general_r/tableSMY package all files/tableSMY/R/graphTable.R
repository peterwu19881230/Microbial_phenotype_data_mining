graphTable <-
function(table){
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

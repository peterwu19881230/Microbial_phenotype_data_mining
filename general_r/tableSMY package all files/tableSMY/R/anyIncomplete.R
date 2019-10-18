anyIncomplete <-
function(table){
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

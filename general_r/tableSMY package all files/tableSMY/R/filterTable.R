filterTable <-
function(table){
  if(!(class(table) %in% c("matrix","data.frame","data.table"))) stop("Input is not a matrix, data frame or data table")
  
  
  colToBeRemoved=apply(table,2,FUN=function(col){
    ifelse(sum(is.na(col)+is.nan(col)+is.null(col)+sapply(col,FUN=function(point){identical(as.character(point),"")}))>=1,T,F) 
  })
  #Note: I don't have to do the same thing for rows because all the undesired values are removed after doing columns
  
  new_table=table[,!colToBeRemoved]
  return(new_table)
}

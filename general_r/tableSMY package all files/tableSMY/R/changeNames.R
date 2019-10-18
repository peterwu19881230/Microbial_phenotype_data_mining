changeNames <-
function(rowOrCol,Table,nameForTable){
  
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

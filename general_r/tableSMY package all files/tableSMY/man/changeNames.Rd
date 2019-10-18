\name{changeNames}
\alias{changeNames}
\title{
Change rol/col names
}
\description{
Change rownames or colnames of a matrix/dataframe/datatable based on another matrix/dataframe/datatable
}
\usage{
changeNames(rowOrCol, Table, nameForTable)
}
\arguments{
  \item{rowOrCol}{
Whether it's the row or col names that need to be changed
}
  \item{Table}{
A matrix/dataframe/datatable
}
  \item{nameForTable}{
A 2 column matrix/dataframe/datatable: 1st column/row: original names; 2nd column/row: new names 
}
}
\author{
Peter I-Fan Wu
}

\examples{
Table=matrix(rnorm(2*3),ncol=2,nrow=3)
  rownames(Table)=c("one","two","three")
  colnames(Table)=c("col_one","col_two")
  Table
  
  
  rowNameForTable=matrix(c("two","one","three","TWO","ONE","THREE"),ncol=2,byrow=FALSE)
  colNameForTable=matrix(c("col_two","col_one","COL_TWO","COL_ONE"),ncol=2,byrow=FALSE)
  
  #newTable=changeNames(rowOrCol="test",Table,nameForTable) #test the error message of the function
  newTable=changeNames(rowOrCol="row",Table,rowNameForTable) #test rownames
  newTable=changeNames(rowOrCol="col",Table,colNameForTable) #test colnames 
}


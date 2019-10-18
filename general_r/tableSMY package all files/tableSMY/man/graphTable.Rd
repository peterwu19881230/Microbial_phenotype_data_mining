\name{graphTable}
\alias{graphTable}
\title{
Draw a heat map of your matrix/dataframe/datatable
}
\description{
Draw a heat map that gives intuitive overview of your matrix/dataframe/datatable.
}
\usage{
graphTable(table)
}
\arguments{
  \item{table}{
A matrix/dataframe/datatable
}
}
\author{
Peter I-Fan Wu
}
\examples{
mat=matrix(c(1,2,3,4,5,6),ncol=2)
graphTable(mat)
  
set.seed(101)
random.matrix=matrix(runif(500, min = -1, max = 1), nrow = 50)
graphTable(random.matrix)
  
set.seed(101)
random.matrix[sample(1:50,10),sample(1:10,2)]=NA
graphTable(random.matrix)

}
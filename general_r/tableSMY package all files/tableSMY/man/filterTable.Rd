\name{filterTable}
\alias{filterTable}
\title{
Generate quick visualization of your matrix/dataframe and filter any NA/NULL/""
}
\description{
This package aims to provide neat small functions that can be used to clean and visualize tabular data. 
}
\usage{
filterTable(table)
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
set.seed(101)
random.matrix=matrix(runif(500, min = -1, max = 1), nrow = 50)

set.seed(101)
random.matrix[sample(1:50,10),sample(1:10,2)]=NA

filtered_random.matrix=filterTable(random.matrix)
str(filtered_random.matrix)

}
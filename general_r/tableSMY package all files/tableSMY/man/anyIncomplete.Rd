\name{anyIncomplete}
\alias{anyIncomplete}
\title{
Check Incompletion
}
\description{
Check if your matrix/dataframe/datatable has any Incompletion (NA, NULL, NaN, "") and return the statistics of them
}
\usage{
anyIncomplete(table)
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
graphTable(random.matrix)

set.seed(101)
random.matrix[sample(1:50,10),sample(1:10,2)]=NA
graphTable(random.matrix)

anyIncomplete(random.matrix)
}
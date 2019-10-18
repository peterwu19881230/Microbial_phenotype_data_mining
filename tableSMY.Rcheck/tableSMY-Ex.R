pkgname <- "tableSMY"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "tableSMY-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('tableSMY')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("anyIncomplete")
### * anyIncomplete

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: anyIncomplete
### Title: Check Incompletion
### Aliases: anyIncomplete

### ** Examples

set.seed(101)
random.matrix=matrix(runif(500, min = -1, max = 1), nrow = 50)
graphTable(random.matrix)

set.seed(101)
random.matrix[sample(1:50,10),sample(1:10,2)]=NA
graphTable(random.matrix)

anyIncomplete(random.matrix)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("anyIncomplete", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("changeNames")
### * changeNames

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: changeNames
### Title: Change rol/col names
### Aliases: changeNames

### ** Examples

Table=matrix(rnorm(2*3),ncol=2,nrow=3)
  rownames(Table)=c("one","two","three")
  colnames(Table)=c("col_one","col_two")
  Table
  
  
  rowNameForTable=matrix(c("two","one","three","TWO","ONE","THREE"),ncol=2,byrow=FALSE)
  colNameForTable=matrix(c("col_two","col_one","COL_TWO","COL_ONE"),ncol=2,byrow=FALSE)
  
  #newTable=changeNames(rowOrCol="test",Table,nameForTable) #test the error message of the function
  newTable=changeNames(rowOrCol="row",Table,rowNameForTable) #test rownames
  newTable=changeNames(rowOrCol="col",Table,colNameForTable) #test colnames 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("changeNames", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("checkDuplicates_vect")
### * checkDuplicates_vect

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: checkDuplicates_vect
### Title: Check items that occur more than once
### Aliases: checkDuplicates_vect

### ** Examples

checkDuplicates_vect(c(1,1,2,3,4,4,4,5,6,7,8,9,10))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("checkDuplicates_vect", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("filterTable")
### * filterTable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: filterTable
### Title: Generate quick visualization of your matrix/dataframe and filter
###   any NA/NULL/""
### Aliases: filterTable

### ** Examples

set.seed(101)
random.matrix=matrix(runif(500, min = -1, max = 1), nrow = 50)

set.seed(101)
random.matrix[sample(1:50,10),sample(1:10,2)]=NA

filtered_random.matrix=filterTable(random.matrix)
str(filtered_random.matrix)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("filterTable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("graphTable")
### * graphTable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: graphTable
### Title: Draw a heat map of your matrix/dataframe/datatable
### Aliases: graphTable

### ** Examples

mat=matrix(c(1,2,3,4,5,6),ncol=2)
graphTable(mat)
  
set.seed(101)
random.matrix=matrix(runif(500, min = -1, max = 1), nrow = 50)
graphTable(random.matrix)
  
set.seed(101)
random.matrix[sample(1:50,10),sample(1:10,2)]=NA
graphTable(random.matrix)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("graphTable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

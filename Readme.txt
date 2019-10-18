#Microbial phenotype data mining is one of the core projects to develop the statistical application of OMP
#Nichols' Data mining is used as a gold standard

*All txt, tsv, csv, .RData... are in the Data folder. Commonly used objects are in the sourced folder. The corresponding .rdx and .rdb files (for lazy loading) will be created from the .RData file by Nichols_preload.R if they don't exist yet.
**If a .RData file needs to be remade, the corresponding .rdx and .rdb files should also be regenerated


*Nichols_preload.R is the file that loads some important stuff for the project (Dataset , functions and important objects I created)

*usedLibraries.R loads the libraries I use often
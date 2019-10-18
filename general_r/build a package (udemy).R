#1. Make sure functions are sourced

#2. create the package folder with files in the getwd() directory (the directory Rstuio has been set to)
package.skeleton(list=c("anyIncomplete","filterTable","graphTable","checkDuplicates_vect","changeNames","relate"),"tableSMY")

#3.
Sys.which("pdflatex")

#Edit this text file: DESCRIPTION (For the format of the line "description" see this: https://stackoverflow.com/questions/34955247/error-installing-a-r-package)
#Edit every file under this folder: man

#open terminal and run:
##R CMD build (directory name of the package)   ##(Package will be built based on the current directory. Specifying a second path doesn't work)
##R CMD check  --as-cran (directory name of the package)

#Read and delte this file: Read-and-delete-me

#Submit the built .tar.gz to CRAN: https://cran.r-project.org/submit.html

listOfPackages=c("googlesheets", "tidyverse")
new_pack=listOfPackages[!(listOfPackages %in% installed.packages()[,"Package"])]
if(length(new_pack)) install.packages(new_pack)


library(googlesheets); library(tidyverse)

gs_ls() #This will open a browser and need user's authorization 
sheet=gs_title("Nichols_FinalListOfStrainsForMakingPagesInOMPwiki") #Read the sheet identified by the name
gs_ws_ls(sheet)

strainInfo=gs_read(sheet,ws="2018-07-14_ReadyToMakeStrainPages")

dim(strainInfo)

#Create and play with a new sheet (https://cran.r-project.org/web/packages/googlesheets/googlesheets.pdf)
yo <- gs_new("yo") #create a new sheet
yo <- gs_edit_cells(yo, input = head(iris), trim = TRUE)
gs_read(yo)

yo <- gs_ws_new(yo, ws = "byrow_FALSE")
yo <- gs_edit_cells(yo, ws = "byrow_FALSE",
                    input = LETTERS[1:5], anchor = "A8")
gs_read_cellfeed(yo, ws = "byrow_FALSE", range = "A8:A12") %>% gs_simplify_cellfeed()
yo <- gs_ws_new(yo, ws = "byrow_TRUE")
yo <- gs_edit_cells(yo, ws = "byrow_TRUE", input = LETTERS[1:5],
                    anchor = "A8", byrow = TRUE)
gs_read_cellfeed(yo, ws = "byrow_TRUE", range = "A8:E8") %>% gs_simplify_cellfeed()
yo <- gs_ws_new(yo, ws = "col_names_FALSE")
yo <- gs_edit_cells(yo, ws = "col_names_FALSE", input = head(iris),
                    trim = TRUE, col_names = FALSE)
gs_read_cellfeed(yo, ws = "col_names_FALSE") %>% gs_reshape_cellfeed(col_names = FALSE)

#gs_delete(yo) #this is to delete the entire sheet


#Create and play with a new sheet (another try)
test=gs_new("test") #create a new sheet
dataToPut=gs_edit_cells(test, ws = "Sheet1",input = data.frame(ECK=strainInfo$`ECK ID`[1:7]), anchor = "C1")






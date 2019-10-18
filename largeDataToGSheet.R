listOfPackages=c("googlesheets", "tidyverse")
new_pack=listOfPackages[!(listOfPackages %in% installed.packages()[,"Package"])]
if(length(new_pack)) install.packages(new_pack)


library(googlesheets); library(tidyverse)

MHD4=pwy_NAimputed.listOfResults[[4]]

gs_ls() #This will open a browser and need user's authorization 
MHD4_tab=gs_new("MHD4_tab") #create a new sheet

##URL for created sheet: https://docs.google.com/spreadsheets/d/1StYLM55Tiv3ZT7ZQP3wWrM1FqzmSBcfORXzfbEkOhuo/edit#gid=0

gs_edit_cells(MHD4_tab,input = MHD4[1:10000,]) #This takes very long (more than an hr) and I stopped it.
#(!) When uploading the whole table it gives me this: 
#Error in function_list[[k]](value) : Bad Request (HTTP 400).
#scrape KEGG modules + map to Nichols' ids
#Refs:
#https://www.analyticsvidhya.com/blog/2017/03/beginners-guide-on-web-scraping-in-r-using-rvest-with-hands-on-knowledge/
#http://bradleyboehmke.github.io/2015/12/scraping-html-text.html

library(rvest)
library(stringr)
library(dplyr)

## Get url for all Nichols' strains (genes)
baseURL="https://www.kegg.jp/dbget-bin/www_bget?eco:"
ids_bNumber=ECK_1st_table[,c("ids","bNumber")] %>% distinct #This is for mapping to Nichols' ids
  
## Start scraping for modules of each gene
start.time=Sys.time()

modules=list()
counter=1
for(bNumber in ids_bNumber$bNumber){
  url=paste(baseURL,bNumber,sep="")
  webpage=read_html(url)
  table_text=webpage %>% html_nodes("table") %>% html_text() #Modules can be obtained from tables
  info=str_extract(table_text,"^eco_M[0-9]{1,6}(.*)") # (.*) captures for anything that follows
  modules[[counter]]=info[!is.na(info)]
  counter=counter+1
}  

end.time=Sys.time()
end.time-start.time #Time difference of 1.512027 hours

names(modules)=ids_bNumber$ids #This is mapping to Nichols' ids
KEGGmodulesForNichols=modules
#save(KEGGmodulesForNichols,file="Data/sourced/KEGGmodulesForNichols.RData")


#How many modules are pathways? (Some of them are protein complexes)
modules=unlist(KEGGmodulesForNichols) %>% unique #total No. = 196




















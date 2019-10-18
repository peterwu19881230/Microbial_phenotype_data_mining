#useful URLs
#http://yetanothermathprogrammingconsultant.blogspot.com/2016/02/r-lazy-load-db-files.html
#https://stackoverflow.com/questions/21583382/r-how-to-lazyload-variables-from-inst-extdata-in-r-package
#http://adv-r.had.co.nz/Environments.html

e_test=new.env(parent=globalenv())
e_test$X = c(1,2,3)
e_test$Y = c(4,5,6)


tools:::makeLazyLoadDB(e_test,"test_DBNAME") #This makes the 2 database files (.rdb and .rdx)
lazyLoad("test_DBNAME")

#save.image(file="testEnv.RData") #just a short-cut for ‘save my current workspace’
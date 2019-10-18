#Goal: instead of using php to upload all the .csv, maybe I can just use RMySQL?

library(RMySQL)


#I used this to change MySQL password after loging inro MySQL: 
#ALTER USER 'root'@'localhost' IDENTIFIED BY 'MyNewPass';

#Have to change the password. Dangerous
#Connect DB 
mydb = dbConnect(MySQL(), user='root', password='1234', dbname='chemgen', host='localhost')

#Try to create a new table in the DB
dbWriteTable(mydb, name='testTABLE', value=mtcars)


#Get data from SQL
rs = dbSendQuery(mydb, "select * from testTABLE")
#transform the data into tidy format (Will be data.frame but instead of having row names, the row names will be in the first column)
data = fetch(rs, n=-1)

#delete the table
dbSendQuery(mydb, "DROP TABLE IF EXISTS testTABLE")

#What will dbWriteTable() do if the input is a complicated matrix like Nichols' Strains-Conditions data?
dbWriteTable(mydb, name='testTABLE', value=All_Data)
#=> All conditions will be separate variables

#Use my function matrix_to_SQL to transform part of All_Data into the SQL table format. Estimated run time for all All_Data is 22 min (haven't really tried)
SQL_partial_Data<-matrix_to_SQL(All_Data[1:100,1:100],"Strain","Condition","Growth Score")
dbWriteTable(mydb, name='testTABLE', value=SQL_partial_Data)


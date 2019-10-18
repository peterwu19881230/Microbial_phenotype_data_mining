#Goal: generate files to run 39 jobs on Ada
#0. set the path 
#1. Generate .R files
#2. Generate .lsf files
#3. Generate a txt file that contains all "bsub < a_job.lsf"

#set the path
path="clustering/hamming distanceNo0_corrected on ada/"


#Read the template file => Assign 1~39 for "part" variable and output 39 .R files
conn=file(paste(path,"/hamming_No0_corrected_template.R",sep=""))
txt <-readLines(conn)
close(conn)
#How to solve the warning "incomplete final line found on '......'":
#https://stackoverflow.com/questions/5990654/incomplete-final-line-warning-when-trying-to-read-a-csv-file-into-r


for(i in 1:39){
  txt.out=sub("part=",paste("part=",i,sep=""),txt) #assign 1~39 for "part" variable
  
  conn<-file(paste(path,"hamming_No0_corrected_part_",i,".R",sep=""))
  writeLines(txt.out, conn)
  close(conn)
}


#Read the template. lsf file => assign job No.1~39 and output 39 .lsf files
conn=file(paste(path,"/hamming_No0_corrected_template.lsf",sep=""))
txt <-readLines(conn)
close(conn)
#How to solve the warning "incomplete final line found on '......'":
#https://stackoverflow.com/questions/5990654/incomplete-final-line-warning-when-trying-to-read-a-csv-file-into-r


for(i in 1:39){
  txt.out=sub("Rscript",paste("Rscript hamming_No0_corrected_part_",i,".R",sep=""),txt) #assign 1~39 for "part" variable
  conn<-file(paste(path,"hamming_No0_corrected_part_",i,".lsf",sep=""))
  writeLines(txt.out, conn)
  close(conn)
}

#Generate a txt file that contains all "bsub < a_job.lsf"
txt.out=c()
for(i in 1:39){
  txt.out=c(txt.out, paste("bsub < hamming_No0_corrected_part_",i,".lsf",sep="")) ##\n is not needed here. (But I don't know why)
  
}

conn<-file(paste(path,"hamming_No0_corrected_allJobs.txt"))
writeLines(txt.out, conn)
close(conn)









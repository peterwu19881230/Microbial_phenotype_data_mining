#improt the dataset and convert types of columns
library(readxl)
Norm_NR <- read_excel("Documents/SBMI courses/BMI 5330/Final project/Norm_NR.xlsx", 
                        col_types = c("numeric", "text", "text", 
                                      "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric"))

#checl the class od eahc column
sapply(Norm_NR, class)

#get the mean of control columns (CDA) & exp columns (PDA)
Norm_NR$meanCDA <- apply(Norm_NR[,4:16], 1, mean)
Norm_NR$meanPDA <- apply(Norm_NR[,17:38], 1, mean)

#count how many NAs of each row
apply(Norm_NR, 2, function(x) length(which(!is.na(x))))

#remove NAs
Norm_NR_narm <- Norm_NR[complete.cases(Norm_NR), ]

#grab gene name, ID, and two means (CDA & PDA) into a dataframe
library("dplyr")
df <- select(Norm_NR_narm, 2, 3, meanCDA, meanPDA)

#convert rows into columns and columns into rows
final_df <- as.data.frame(t(df))



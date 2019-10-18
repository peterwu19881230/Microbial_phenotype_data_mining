#test whether ifelse works on matrix and data frame


##4 X 6 matrix
mat_small = matrix(c(1,2,4,5,6,3,5,7,0,0,3,5),nrow=4,ncol=3)
mat_small

greaterThanTwo_mat=ifelse(mat_small>2,1,0)
greaterThanTwo_mat


## 4 X 6 dataframe
df_small=as.data.frame(mat_small)
df_small

greaterThanTwo_df=ifelse(df_small>2,1,0)
greaterThanTwo_df

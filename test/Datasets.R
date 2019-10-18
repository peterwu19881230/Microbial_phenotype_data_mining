#Toy data for experiments:
#This file is by default sourced

##Phenotypes of synthetic genes a,b,c. They have the same number of 0s and 1s
a<-c(1,1,0,0,0)
b<-c(0,0,0,1,1)
c<-c(1,0,1,0,0)
a_b_c<-rbind(a,b,c)
colnames(a_b_c)<-c("NaCl","glucose","Ampicillin","Tabasco","UV")


##Random matrices that will help me experimenting:
##12 X 10
mat = cbind(rbind(matrix(rnorm(16, -1), 4), matrix(rnorm(32, 1), 8)),
            rbind(matrix(rnorm(24, 1), 4), matrix(rnorm(48, -1), 8)))

##4 X 6
mat_small<-matrix(c(1,2,4,5,6,3,5,7,0,0,3,5),nrow=4,ncol=3)

##5 X 5
mat55<-matrix(c(0,0,0,0,0,
                0,0,0,0,0,
                1,-1,0,0,1,
                1,1,1,0,-1,
                -1,1,-1,0,0
),nrow=5,ncol=5)

# A synthetic binary data (12 X 10)
synthetic_bibary<-matrix(
c(
c(1,1,1,1,1,1,1,1,1,1),
c(0,0,0,0,0,0,0,0,0,0),
c(1,1,1,1,1,0,0,0,0,0),
c(1,0,1,0,1,0,1,0,1,0),
c(1,1,1,0,0,0,0,0,0,0),
c(1,1,1,1,1,1,1,0,0,0),
c(1,1,0,0,1,1,0,0,1,1),
c(0,0,1,1,0,0,1,1,0,0),
c(0,1,0,0,0,1,1,1,0,1),
c(1,0,0,0,1,1,1,1,0,1),
c(1,1,0,1,0,1,1,0,0,1),
c(1,1,0,1,0,1,1,1,1,1)
),
nrow=12,
ncol=10
)

#1000 X 100 binary data (This runs like 10 sec on my mac)
large_binary<-c()
for(i in 1:1000){
  large_binary<-rbind(large_binary,sample(x=c(1,0),size=100,replace=TRUE))
}


##Bottom up experiment test data frame and hclustering obj
btUpTestDF=data.frame(id=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"),label=c("A","A","A","A","B","B","B","C","C","C"))
set.seed(101)
btUpTesScore=matrix(c(1+0.1*rnorm(10*4),2+0.1*rnorm(10*3),3+0.1*rnorm(10*3)),ncol=10,nrow=10,byrow=T)
rownames(btUpTesScore)=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
btUpTestHclust=hclust(dist(btUpTesScore))






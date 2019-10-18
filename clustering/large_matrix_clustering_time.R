#Goal: Test how long it takes when running large matrix clustering

mat_10<-matrix(rnorm(10^2,mean=0,sd=1),nrow=10,ncol=10)

mat_100<-matrix(rnorm(100^2,mean=0,sd=1),nrow=100,ncol=100)

mat_1000<-matrix(rnorm(1000^2,mean=0,sd=1),nrow=1000,ncol=1000)

mat_10000<-matrix(rnorm(10000^2,mean=0,sd=1),nrow=10000,ncol=10000)

mat_100000<-matrix(rnorm(100000^2,mean=0,sd=1),nrow=100000,ncol=100000)

#I don't think our strains and conditions can be more than this:
system.time(mat_1000000<-matrix(rnorm(1000000^2,mean=0,sd=1),nrow=1000000,ncol=1000000))

t1<-system.time(hclust_10<-hclust(dist(mat_10)))
t2<-system.time(hclust_100<-hclust(dist(mat_100)))
t3<-system.time(hclust_1000<-hclust(dist(mat_1000)))
t4<-system.time(hclust_10000<-hclust(dist(mat_10000)))
t5<-system.time(hclust_100000<-hclust(dist(mat_100000)))
t6<-system.time(hclust_1000000<-hclust(dist(mat_100000)))
                
dput(c(t1,t2,t3,t4,t5,t6),"time.txt")
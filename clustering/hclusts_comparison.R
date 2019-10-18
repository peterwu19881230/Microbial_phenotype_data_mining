#Compare hclusts (trees of clustering results)

#synthetic dataset -> scramble the rows


#tree generation



#Test codes
library("gtools")
hs<-c(1,1,1,1,100,100) #All hclusts in a vector
i<-length(hs) #change i to give more iterations
Identical<-c()
names<-c()
mat<-combinations(i,2,1:i) #combinations here returns a matrix

#Examine all (h[i],h[j] pairs)
for(j in 1:nrow(mat)){
    #comb_value[j]<-hs[mat[j,1]]*hs[mat[j,2]]
    Identical[j]<-identical(hs[mat[j,1]],hs[mat[j,2]])
    names[j]<-paste(c(mat[j,1],mat[j,2]),collapse="-")
}
names(Identical)<-names



#Grouping and showing stats
group<-list()
m<-1
for(i in 1:length(hs)){
    
    
    for(j in 1:length(hs)){
        group[m]<-c(Identical[Identical[TRUE]]) #put the same thing into the same group
        #if() to exit the for if there are 0 element to be compared
        Identical<-Identical[Identical[Identical[FALSE]]] #exclude the elements that were grouped
    m<-m+1
    }
}


#Demonstrate some ways to do basic plotting in R

#plot(): 1. give a scatter plot 2. Can be used to plot a dendrogram when input is a hclust or dendrogram obj
plot(1:10,rnorm(10))



#hist(): histogram
hist(rnorm(100))

#barplot(): barplot
bars=c(1,3,5)
barplot(bars)


#------Things that can be added to the above graphs------

y=rnorm(10)
plot(1:10,y)
lines(1:10,y) #lines(): segments connecting points
lines(1:10,y+0.05,col="grey") #col: specifies colors
lines(1:10,y-0.05,lwd="3") #lwd" specifies width
points(1:10,y+0.1,col="blue") #This adds points that have the same x coordinate but different y coordinate

##type=: p: points, l: lines, o: overplotted points and lines
plot(1:10,y,type="p")
plot(1:10,y,type="l")
plot(1:10,y,type="o")

#cex.lab: size of x,y text (default=1)
plot(1:10,y,cex.lab=2) #in a barplot cex.names should be used instead of cex.lab

#las: rotate the x-axis & y-axis
plot(1:10,y,las=2)


#lty: line type
plot(1:10,y,type="o",lty=2)
plot(1:10,y,type="o",lty=3)
plot(1:10,y,type="o",lty=4)
plot(1:10,y,type="o",lty=5)
plot(1:10,y,type="o",lty=6)

#legend(): add the legend
plot(1:10,y)
lines(1:10,y) #lines(): segments connecting points
lines(1:10,y+0.05,col="grey") 
lines(1:10,y-0.05,col="blue",lty=2) 
legend(x=1,col=c("grey","blue"),legend=c("grey line","blue line"),lty=c(1,2))



#ggplot using for loop

##data generation
df1=data.frame(x=rnorm(20),y=rnorm(20))
df2=data.frame(x=rnorm(20),y=rnorm(20))
df3=data.frame(x=rnorm(20),y=rnorm(20))




##old way 
p=ggplot()+geom_line(data=df1,aes(x,y),color="blue")+
          geom_line(data=df2,aes(x,y),color="red")+
          geom_line(data=df3,aes(x,y),color="black")


##new way
df_list=list(df1,df2,df3)
colors_=c("blue","red","black")
p_new=ggplot()
for(i in 1:length(df_list)){
  p_new=p_new+geom_line(data=df_list[[i]],aes(x,y),color=colors_[i])
}



##display the graphs
gridExtra::grid.arrange(p, p_new, ncol=2)



##Note: I don't know how to use for loop + display legend properly. Eg. the following doesn't show figure legend:
df_list=list(df1,df2,df3)
colors_=c("blue","red","black")
p_new=ggplot()
for(i in 1:length(df_list)){
  p_new=p_new+geom_line(data=df_list[[i]],aes(x,y,color=colors_[i]),color=colors_[i]) ##color argument inside aes() should draw legend but it doesn't
}

p_new






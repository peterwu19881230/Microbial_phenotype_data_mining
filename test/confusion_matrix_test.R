




df=melted_metrics_df



ggplot(data=df,aes(x=Order,y=value,group=variable),size=size,alpha=alpha)+geom_line()+theme_minimal()+theme(text = element_text(size=18))
  


melted_Pwy_metrics=melt(cbind(id=1:7914231,Pwy_metrics),id.vars="id")



#Function to get the graph. 
plot_metrics=function(samples,alpha,size,metric_df,x_lab="high |PCC| -- ranked pairs -- low |PCC|",y=""){
  
  
  melted_metrics_df=melt(cbind(Order=samples,metric_df[samples,]),id.vars="Order")
  
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot(data=df,aes(x=Order,y=value,group=variable),size=size,alpha=alpha)+geom_line()+theme_minimal()+theme(text = element_text(size=18))+
    labs(x = x_lab,y=y,aesthetic='custom text',color="Metrics")+
    scale_x_continuous(expand = c(0,0),labels=comma)
  
  #geom_line(data=cbind(no=samples,melted_rand_df[samples,]),aes(x=Order,y=value,group=variable),size=size,alpha=alpha)+
  #geom_line(data=data.frame(no=samples,specificity=metric_df$specificity[samples]),aes(no,specificity,color="specificity"),size=size,alpha=alpha)+
  #geom_line(data=data.frame(no=samples,precision=metric_df$precision[samples]),aes(no,precision,color="precision"),size=size,alpha=alpha)+
  #geom_line(data=data.frame(no=samples,accuracy=metric_df$accuracy[samples]),aes(no,accuracy,color="accuracy"),size=size,alpha=alpha)+
  
  #geom_line(data=data.frame(no=samples,sensitivity=rand_df$sensitivity[samples]),aes(no,sensitivity,color="sensitivity"),size=size,alpha=alpha)+
  #geom_line(data=data.frame(no=samples,specificity=rand_df$specificity[samples]),aes(no,specificity,color="specificity"),size=size,alpha=alpha)+
  #geom_line(data=data.frame(no=samples,precision=rand_df$precision[samples]),aes(no,precision,color="precision"),size=size,alpha=alpha)+
  #geom_line(data=data.frame(no=samples,accuracy=rand_df$accuracy[samples]),aes(no,accuracy,color="accuracy"),size=size,alpha=alpha)+
  
  
}
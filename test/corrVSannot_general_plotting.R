#Use ggplot to do correlation VS annotation


#The following is a mess. Have to redo it 


#(!)Have to solve this problem later: ranking for the above distances (pcc, spearman, euclidean) should be recoded because there are many same distances
plot_corrVSannot=function(cumSum,samples=1:4000){
  
  slope_negativeControl=tail(cumSum,n=1)/length(cumSum)
  
  
  p=ggplot()+theme_minimal()+theme(text = element_text(size=18))+ 
    geom_line(data=data.frame(no=samples,PCC_based_distance=cumSum$pcc[samples]),aes(no,PCC_based_distance,color="PCC based distance 1"),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Spearman_based_distance=cumSum$spearman[samples]),aes(no,Spearman_based_distance,color="Spearman based distance"),size=size,alpha=alpha)+
    geom_line(data=data.frame(no=samples,Euclidean_based_distance=cumSum$euclidean[samples]),aes(no,Euclidean_based_distance,color="Euclidean based distance"),size=size,alpha=alpha)+
    
    ##I tried to show legend for geom_abline but even the following failed
    ##https://groups.google.com/forum/#!topic/ggplot2/HtMK2ynltyw
    geom_abline(intercept=0,slope=slope_negativeControl,color="#D76E6E",linetype='dashed',size=1)+ #negative control
    geom_abline(intercept=0, slope=1,color="#1F89C2",linetype='dashed',size=1)+ #positive control
    scale_y_continuous(expand = c(0,0),labels=comma)+scale_x_continuous(expand = c(0,0),labels=comma)+
    labs(x = "Ranking of Distance",y="Cumulative No. of co-annotated pairs",aesthetic='custom text')+
    scale_colour_manual("Distance", 
                        breaks = c("PCC based distance", "Spearman based distance", "Euclidean based distance"),
                        values = c("PCC based distance"="#1aa548", "Spearman based distance"="#bc2ec9", "Euclidean based distance"="#baa118"))
  
  #+guides(colour=guide_legend(override.aes=list(shape=c(NA,NA,NA,NA,16,NA)))) #ref: https://stackoverflow.com/questions/15059832/problems-with-ggplot-scale-color-and-linetype-parameters
  
  
  
  return(p)
}
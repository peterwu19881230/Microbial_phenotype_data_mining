#GO boxplot

get_boxplot=function(similarity){
  
  df=strain1strain2_allAnnotations_allDistances
  mean(is.na(df$Wang_BP)) # ~56% of the similarity are NA
  
  ##the first box
  #------------------------------------------------
  sorted_WangBP=df$Wang_BP[order(df[[similarity]])]
  sorted_WangBP_naRemoved=sorted_WangBP[!is.na(sorted_WangBP)]
  #------------------------------------------------
  
  ##the second box
  #------------------------------------------------
  cutoffIndex=1946 #optimal cutoff determined by using all 5 annots (|PCC|=0.67)
  sorted_WangBP_subset=sorted_WangBP[1:cutoffIndex]
  sorted_WangBP_subset_naRemoved=sorted_WangBP_subset[!is.na(sorted_WangBP_subset)]
  #------------------------------------------------
  
  
  df_all=data.frame(All=sorted_WangBP_naRemoved)
  df1=data.frame(aboveCutoff=sorted_WangBP_subset_naRemoved)
  
  
  xlabs=c("All","|PCC|>0.67 ")
  
    ggplot() +
    geom_violin(data = df_all,aes(xlabs[1],All)) +
    geom_boxplot(data = df_all,aes(xlabs[1],All),width=0.1,outlier.shape = NA)+ #outlier.shape decides the shape of outliers. Here I don't let them show
    geom_violin(data = df1,aes(xlabs[2],aboveCutoff)) +
    geom_boxplot(data = df1,aes(xlabs[2],aboveCutoff),width=0.1,outlier.shape = NA)+
    scale_x_discrete("",limits=xlabs)+ 
    #I want the x axis to be empty. And if I don't use this, the order is not right
    scale_y_continuous("Semantic similarity")+
    theme(text=element_text(size=15), #if this is 25 some labels will overlap each other
          axis.text.y=element_text(size=25),
          axis.title=element_text(size=25))
}


#Using pcc
p1=get_boxplot("pcc")

#Using ternary MI
p2=get_boxplot("mi_ternary")

#Using ternary MI
p3=get_boxplot("mi_ternary_collapsedCond")

dir_of_workingScript=dirname(rstudioapi::getSourceEditorContext()$path)
ggsave(file=paste(dir_of_workingScript,"GOsim_boxplot_pcc.pdf",sep="/"),plot=p1,device="pdf",width=3*4,height=2*4)
ggsave(file=paste(dir_of_workingScript,"GOsim_boxplot_mi_ternary.pdf",sep="/"),plot=p2,device="pdf",width=3*4,height=2*4)
ggsave(file=paste(dir_of_workingScript,"GOsim_boxplot_mi_ternary_collapsedCond.pdf",sep="/"),plot=p3,device="pdf",width=3*4,height=2*4)

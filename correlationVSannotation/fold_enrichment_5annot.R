#Function to get the graph. 
#========================================================================================================================
plot_foldEnrichment=function(samples,alpha,size,four_annot_list,x_lab="high |PCC| -- ranked pairs -- low |PCC|",y,legend_title="Annotation set(s)"){
  
  ##Cols=c("#56B4E9","#F3518A","#BE70EA","#09D38A") 
  
  ##It's so annoying that ggplot() keeps reordering my list for the legend and the manually assigned colors follows after reordering
  Cols=c("#09D38A","#56B4E9","#BE70EA","#F3518A")  
  
  
  library(scales) #This is for label=comma (not using scientific notation and put commas) to work
  ggplot()+theme_minimal()+theme(legend.title=element_text(size=15),
                                 legend.text=element_text(size=15),
                                 axis.text.x = element_text(size=15),
                                 axis.text.y=element_text(size=15),
                                 axis.title=element_text(size=20))+ 
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[1]][samples]),aes(no,Metric,color=names(four_annot_list[1])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[2]][samples]),aes(no,Metric,color=names(four_annot_list[2])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[3]][samples]),aes(no,Metric,color=names(four_annot_list[3])),size=size)+
    geom_line(data=data.frame(no=samples,Metric=four_annot_list[[4]][samples]),aes(no,Metric,color=names(four_annot_list[4])),size=size)+
    
    
    labs(x = x_lab,y=y,aesthetic='custom text',color=legend_title)+
    scale_x_continuous(expand = c(0,0),labels=comma)+
    scale_color_manual(values=rep(Cols,2)) 
}
#========================================================================================================================


#pwy
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcc")]
Pwy=cumsum(df$Pwy[order(df$pcc)])
Pwy_foldEnrichment=Pwy/((1:7914231)*mean(df$Pwy))



#pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("pcomplex","pcc")]
pcomplex=cumsum(df$pcomplex[order(df$pcc)])
pcomplex_foldEnrichment=pcomplex/((1:7914231)*mean(df$pcomplex))



#pwy+pcomplex
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","pcc")]
TF=( (df$Pwy+df$pcomplex)==2 )
sum(TF) #230
PwyANDpcomplex=cumsum( TF[order(df$pcc)])
PwyANDpcomplex_foldEnrichment=PwyANDpcomplex/((1:7914231)*mean(TF))



#pwy + pwcomplex + KEGG modules + regulon + operon
df=strain1strain2_allAnnotations_allDistances[,c("Pwy","pcomplex","operon","regulator","kegg_modules","pcc")]
TF=( rowSums(df[,-6])==5 )
sum(TF) #106
All_annotSet=cumsum( TF[order(df$pcc)])
All_annotSet_foldEnrichment=All_annotSet/((1:7914231)*mean(TF))




four_annot_list=list("Same Pathways"=Pwy_foldEnrichment,
                     "Same protein complexes"=pcomplex_foldEnrichment,
                     "Same pathways & protein complexes"=PwyANDpcomplex_foldEnrichment,
                     "Same in all 5 sets"=All_annotSet_foldEnrichment)


##Set no. of samples to be plotted
##samples=1:200
samples=1:4000
##samples=seq(from=0,to=7914231,by=1989); samples[1]=1

plot_foldEnrichment(samples=samples,alpha=0.5,size=1,four_annot_list,x_lab="high |PCC| -- ranked pairs -- low |PCC|",y="Fold Enrichment")

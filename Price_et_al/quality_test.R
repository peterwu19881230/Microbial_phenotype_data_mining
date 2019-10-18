#Evaluate E. coli keio collection phenotypes done by this paper: https://www.nature.com/articles/s41586-018-0124-0

#Place to download the data: http://genomics.lbl.gov/supplemental/bigfit/
new_dat=read.table("Data/fit_logratios_good.tab.txt",quote="",fill=T,sep="\t",header = T)
new_dat_scoresOnly=new_dat[,-(1:4)]


#gcvP, H, T
cor(t(All_Data[c(1628,1720,1809),])) #Nichols' gcv
cor(t(new_dat_scoresOnly[2350:2352,])) #new data's gcv

#(?) Do 2 papers' conditions overlap?
#gcvP, H, T
#Nichols' id: gcvP:1628,gcvT: 1809, gcvH: 1720 
#new data id:(b2903: gcvP  b2904: gcvH b2905: gcvT)
combined_gcv=cbind(All_Data[c(1628,1720,1809),],new_dat_scoresOnly[2350:2352,])

#3061(b3774), 3107(b3959) is of the highest pcc in Nichols'(>0.95)
#=>b3774: 2999 b3959: 3125 in the new data
cor(t(new_dat_scoresOnly[c(2999,3125),])) #0.9327354 #Not bad

#2700(b2508) 3222(b3643) pcc=-0.64 (lowest in Nichols')
#=>b2508 is not in the new_data

#3643(b0128) 3644(b3961) pcc=-0.6251828
#=>b0128: 88 b3961: 3127
cor(t(new_dat_scoresOnly[c(88,3127),])) #-0.19167 #Not terrible but not very good

#398(b0812)    2571(b3476)  1.458438e-08
#=>b0812: 634  b3476: 2768
cor(t(new_dat_scoresOnly[c(634,2768),])) #0.02904758 #Not bad







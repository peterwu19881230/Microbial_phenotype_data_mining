#purD: 2170,purL: 2618 
purD=as.numeric(All_Data[2170,]) #(!Carefule) The last condition is NA
purL=as.numeric(All_Data[2618,])
plot(purD,purL,col = ifelse((purD >3.463|purD <(-3.463)) | (purL >3.463|purL <(-3.463))
                            ,'blue','black'),main="Fitness scores of purD V.S. purL's")
fit=lm(purL~purD)
abline(fit,col="blue")
abline(h=3.463,col="red")
abline(v=3.463,col="red")
abline(h=-3.463,col="red")
abline(v=-3.463,col="red")
cor(purD,purL,use="pairwise.complete.obs",method="pearson")

res=resid(fit)
abs_res=abs(res)
sorted_abs_res=sort(abs_res)

#Get the indices of conditions that have the smallest residuals (get the first 20)
sorted_abs_res[1:20]

#231         220         299         225         318         198 
#0.002631403 0.007657765 0.008553569 0.015648472 0.019410699 0.021028914 
#285         164          79         254         166          83 
#0.021602699 0.024519728 0.025062569 0.029307387 0.032374655 0.033072707 
#72         210         120         209         200         319 
#0.036221645 0.036235325 0.037937707 0.044354843 0.071369478 0.075951188 
#308         126 
#0.076904756 0.086321600 

#Find the corresponding conditions
colnames(All_Data)[c(231,220, 299,225,318,198,
                   285,164,79 ,254,166,83,
                   72,210,120,209,200,319,
                   308,126)]

#[1] "VERAPAMIL.1.0...UNSPECIFIED"      "TRITONX.0.01....UNSPECIFIED"     
#[3] "POLYMYXINB.2.0...UNSPECIFIED"     "TUNICAMYCIN.7.5...UNSPECIFIED"   
#[5] "TRIMETHOPRIM.0.3...UNSPECIFIED"   "SDS.3.0....UNSPECIFIED"          
#[7] "NOVOBIOCIN.6...UNSPECIFIED"       "PROPIDIUMIODIDE.1...UNSPECIFIED" 
#[9] "CHLOROPROMAZINE.6...UNSPECIFIED"  "CYCLOSERINED.16...UNSPECIFIED"   
#[11] "PROPIDIUMIODIDE.50...UNSPECIFIED" "CHOLATE.2.0....UNSPECIFIED"      
#[13] "CERULENIN.1.0...UNSPECIFIED"      "THEOPHYLLINE.10...UNSPECIFIED"   
#[15] "FOSFOMYCIN.1.0...UNSPECIFIED"     "TAUROCHOLATE.1.0....UNSPECIFIED" 
#[17] "STREPTOMYCIN.0.05...UNSPECIFIED"  "TRIMETHOPRIM.0.4...UNSPECIFIED"  
#[19] "SPIRAMYCIN.5...UNSPECIFIED"       "GLUFOSFOMYCIN.0.05...UNSPECIFIED"




#Do this in the spyder console first: !pip install pyreadr
import pyreadr
import os
#exec(open("/Users/peterwu/Google Drive File Stream/Shared Drives/HuLab/Nichols_Data_mining/machine_learning/NeuroNetwork/NeuroNetworkPipeline.py").read()) #import the function I defined
exec(open("G://Shared Drives/HuLab/Nichols_Data_mining/machine_learning/NeuroNetwork/NeuroNetworkPipeline.py").read()) #import the function I defined

#os.chdir("/Users/peterwu/Google Drive File Stream/Shared Drives/HuLab/Nichols_Data_mining/")
os.chdir("G://Shared Drives/HuLab/Nichols_Data_mining")


#Use label matrix
#===============================================================================


#Load the .RData (based on: https://stackoverflow.com/questions/21288133/loading-rdata-files-into-python)
objects = pyreadr.read_r('Data/sourced/All_Data_NAimputed.RData')
NewPheno=objects["All_Data_NAimputed"]


objects = pyreadr.read_r('Data/pwyLabel_df.RData')
Label=objects["pwyLabel_df"]




type(NewPheno) #pandas.core.frame.DataFrame
type(Label) #pandas.core.frame.DataFrame

NewPheno.shape 
Label.shape 





metrics_df=NeuroNetworkPipeline(NewPheno,Label)
metrics_df.to_csv("machine_learning/NeuroNetwork/pathway.csv") #doesn't output if run using tensorflow backend on my pc

#Do this in terminal: pip install pyreadr
import pyreadr
import os
exec(open("/Users/peterwu/Google Drive File Stream/Shared Drives/HuLab/Nichols_Data_mining/test/machine_learning/SupportVectorMachine/SupportVectorMachinePipeline.py").read()) #import the function I defined

os.chdir("/Users/peterwu/Google Drive File Stream/Shared Drives/HuLab/Nichols_Data_mining/")



#Use label matrix
#===============================================================================


#Load the .RData (based on: https://stackoverflow.com/questions/21288133/loading-rdata-files-into-python)
objects = pyreadr.read_r('Data/sourced/All_Data_NAimputed.RData')
NewPheno=objects["All_Data_NAimputed"]


objects = pyreadr.read_r('Data/operonLabel_df.RData')
Label=objects["operonLabel_df"]




type(NewPheno) #pandas.core.frame.DataFrame
type(Label) #pandas.core.frame.DataFrame

NewPheno.shape 
Label.shape 





metrics_df=SupportVectorMachinePipeline(NewPheno,Label)
metrics_df.to_csv("test/machine_learning/SupportVectorMachine/operon.csv")

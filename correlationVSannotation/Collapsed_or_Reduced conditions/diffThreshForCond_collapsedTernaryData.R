#Create necessary objs for further CorrVSAnnot analysis

##All genes using Nichols' fig.S1 pwys


#MHD3
dat=Ternary_Data_324cutff_condCollapsed
attribute_list=attr_list(ECK.Pathway$ids,ECK.Pathway$Data)

start.time = Sys.time()
pwy_324cutff_collapsedCond.MHD3=dist_TF_cumsum_MHD3(dat,attribute_list=attribute_list) #Have to give the argument name otherwise attribute_list is assigned to the second argument "ternary"
end.time = Sys.time()
end.time - start.time #Time difference of 27.32015 mins
#save(pwy_324cutff_collapsedCond.MHD3,file="Data/sourced/pwy_324cutff_collapsedCond.MHD3.RData")

#MHD4
dat=Ternary_Data_324cutff_condCollapsed
attribute_list=attr_list(ECK.Pathway$ids,ECK.Pathway$Data)

start.time = Sys.time()
pwy_324cutff_collapsedCond.MHD4=dist_TF_cumsum_MHD4(dat,attribute_list=attribute_list) #Have to give the argument name otherwise attribute_list is assigned to the second argument "ternary"
end.time = Sys.time()
end.time - start.time #Time difference of 34.84176 mins
#save(pwy_324cutff_collapsedCond.MHD4,file="Data/sourced/pwy_324cutff_collapsedCond.MHD4.RData")














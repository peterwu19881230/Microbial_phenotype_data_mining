attribute_list=attribute_list=attr_list(ECK.Pathway$ids,ECK.Pathway$Data)

start.time = Sys.time()
nonRedun.140cond.Ternary_Data_MHD3=dist_TF_cumsum_MHD3(nonRedun.140cond.Ternary_Data,attribute_list=attribute_list)
end.time = Sys.time()
end.time - start.time #Time difference of 15.68224 mins
#save(nonRedun.140cond.Ternary_Data_MHD3,file="Data/sourced/nonRedun.140cond.Ternary_Data_MHD3.RData")
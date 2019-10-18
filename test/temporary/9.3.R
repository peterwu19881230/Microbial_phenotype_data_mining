###(Updated 9/3/2019) manually correct the fake ECKs made by Nichols (search for gene names in EcoCyc and put back the right ECK)

###====================================================================================================================================================
###The list of 22 fake ECKs (I pulled them out by manually looking at the sorted ECK of the original data file that Nichols' provides) -> correct ECK
###ECK4466-MOKC ECK0018
###ECK4472-YOAI ECK1786
###ECK5000-SROH ECK4505
###ECK5001-SGRT ECK4477
###ECK5002-ISTR-1 G0-10202 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5003-RYBD ECK4621
###ECK5004-RYEF ECK4574
###ECK5005-TP2 G0-8894 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5006-TPKE70 G0-8906 This doesn't have ECK number. Use this as the EcoCyc ID
###ECK5007-YKGR ECK4486
###ECK5008-YMIB ECK4487
###ECK5009-YMJD ECK4488
###ECK5010-YNBG ECK4489
###ECK5011-YOAJ ECK4490
###ECK5012-YOAK ECK4491
###ECK5013-YOBI ECK4492
###ECK5014-YOEI ECK4493
###ECK5015-YOHP ECK4494
###ECK5016-YPDK ECK4495
###ECK5017-YQCG ECK4497
###ECK5018-YQEL ECK4498
###ECK5019-YQFG ECK4499


###correct the rows of ECK_1st_table for the above strains

##correct the 2nd, 3rd columns of "id_ECKs_CorrectedECKs_AssociatedGeneNames" generated at the end of clean_names.R (associated_gene_names shouldn't have to change)
fake_ECK_genes=c("ECK4466-MOKC","ECK4472-YOAI","ECK5000-SROH","ECK5001-SGRT","ECK5002-ISTR-1","ECK5003-RYBD","ECK5004-RYEF","ECK5005-TP2","ECK5006-TPKE70",
                "ECK5007-YKGR","ECK5008-YMIB","ECK5009-YMJD","ECK5010-YNBG","ECK5011-YOAJ","ECK5012-YOAK","ECK5013-YOBI","ECK5014-YOEI","ECK5015-YOHP",
                "ECK5016-YPDK","ECK5017-YQCG","ECK5018-YQEL","ECK5019-YQFG")

corrected_ECK=c("ECK0018","ECK1786","ECK4505","ECK4477","","ECK4621","ECK4574","","","ECK4486","ECK4487","ECK4488","ECK4489","ECK4490","ECK4491","ECK4492","ECK4493","ECK4494","ECK4495",
                "ECK4497","ECK4498","ECK4499")

temp=id_ECKs_CorrectedECKs_AssociatedGeneNames
for(i in 1:length(fake_ECK_genes)){
  index_=grep(fake_ECK_genes[[i]],id_ECKs_CorrectedECKs_AssociatedGeneNames$originalECKs)
  temp[index_,2]=str_replace(fake_ECK_genes[[i]],"^ECK[0-9]{4}",corrected_ECK[[i]])
  temp[index_,3]=temp[index_,2]
}

id_ECKs_CorrectedECKs_AssociatedGeneNames=temp
###====================================================================================================================================================








#All_Data obj contains the conditions in the same order as in allData.csv (allData.csv is directly derived from Nichols' supplementary data that contain scores)
condName=names(All_Data)

#Manual collapse of conditions by indices
##Note: the names of those indices are exactly the same as those in this Nichols' supplementary: TableS1_ListOfConditions.xls
uniqueChemIndex=list(
  'SDS+EDTA'=1:3,'Cold shock'=4:6,'Heat shock'=7:10,
  A22=11:14,'Actinomycin D'=15:18,'Acetate (M9)'=19,
  Acriflavine=20:21,Amikacin=22:24,Amoxicillin=25:28,
  Ampicillin=29:32,Anaerobic=33,Azidothymidine=34:36,
  Azithromycin=37:39,Benzalkonium=40:42,Bicyclomycin=43:44,
  'Bile salts'=45:48,Bleomycin=49:52,'CCCP (Carbonyl cyanide 3-chlorophenylhydrazone)'=53:55,
  'CHIR-090'=56:60,Cefoxitin=61:64,'Calcofluor (F3543    
Fluorescent Brightener 28) '=65, #I intended to leave this guy like this (yes, an Enter after F3543) to match the name of the original excel file
  'Cecropin B'=66:67,Cefsulodin=68:71,Cerulenin=72:75,Chlorpromazine=76:79,Cholate=80:83,
  Cisplatin=84:86,Clarythromycin=87:90,'Cobalt stress-CoCl2'=91:92,
  'Copper stress-CuCl2'=93:95,Deoxycholate=96:98,Dibucaine=99:101,
  Doxorubicin=102:103,EDTA=104:106,'Epigallocatechin gallate (EGCG)'=107:109,
  Epinephrine=110:112,Erythromycin=113:116,'Ethidium Bromide'=117:119,
  Fosfomycin=120,'Fusidic acid'=121:124,'N-acetyl Glucosamine'=125,
  'Fosfomycin +Glucose 6P'=126:127,'Glucosamine (M9)'=128,'Glucose (M9)'=129,
  'Glycerol (M9)'=130,Hydroxyurea=131:133,Indolicidin=134,Isoniazid=135:137,
  Levofloxacin=138,MMS=139,'Maltose (M9)'=140, Mecillinam=141:144,
  Methotrexate=145:146,Minocycline=147:149,'Mitomycin C'=150,'NH4Cl (MOPS)'=151,
  'Nickel stress-NiCl2'=152:153,Nigericin=154:156,Nitrofurnatoin=157:161, #Their typo: Nitrofurnatoin -> Nitrofurantoin
  Norepinephrine=162:163,'Propidium iodide'=164:166,'Phenazine methosulfate (PMS)'=167:169,
  'Paraquat dichloride'=170:174,'Hydrogen peroxide'=175:178,Phleomycin=179:181,Procaine=182:185,
  Puromycin=186:188,Pyocyanin=189:191,Radicicol=192:194,
  SDS=195:199,Streptomycin=200,Streptonigrin=201:203,
  'Succinate (M9)'=204,Sulfamonomethoxine=205:206,Taurocholate=207:209,
  Theophylline=210:211,Thiolactomycin=212:214,Tobramycin=215:218,
  'Triclosan/Irgasan'=219,'Triton X-100'=220:222,Tunicamycin=223:225,
  Vancomycin=226:228,Verapamil=229:231,Aztreonam=232:233,
  Bacitracin=234:236,Carbenicillin=237:239,'Cefsulodin + Mecillinam'=240,
  Cefaclor=241:243,Ceftazidime=244:246,Chloramphenicol=247:250,
  Ciprofloxacin=251:253,'Cycloserine D'=254,Doxycycline =255:258,
  EGTA=259:262,EtOH=263:265,Gentamicin=266:267,
  'Iron excess-FeSO4'=268,'Iron starvation-FeSO4'=269,NaCl=270:273,
  'Nalidixic acid'=274:277,Norfloxacin=278:280,Novobiocin=281:286,
  Oxacillin=287:289,'basic pH (TAPS)'=c(290,295:297),'acidic pH (MES-HOMOPIPES)'=291:294,
  'Polymyxin B'=298:301,Rifampicin=302:303,Spectinomycin=304:305,
  Spiramycin=306:308,Sulfamethizole=309:311,Tetracycline=312:315,
  Trimethoprim=316:319,'Trimethoprim + Sulfamethizole'=320,UV=321:324
  )

save(uniqueChemIndex,file="Data/sourced/uniqueChemIndex.RData")


##Code to verify that there is no error:
length(uniqueChemIndex) #Should be 114 (as described in the papaer and listed in TableS1_ListOfConditions.xls)
sapply(uniqueChemIndex,length) %>% sum #Should be 324

##The distribution of No. of conditions for all treatments
table(sapply(uniqueChemIndex,length) %>% as.character)

### Histogram
NoOfCond=data.frame(No_Of_Conditions=sapply(uniqueChemIndex,length))
ggplot(NoOfCond,aes(No_Of_Conditions))+
  geom_histogram(bins=6)+
  theme_minimal()
  










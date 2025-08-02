# WGCNA-PhosphProteomeData
WGCNA is used to analyze phospho-proteome data instead of gene expression data, the aim is to looking for phosphorylation that co-occurred across a variety of samples. 
Data preprocessing: Mass spec intensity values -> normalization by substration method -> log transformation -> filter rows based on 70% valid value -> impute missing values based on normalization distribution
The preprocessed data then go into WGCNA workflow as show in the script file

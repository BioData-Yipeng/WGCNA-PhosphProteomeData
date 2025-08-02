# WGCNA-PhosphProteomeData
WGCNA is used to analyze phospho-proteome data instead of gene expression data, the aim is to looking for phosphorylation that co-occurred across a variety of samples. 


This is the dataset from our phospho-proteomic experiments of the Rett syndromne animal efficacy study, which contains 4 groups, wt, het(disease), het+low dose compound, het+high dose compound. Samples were mouse whole brain tissues snap freezed in liquid nitrogen.

I found these dataset can also be used to do a WGCNA analysis instead of the classic proteomic workflows. 

The WGCNA analysis is useful to identify phospho site what are co-regulated(in the same module), and relationship between each module or genes to clinical traits such as isoprostan levels we have data from the same individual animals.

Data preprocessing: Mass spec intensity values -> normalization by substration method -> log transformation -> filter rows based on 70% valid value -> impute missing values based on normalization distribution

The preprocessed data then go into WGCNA workflow as show in the script file

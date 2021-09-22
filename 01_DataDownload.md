---
title: "Data Download"
output: html_notebook
---


```r
# libraries
library(tidyverse)
library(TCGAbiolinks)
library(kableExtra)
```

```r
# retrieve a list of TCGA projects
getGDCprojects() %>% 
    dplyr::select("project_id", "tumor", "name") %>%
    filter(str_detect(project_id, "TCGA")) %>% 
    kbl() %>%
    kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> project_id </th>
   <th style="text-align:left;"> tumor </th>
   <th style="text-align:left;"> name </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> TCGA-UCEC </td>
   <td style="text-align:left;"> UCEC </td>
   <td style="text-align:left;"> Uterine Corpus Endometrial Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-LGG </td>
   <td style="text-align:left;"> LGG </td>
   <td style="text-align:left;"> Brain Lower Grade Glioma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-SARC </td>
   <td style="text-align:left;"> SARC </td>
   <td style="text-align:left;"> Sarcoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-PAAD </td>
   <td style="text-align:left;"> PAAD </td>
   <td style="text-align:left;"> Pancreatic Adenocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-ESCA </td>
   <td style="text-align:left;"> ESCA </td>
   <td style="text-align:left;"> Esophageal Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-PRAD </td>
   <td style="text-align:left;"> PRAD </td>
   <td style="text-align:left;"> Prostate Adenocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-LAML </td>
   <td style="text-align:left;"> LAML </td>
   <td style="text-align:left;"> Acute Myeloid Leukemia </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-KIRC </td>
   <td style="text-align:left;"> KIRC </td>
   <td style="text-align:left;"> Kidney Renal Clear Cell Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-PCPG </td>
   <td style="text-align:left;"> PCPG </td>
   <td style="text-align:left;"> Pheochromocytoma and Paraganglioma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-HNSC </td>
   <td style="text-align:left;"> HNSC </td>
   <td style="text-align:left;"> Head and Neck Squamous Cell Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-OV </td>
   <td style="text-align:left;"> OV </td>
   <td style="text-align:left;"> Ovarian Serous Cystadenocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-GBM </td>
   <td style="text-align:left;"> GBM </td>
   <td style="text-align:left;"> Glioblastoma Multiforme </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-UCS </td>
   <td style="text-align:left;"> UCS </td>
   <td style="text-align:left;"> Uterine Carcinosarcoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-MESO </td>
   <td style="text-align:left;"> MESO </td>
   <td style="text-align:left;"> Mesothelioma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-TGCT </td>
   <td style="text-align:left;"> TGCT </td>
   <td style="text-align:left;"> Testicular Germ Cell Tumors </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-KICH </td>
   <td style="text-align:left;"> KICH </td>
   <td style="text-align:left;"> Kidney Chromophobe </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-READ </td>
   <td style="text-align:left;"> READ </td>
   <td style="text-align:left;"> Rectum Adenocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-UVM </td>
   <td style="text-align:left;"> UVM </td>
   <td style="text-align:left;"> Uveal Melanoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-THCA </td>
   <td style="text-align:left;"> THCA </td>
   <td style="text-align:left;"> Thyroid Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-LIHC </td>
   <td style="text-align:left;"> LIHC </td>
   <td style="text-align:left;"> Liver Hepatocellular Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-THYM </td>
   <td style="text-align:left;"> THYM </td>
   <td style="text-align:left;"> Thymoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-CHOL </td>
   <td style="text-align:left;"> CHOL </td>
   <td style="text-align:left;"> Cholangiocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-DLBC </td>
   <td style="text-align:left;"> DLBC </td>
   <td style="text-align:left;"> Lymphoid Neoplasm Diffuse Large B-cell Lymphoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-KIRP </td>
   <td style="text-align:left;"> KIRP </td>
   <td style="text-align:left;"> Kidney Renal Papillary Cell Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-BLCA </td>
   <td style="text-align:left;"> BLCA </td>
   <td style="text-align:left;"> Bladder Urothelial Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-BRCA </td>
   <td style="text-align:left;"> BRCA </td>
   <td style="text-align:left;"> Breast Invasive Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-COAD </td>
   <td style="text-align:left;"> COAD </td>
   <td style="text-align:left;"> Colon Adenocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-CESC </td>
   <td style="text-align:left;"> CESC </td>
   <td style="text-align:left;"> Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-LUSC </td>
   <td style="text-align:left;"> LUSC </td>
   <td style="text-align:left;"> Lung Squamous Cell Carcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-STAD </td>
   <td style="text-align:left;"> STAD </td>
   <td style="text-align:left;"> Stomach Adenocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-SKCM </td>
   <td style="text-align:left;"> SKCM </td>
   <td style="text-align:left;"> Skin Cutaneous Melanoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-LUAD </td>
   <td style="text-align:left;"> LUAD </td>
   <td style="text-align:left;"> Lung Adenocarcinoma </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TCGA-ACC </td>
   <td style="text-align:left;"> ACC </td>
   <td style="text-align:left;"> Adrenocortical Carcinoma </td>
  </tr>
</tbody>
</table>
The focus of this analysis is on two subtypes of lung cancer, namely:

1) Lung Squamous Cell Carcinoma which is identified in TCGA as **LUSC**.
2) Lung Adenocarcinoma which is identified in TCGA as **LUAD**.

Let's start with LUSC and download the gene expression data (raw count)




```r
TCGAbiolinks:::getProjectSummary("TCGA-LUSC")$case_count
```

```
## [1] 504
```

TCGA LUSC cohort consists of **504** cases and include the following data categories:

- Sequencing Reads
- Transcriptome
- Simple Nucleotide Variation (SNV)
- Copy Number Variation (CNV)
- DNA Methylation
- Clinical and
- Biospecimen


```r
TCGAbiolinks:::getProjectSummary("TCGA-LUSC")$data_categories %>%
    kbl() %>%
    kable_classic(full_width = FALSE, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> file_count </th>
   <th style="text-align:right;"> case_count </th>
   <th style="text-align:left;"> data_category </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 2248 </td>
   <td style="text-align:right;"> 504 </td>
   <td style="text-align:left;"> Sequencing Reads </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 577 </td>
   <td style="text-align:right;"> 504 </td>
   <td style="text-align:left;"> Clinical </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3147 </td>
   <td style="text-align:right;"> 504 </td>
   <td style="text-align:left;"> Copy Number Variation </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2630 </td>
   <td style="text-align:right;"> 504 </td>
   <td style="text-align:left;"> Biospecimen </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 573 </td>
   <td style="text-align:right;"> 503 </td>
   <td style="text-align:left;"> DNA Methylation </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2699 </td>
   <td style="text-align:right;"> 504 </td>
   <td style="text-align:left;"> Transcriptome Profiling </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4494 </td>
   <td style="text-align:right;"> 497 </td>
   <td style="text-align:left;"> Simple Nucleotide Variation </td>
  </tr>
</tbody>
</table>




```r
# view the results of the query
head(getResults(query_LUSC.exp))
```

```
##                                     id data_format access
## 1 a4428e1e-de4c-4011-abc1-9ecf1262a51d         TXT   open
## 2 29f8cd8d-6125-47a2-9bb3-d7fab9220fae         TXT   open
## 3 9cc0a963-d16f-46b0-a41c-55c3d1cdacf9         TXT   open
## 4 8fb3a2ca-39af-41f8-94a5-4954f92a83cc         TXT   open
## 5 78a55e36-22a6-4204-967c-0c62b61fd88b         TXT   open
## 6 5f4438e3-4e34-4528-a8e5-285364a12b81         TXT   open
##                          cases
## 1 TCGA-94-7943-01A-11R-2187-07
## 2 TCGA-68-8251-01A-11R-2296-07
## 3 TCGA-33-A5GW-01A-11R-A27Q-07
## 4 TCGA-85-8070-01A-11R-2247-07
## 5 TCGA-85-8479-01A-11R-2326-07
## 6 TCGA-22-5482-01A-01R-1635-07
##                                              file_name           data_category
## 1 3bdbc229-cbc8-401b-836c-90e76ff3866b.htseq.counts.gz Transcriptome Profiling
## 2 db2ba932-0f99-4e43-ab09-355bed4651e7.htseq.counts.gz Transcriptome Profiling
## 3 e30a275e-eaf3-43d4-b29c-dc75d77a0c9e.htseq.counts.gz Transcriptome Profiling
## 4 6219a10b-0637-4332-a4f0-118b5d434d74.htseq.counts.gz Transcriptome Profiling
## 5 0ad0b568-345f-4be5-a093-c6b08ed5e953.htseq.counts.gz Transcriptome Profiling
## 6 84f304ea-8b50-41e9-89c8-bfca6668b746.htseq.counts.gz Transcriptome Profiling
##                                 submitter_id            type
## 1 3bdbc229-cbc8-401b-836c-90e76ff3866b_count gene_expression
## 2 db2ba932-0f99-4e43-ab09-355bed4651e7_count gene_expression
## 3 e30a275e-eaf3-43d4-b29c-dc75d77a0c9e_count gene_expression
## 4 6219a10b-0637-4332-a4f0-118b5d434d74_count gene_expression
## 5 0ad0b568-345f-4be5-a093-c6b08ed5e953_count gene_expression
## 6 84f304ea-8b50-41e9-89c8-bfca6668b746_count gene_expression
##                   created_datetime file_size                           md5sum
## 1 2016-05-29T10:56:22.433827-05:00    256920 2dc1aa6b8731cc72689a3c992da1a511
## 2 2016-05-30T18:52:42.175781-05:00    254753 73a0bc107c45e693c8987693fd3b590a
## 3 2016-05-29T10:50:02.552154-05:00    251847 bca8e5a8bb775f31a19a2bcce093a712
## 4 2016-05-29T10:41:18.875484-05:00    251372 e49bcbee663e5bb33107b93465419e47
## 5 2016-05-30T18:26:52.287123-05:00    252953 491d89645cf14389369ba2b6e6e049dc
## 6 2016-05-29T11:08:45.032411-05:00    255061 6ceda4104bd4d689a2f3f951ca83dd23
##                   updated_datetime                              file_id
## 1 2018-09-06T16:31:22.449726-05:00 a4428e1e-de4c-4011-abc1-9ecf1262a51d
## 2 2018-09-06T16:31:22.449726-05:00 29f8cd8d-6125-47a2-9bb3-d7fab9220fae
## 3 2018-09-06T16:31:22.449726-05:00 9cc0a963-d16f-46b0-a41c-55c3d1cdacf9
## 4 2018-09-06T16:31:22.449726-05:00 8fb3a2ca-39af-41f8-94a5-4954f92a83cc
## 5 2018-09-06T16:31:22.449726-05:00 78a55e36-22a6-4204-967c-0c62b61fd88b
## 6 2018-09-06T16:31:22.449726-05:00 5f4438e3-4e34-4528-a8e5-285364a12b81
##                        data_type    state experimental_strategy version
## 1 Gene Expression Quantification released               RNA-Seq       1
## 2 Gene Expression Quantification released               RNA-Seq       1
## 3 Gene Expression Quantification released               RNA-Seq       1
## 4 Gene Expression Quantification released               RNA-Seq       1
## 5 Gene Expression Quantification released               RNA-Seq       1
## 6 Gene Expression Quantification released               RNA-Seq       1
##   data_release   project                          analysis_id analysis_state
## 1  12.0 - 29.0 TCGA-LUSC b5e89bf8-d38a-4053-aa25-5d751c825fc3       released
## 2  12.0 - 29.0 TCGA-LUSC 19825d40-ec5f-4d90-a231-173ac23187a5       released
## 3  12.0 - 29.0 TCGA-LUSC 8de1a986-1b2c-4dbb-85d9-294ec7b5a651       released
## 4  12.0 - 29.0 TCGA-LUSC b2eaff3c-04ea-4e78-b628-55f3fdf01a29       released
## 5  12.0 - 29.0 TCGA-LUSC 029f2232-936e-4be4-a9a8-c610f9307517       released
## 6  12.0 - 29.0 TCGA-LUSC 12f816ff-ae0a-44af-8a42-4a2a2906b777       released
##                        analysis_submitter_id
## 1 3bdbc229-cbc8-401b-836c-90e76ff3866b_count
## 2 db2ba932-0f99-4e43-ab09-355bed4651e7_count
## 3 e30a275e-eaf3-43d4-b29c-dc75d77a0c9e_count
## 4 6219a10b-0637-4332-a4f0-118b5d434d74_count
## 5 0ad0b568-345f-4be5-a093-c6b08ed5e953_count
## 6 84f304ea-8b50-41e9-89c8-bfca6668b746_count
##                 analysis_workflow_link analysis_workflow_type
## 1 https://github.com/NCI-GDC/htseq-cwl         HTSeq - Counts
## 2 https://github.com/NCI-GDC/htseq-cwl         HTSeq - Counts
## 3 https://github.com/NCI-GDC/htseq-cwl         HTSeq - Counts
## 4 https://github.com/NCI-GDC/htseq-cwl         HTSeq - Counts
## 5 https://github.com/NCI-GDC/htseq-cwl         HTSeq - Counts
## 6 https://github.com/NCI-GDC/htseq-cwl         HTSeq - Counts
##   analysis_workflow_version   sample_type is_ffpe cases.submitter_id
## 1                        v1 Primary Tumor   FALSE       TCGA-94-7943
## 2                        v1 Primary Tumor   FALSE       TCGA-68-8251
## 3                        v1 Primary Tumor   FALSE       TCGA-33-A5GW
## 4                        v1 Primary Tumor   FALSE       TCGA-85-8070
## 5                        v1 Primary Tumor   FALSE       TCGA-85-8479
## 6                        v1 Primary Tumor   FALSE       TCGA-22-5482
##   sample.submitter_id
## 1    TCGA-94-7943-01A
## 2    TCGA-68-8251-01A
## 3    TCGA-33-A5GW-01A
## 4    TCGA-85-8070-01A
## 5    TCGA-85-8479-01A
## 6    TCGA-22-5482-01A
```


```r
# list all columns in the results DF
colnames(getResults(query_LUSC.exp))
```

```
##  [1] "id"                        "data_format"              
##  [3] "access"                    "cases"                    
##  [5] "file_name"                 "data_category"            
##  [7] "submitter_id"              "type"                     
##  [9] "created_datetime"          "file_size"                
## [11] "md5sum"                    "updated_datetime"         
## [13] "file_id"                   "data_type"                
## [15] "state"                     "experimental_strategy"    
## [17] "version"                   "data_release"             
## [19] "project"                   "analysis_id"              
## [21] "analysis_state"            "analysis_submitter_id"    
## [23] "analysis_workflow_link"    "analysis_workflow_type"   
## [25] "analysis_workflow_version" "sample_type"              
## [27] "is_ffpe"                   "cases.submitter_id"       
## [29] "sample.submitter_id"
```


```r
# how many cases in total?
length(getResults(query_LUSC.exp)$cases)
```

```
## [1] 551
```


```r
# how many unique cases?
length(unique(getResults(query_LUSC.exp)$cases.submitter_id))
```

```
## [1] 501
```


```r
# let's check the sample types and the number of cases in each type
table(getResults(query_LUSC.exp)$sample_type)
```

```
## 
##       Primary Tumor Solid Tissue Normal 
##                 502                  49
```


```r
# let's download the data
GDCdownload(query_LUSC.exp,
            method = "api",
            directory = "../Data/GDCdata")
```

```
## Downloading data for project TCGA-LUSC
## Of the 551 files for download 551 already exist.
## All samples have been already downloaded
```

Let's look up another modality, taking into account the available resources on my machine





```r
# let's check the common patients between the gene expression and methylation data
LUSC.common <- intersect(substr(getResults(query_LUSC.exp, cols = "cases"), 1, 12),
                         substr(getResults(query_LUSC.met, cols = "cases"), 1, 12))
length(LUSC.common)
```

```
## [1] 372
```


```r
# let's repeat the queries using only the common patients
query_LUSC.exp <- GDCquery(project = 'TCGA-LUSC',
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type = "HTSeq - Counts",
                           barcode = LUSC.common)
```

```r
query_LUSC.met = GDCquery(project = "TCGA-LUSC",
                          data.category = "DNA Methylation",
                          platform = c("Illumina Human Methylation 450"),
                          barcode = LUSC.common)
```


```r
# download the gene expression data
GDCdownload(query_LUSC.exp,
            method = "api",
            directory = "../Data/GDCdata")
```

```
## Downloading data for project TCGA-LUSC
## Of the 419 files for download 419 already exist.
## All samples have been already downloaded
```

```r
# download the DNA methylation data
GDCdownload(query_LUSC.met,
            method = "api",
            directory = "../Data/GDCdata")
```

The methylation data size is in GBs, so probably a bad idea. Let's try something simpler.

First, let's query the gene expression data once more, without subsetting to `LUSC.common`


```r
# download the gene expression data
GDCdownload(query_LUSC.exp,
            method = "api",
            directory = "../Data/GDCdata")
```

```
## Downloading data for project TCGA-LUSC
## Of the 551 files for download 551 already exist.
## All samples have been already downloaded
```



```r
# prepare the data by importing the downloaded matrices (.htseq.counts) into RangedSummarisedExperiement
LUSC.exp <- GDCprepare(query_LUSC.exp,
                       directory = "../Data/GDCdata")
```



```r
# let's save the prepared object
saveRDS(object = LUSC.exp,
        file = "../Out/RData/LUSC_exp.RDS",
        compress = FALSE)
```



```r
sessionInfo()
```

```
## R version 4.1.1 (2021-08-10)
## Platform: x86_64-apple-darwin20.5.0 (64-bit)
## Running under: macOS Big Sur 11.6
## 
## Matrix products: default
## BLAS/LAPACK: /usr/local/Cellar/openblas/0.3.17/lib/libopenblasp-r0.3.17.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] kableExtra_1.3.4            TCGAbiolinks_2.20.0        
##  [3] viridis_0.6.1               viridisLite_0.4.0          
##  [5] DESeq2_1.32.0               SummarizedExperiment_1.22.0
##  [7] Biobase_2.52.0              MatrixGenerics_1.4.3       
##  [9] matrixStats_0.60.1          GenomicRanges_1.44.0       
## [11] GenomeInfoDb_1.28.4         IRanges_2.26.0             
## [13] S4Vectors_0.30.0            BiocGenerics_0.38.0        
## [15] forcats_0.5.1               stringr_1.4.0              
## [17] dplyr_1.0.7                 purrr_0.3.4                
## [19] readr_2.0.1                 tidyr_1.1.3                
## [21] tibble_3.1.4                ggplot2_3.3.5              
## [23] tidyverse_1.3.1            
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_2.0-2            ellipsis_0.3.2             
##   [3] visdat_0.5.3                rprojroot_2.0.2            
##   [5] XVector_0.32.0              fs_1.5.0                   
##   [7] rstudioapi_0.13             farver_2.1.0               
##   [9] bit64_4.0.5                 AnnotationDbi_1.54.1       
##  [11] fansi_0.5.0                 lubridate_1.7.10           
##  [13] xml2_1.3.2                  splines_4.1.1              
##  [15] R.methodsS3_1.8.1           cachem_1.0.6               
##  [17] geneplotter_1.70.0          knitr_1.34                 
##  [19] jsonlite_1.7.2              broom_0.7.9                
##  [21] annotate_1.70.0             dbplyr_2.1.1               
##  [23] png_0.1-7                   R.oo_1.24.0                
##  [25] compiler_4.1.1              httr_1.4.2                 
##  [27] backports_1.2.1             assertthat_0.2.1           
##  [29] Matrix_1.3-4                fastmap_1.1.0              
##  [31] cli_3.0.1                   htmltools_0.5.2            
##  [33] prettyunits_1.1.1           tools_4.1.1                
##  [35] gtable_0.3.0                glue_1.4.2                 
##  [37] GenomeInfoDbData_1.2.6      rappdirs_0.3.3             
##  [39] Rcpp_1.0.7                  cellranger_1.1.0           
##  [41] vctrs_0.3.8                 Biostrings_2.60.2          
##  [43] svglite_2.0.0               xfun_0.26                  
##  [45] rvest_1.0.1                 lifecycle_1.0.0            
##  [47] XML_3.99-0.7                zlibbioc_1.38.0            
##  [49] scales_1.1.1                hms_1.1.0                  
##  [51] RColorBrewer_1.1-2          curl_4.3.2                 
##  [53] memoise_2.0.0               gridExtra_2.3              
##  [55] downloader_0.4              biomaRt_2.48.3             
##  [57] stringi_1.7.4               RSQLite_2.2.8              
##  [59] highr_0.9                   genefilter_1.74.0          
##  [61] filelock_1.0.2              BiocParallel_1.26.2        
##  [63] systemfonts_1.0.2           rlang_0.4.11               
##  [65] pkgconfig_2.0.3             bitops_1.0-7               
##  [67] evaluate_0.14               TCGAbiolinksGUI.data_1.12.0
##  [69] lattice_0.20-44             labeling_0.4.2             
##  [71] bit_4.0.4                   tidyselect_1.1.1           
##  [73] here_1.0.1                  plyr_1.8.6                 
##  [75] magrittr_2.0.1              R6_2.5.1                   
##  [77] generics_0.1.0              DelayedArray_0.18.0        
##  [79] DBI_1.1.1                   pillar_1.6.2               
##  [81] haven_2.4.3                 withr_2.4.2                
##  [83] survival_3.2-13             KEGGREST_1.32.0            
##  [85] RCurl_1.98-1.4              modelr_0.1.8               
##  [87] crayon_1.4.1                utf8_1.2.2                 
##  [89] BiocFileCache_2.0.0         rmarkdown_2.11             
##  [91] tzdb_0.1.2                  progress_1.2.2             
##  [93] locfit_1.5-9.4              grid_4.1.1                 
##  [95] readxl_1.3.1                data.table_1.14.0          
##  [97] blob_1.2.2                  webshot_0.5.2              
##  [99] reprex_2.0.1                digest_0.6.27              
## [101] xtable_1.8-4                R.utils_2.10.1             
## [103] munsell_0.5.0
```


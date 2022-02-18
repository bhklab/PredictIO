# Leveraging Big Data of Immune Checkpoint Blockade Response Identifies Novel Potential Targets

------------------------------------------------
------------------------------------------------ 


*Yacine Bareche, Deirdre Kelly, Farnoosh Abbas-Aghababazadeh, Minoru Nakano, Parinaz Nasr Esfahani, Denis Tkachuk, Hassan Mohammad, Robert Samstein, Chung-Han Lee, Luc G. T. Morris, Philippe L. Bedard, Benjamin Haibe-Kains & John Stagg*

## Summary
<p align="justify">

Several genomics and gene expression signatures have been proposed as predictive biomarkers of clinical response to immune checkpoint blockade (ICB), with questionable reproducibility. We here report a large-scale comparative analysis of candidate biomarkers of ICB responses in a pan-cancer meta-analysis of over 3,500 ICB-treated patients representing 12 different tumor types. A web-application (predictIO.ca) was developed to allow researchers to further interrogate this data compendium. We tested the hypothesis that a de novo pan-cancer gene expression analysis would bring forth critical ICB resistance pathways and novel therapeutic targets. At the genomic level, we confirmed that non-synonymous tumour mutational burden (nsTMB) was significantly associated with ICB responses across tumor types, with the exception of kidney cancer. At the transcriptional level, 21 out of 39 published gene expression signatures were significantly associated with pan-cancer ICB responses. Strikingly, the predictive value of a de novo gene expression signature (referred to as PredictIO_100) composed of the top 100 genes most significantly associated with ICB responses at pan-cancer level was superior to nsTMB and other gene expression biomarkers. Within PredictIO_100, two genes, F2RL1 (encoding protease-activated receptor-2) and RBFOX2 (encoding RNA binding motif protein 9), were concomitantly associated with worse ICB clinical outcomes, T cell dysfunction in ICB-naive patients and resistance to dual PD-1/CTLA-4 blockade in preclinical mouse cancer models. Taken together, our study underlined the relative impact of candidate biomarkers of ICB responses in a large pan-cancer cohort, demonstrated the potential of de novo pan-cancer gene expression signatures and identified F2RL1, previously involved in tumor immune regulation, and RBFOX2, a critical regulator of epithelial-to-mesenchymal transition, as potential therapeutic targets to overcome ICB resistance. 
</p >

------------------------------------------------
------------------------------------------------ 

## System Requirements

### Hardware Requirements
Running the scripts present in this repository only requires a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

### Software Requirements

#### OS Requirements

Scripts should work on any updated operating systems with R version 4.1.0 or higher installed.  

#### Package dependencies

Users should install the following packages prior to running the run_analysis.R script, from an R terminal:

```
install.packages( c( 'BiocManager' , 'R.utils' , 'RColorBrewer' , 'SuppDists' , 'apcluster' , 'beeswarm' , 'calibrate' , 'coin' , 'corrplot' , 'data.table' , 'enrichR' , 'forestplot' , 'ggfortify' , 'ggplot2' , 'grid' , 'gridExtra' , 'meta' , 'metafor' , 'pROC' , 'pheatmap' , 'plotROC' , 'reshape2' , 'survcomp' , 'survival' , 'survminer' , 'tidyverse' , 'vcd' ) )

BiocManager::install( c( 'Biobase' , 'BiocVersion' , GSVA' , genefu' , survcomp' ) )
```

The installation of the dependencies should take approximately few minutes to completely run on a recommended computer.


## Running the Script

Data needs to first be downloaded from https://zenodo.org/record/6142546 and stored in the data directory.
Users can then run the script run_analysis.R present in the code directory :
```
cd code
Rscript run_analysis.R
```

The script should take approximately 2 hours to completely run on a recommended computer.
 

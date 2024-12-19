# WES_project

Scripts for statistics and plotting figures in "Whole Exome Sequencing Reveals an FCGBP Variant Associated with Spontaneous Intraabdominal Hemorrhage in Severe Acute Pancreatitis", published in iMetaOmics 2024.

- /fig2,3,5,6 # raw data and R scripts for each panel in figure
- /WES_SNV_calling # General pipeline for WES variant calling used in this study

wes.yaml file for reproducable environment:
```sh
conda env create -f wes.yaml
```

**Reference bundle** required by GATK can be obtained through:
```sh
##lftp should be installed previously
lftp ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
get mirror ***
```

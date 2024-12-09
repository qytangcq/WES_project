All analysis is based on softwares on Linux system or R.
We suggest to install all the packages and tools in a **conda** environment.


You can install **conda** through:
```sh
get https://mirrors.ustc.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

##customize conda configuration
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
conda config --set show_channel_urls yes

##create environment,-n to specify the name of environment
conda  create -n wes
```

softwares such as  sra-tools,fastqc,multiqc,vcftools,bcftools,qualimap,trim-galore,csvtk,plink(version 1.9 and 2),gatk etc. should be installed by `conda search *** ` and `conda install ***` or installed manually according to official guidance.

We provide .yaml file created from our own conda environment. You can create a environment and install the required packages simultaneously with:
```sh
conda env create -f wes.yaml
```




**Reference bundle** required by GATK can be obtained through:
```sh
##lftp should be installed previously
lftp ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
get mirror ***
```

**Here is an expample of how to download fastq data from SRA database**
 ```sh
 ##fastq-dump reuired
 fastq-dump -v --split-3 --gzip SRR645165

 ##trans sra file to fastq
 fastq-dump SRR645165 --split-3 --gzip -O ./ -v
 ```

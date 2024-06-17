# Phyrvm
![python3.8](https://img.shields.io/badge/python-3.8-brightgreen)



## Model



| ![Image](./info/pipeline.jpg)|
|:--:|
| Phyrvm automates RNA virus detection in three stages: RNA virus discovery, phylogenetic analysis, and phylogeny-based virus characterization. It outputs viral sequences, phylogenetic trees, and detailed reports, offering flexibility and accuracy in identifying putative viral sequences and their host associations. | 

## Installation
- python 3.8
- R 4.2

### Step 1: Install conda and third-party dependencies
Phyrvm requires third-party packages from the conda-forge and bioconda channels

```shell
conda install -c bioconda blast bbmap seqkit  mafft megahit trimal  pplacer phyml taxonkit  bowtie2 cd-hit
conda install taxonkit diamond==2.0.15  bowtie2 samtools==1.16.1
```
**Notes:**

- Version of the tool available for reference:
  - bbduk.sh：bbmap v39.01 ； Seqkit v2.4.0 ；bowtie2 v2.5.1 ；megahit v1.2.9
  - mafft v7.520 ；trimal v1.4.1 ；makeblastdb,blastn,blastp v2.13.0+
  - phyml v3.3.20211231 ；samtools v1.16.1 ；diamond v2.0.15

- The [taxonkit](https://bioinf.shenwei.me/taxonkit/download/) dataset should also be downloaded!


### Step 2: Install Phyrvm via pip

All python packages will be downloaded automatically!

```shell
pip install -i phyrvm
```

### Step 3: Install R and R package

```shell
#install R package
R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
packages=c("tidyverse"，"ggplot2","RColorBrewer","phangorn","networkD3","jsonlite","dplyr","networkD3","jsonlite")
ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))  
        install.packages(new.pkg)
    sapply(pkg, require, character.only = TRUE)
}
ipak(packages)
```
**Notes:**
*tidyverse* is based on *systemfonts* and you may need the following code to install it
```shell
conda install r-systemfonts
```

### Step 4: Downloading and configuring the database.

Phyrvm requires an environment variable named `PHYRVM_DB_PATH`, this is the parent directory for the following databases.

See below for specific database configurations.
```shell
#set PHYRVM_DB_PATH to environment variable
export PHYRVM_DB_PATH=/path/to/the/database/
```



**Notes:**

- The databases take up a lot of space, so make sure you have enough disk space. 
If you already have these databases, you can skip the download step and just configure them.

- The download speed of the database depends on the internet. You can also choose other download methods such as `ascp`.


#### rRNA

- 1.1 Download the file from this [link](https://zenodo.org/records/10435588/files/Phyrvm_rRNA_db.fasta.bz2?download=1&preview=1).

- 1.2 Unzip the file and Using ***[bowtie2](https://github.com/BenLangmead/bowtie2)*** to build the index.
    ```shell
    bunzip2 -cv Phyrvm_rRNA_db.fasta.bz2 > PHYRVM_DB_PATH/rRNA/Phyrvm_rRNA_db.fasta
    bowtie2-build PHYRVM_DB_PATH/rRNA/Phyrvm_rRNA_db.fasta PHYRVM_DB_PATH/rRNA/rRNA_cutout_ref
    ```

#### **PROT_ACC2TAXID**

  ```shell
  #Download the `PROT_ACC2TAXID` file
  wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
  wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5

  #Check for the file integrity
  md5sum -c prot.accession2taxid.gz.md5

  #Unzip the files and onfiguration
  gunzip -c prot.accession2taxid.gz > PHYRVM_DB_PATH/accession2taxid/prot.accession2taxid
  ```

#### [NCBI Non-Redundant Protein Database (NR)](./info/db_NR.md)



#### [NCBI Nucleotide Sequence Database (NT)](./info/db_NT.md)



## Usage
phyrvm *medthod* [options]

- medthod:
  - end_to_end
  - contigs_filter
  - phylogenetic_analysis

### Example

```shell
phyrvm end_to_end  -i 1.fastq -i2 2.fastq \
    -out_dir out_path  --threads 60 -classify_model All --keep-dup
	
phyrvm contigs_filter -i 1.fastq -i2 2.fastq \
    -out_dir out_path  --threads 60 
	
phyrvm phylogenetic_analysis -classify_i test/test_contig.fasta   \
	-out_dir out_path   -classify_model All  --threads 90 --keep-dup
```



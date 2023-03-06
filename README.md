# Welcome
This repository contains analysis scripts for Lee et al. Tumor Intrinsic Mechanisms of Antigen Escape to Anti-BCMA and Anti-GPRC5D Targeted Immunotherapies in Multiple Myeloma. Under review.

![](images/Myeloma_Immunotherapy_Antigen_Escape_Graphical_Abstract.png)

# Summary

## Abstract
B cell maturation antigen (BCMA) target loss is considered to be a rare event that mediates multiple myeloma (MM) resistance to anti-BCMA chimeric antigen receptor T cell (CAR T) or bispecific T cell engager (TCE) therapies. Emerging data report that downregulation of G protein coupled receptor family C group 5 member D (GPRC5D) protein often occurs at relapse after anti-GPRC5D CAR T. To examine the tumor intrinsic factors that promote MM antigen escape, we performed combined bulk and single cell whole genome sequencing/ copy number variation analysis of 25 patients treated with anti-BCMA and/ or -GPRC5D CAR T/ TCE. In two cases, MM relapse post TCE/CAR T was driven by BCMA negative clones harboring focal biallelic deletions at the TNFRSF17 locus at relapse or by selective expansion of pre-existing subclones with biallelic TNFRSF17 loss. In another three cases of MM relapse, we identified non-truncating missense mutations or in-frame deletions in the extracellular domain of BCMA which negates the efficacies of anti-BCMA TCEs, despite detectable surface BCMA protein expression. With respect to GPRC5D, we report the first four cases of MM relapse with biallelic mutations of GPRC5D following anti-GPRC5D TCE, including two cases with convergent evolution where multiple subclones lost GPRC5D through different somatic events. Our data support that immunoselection of BCMA or GPRC5D negative or mutantclones post CAR T/ TCE therapies may be more prevalent than currently perceived. The engineering and selection of immunotherapies in MM should account for target antigen structural and extracellular domain mutations.

# Data

## Single-cell RNA-Seq
NCBI GEO: [GSE226335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226335) <br/>
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE226335/suppl/GSE226335_RAW.tar
tar -xvf GSE226335_RAW.tar
```
NCBI SRA: SRPXXXXXX (provided once public) <br/>
```
source activate sratoolkit
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID
## Specify SRR_ID - obtained using SRA Run selector.
```

## Single-cell CNV-Seq
NCBI GEO: [GSE226327](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226327) <br/>
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE226327/suppl/GSE226327_RAW.tar
tar -xvf GSE226327_RAW.tar
```
NCBI SRA: SRPXXXXXX (provided once public) <br/>
```
source activate sratoolkit
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID
## Specify SRR_ID - obtained using SRA Run selector.
```

# Contact
Dr. Nizar Bahlis, MD (nbahlis@ucalgary.ca) <br/>
Arnie Charbonneau Cancer Institute, University of Calgary

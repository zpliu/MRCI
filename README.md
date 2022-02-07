# MRCI (Mixture model Reciprocal Causation Inference)
---
MRCI is a tool for estimating reciprocal causation between two phenotypes simultaneously using the genome-scale summary statistics of the two phenotypes and reference LD data.
## Installation Notes
### Required packages
Currently MRCI only supports Linux-based system. R needs to be installed (tested under R-3.5.1). 
```
library(roptim)
library(Rcpp)
library(argparse)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(plyr)
```
### LD data
Download the LD data file (~150MB) from [here](https://www.dropbox.com/s/hyi6huw14hg4r2w/LDwindow1MB_cutoff0.1.RData?dl=0) and put it into `LDdata` folder.
## Input format
The input of MRCI are GWAS summary statistic of the two phenotypes. The required column and header names are as follows:
```
CHR SNP BP NMISS BETA SE P A1 A2 MAF
```
- `CHR` & `BP`  chromosome and position (GRCh37).
- `SNP` rsid of SNPs.
- `NMISS`   sample sizes for each SNP.
- `BETA`    effec sizes from GWAS. Use lnOR values for a binary trait.
- `SE`  standard errors from GWAS. Use SE for lnOR for a binary trait. 
- `P`   P-values from GWAS.
- `A1`&`A2` alleles of SNPs. `A1` is the effect allele.
- `MAF` allele frequency for either `A1` or `A2`.

The two GWAS data could have sample overlap but must come from the same population (Current version only supports European-ancestry). 
Two example GWAS summary data are prepared in `example/`. 
## Run MRCI
Download the scripts and configure parameters in `pipeline.sh`
```
# +++++++++++++  Edit parameters from here  ++++++++++++++++++++++
WORKINGDIR="FOLDER-PATH-TO-SAVE-RESULTS"

threadnum=20 # thread number
Rscriptcmd=/PATH/TO/R_X.X.X/bin/Rscript   # In case of multiple R versions, assign the corresponding Rscript version

### Location of MRCI script folder
MRCI_scriptDIR=/PATH/TO/MRCI

### Location of GWAS summmary statistic data
gwas_Y1_file=$MRCI_scriptDIR/example/Y1.gwas.gz
gwas_Y2_file=$MRCI_scriptDIR/example/Y2.gwas.gz

### For a quantitative trait, set as '-100'; 
### For a binary trait, please change accordingly
prevalence_Y1=-100 # prevalence for trait Y1
caseProp_Y1=-100 # case proportion in total samples for trait Y1
prevalence_Y2=-100 # prevalence for trait Y2
caseProp_Y2=-100 # case proportion in total samples for trait Y2

### Give a human-readable name for each trait
Y1="Y1"
Y2="Y2"
# ++++++++++++++  Edit stop  +++++++++++++++++++
```
Note: Since MRCI uses whole-genome scale summary data, using multi-thread is always recommended in order to improve the calculation efficiency. Memory comsumption is often not demmanding unless the GWAS summary data is extremely large.

MRCI has four main steps for estimation:
- step 1: Harmonize two GWAS data
- step 2: Run EM estimation under each model
- step 3: Model averaging
- step 4: Output final estimation and plot major results

## Output files
MRCI will generate some output files, most of which are used as log files. 
The final results can be found in three files:
- `final.estimate.txt`  This contains the final estimates of all parameters. In this file, `delta12` is causal direction from phenotype `Y2 -> Y1`, `delta21` is the other causal direction.
- `pltnus.png`  This file shows the plots of nuisance parameter estimates of each model.
- `pltsum.png`  This file shows the plots of causal estimation as well as some summary results between full model and averaged model.
## More information
Authors: Zipeng Liu (The University of Hong Kong) and Yiming Qin (The University of Hong Kong)




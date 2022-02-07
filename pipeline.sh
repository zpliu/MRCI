#!/bin/bash

# +++++++++++++  Edit parameters from here  +++++++++++++++++++++++++++++++++++++
WORKINGDIR="FOLDER-PATH-TO-SAVE-RESULTS"

threadnum=20 # thread number
Rscriptcmd=/PATH/TO/R_X.X.X/bin/Rscript   # In case of multiple R versions, assign the corresponding Rscript version

### Location of MRCI script folder
MRCI_scriptDIR=/PATH/OF/MRCI_v1.0

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
# ++++++++++++++  Stop editing +++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [ ! -e $WORKINGDIR ]; then
	mkdir $WORKINGDIR
fi

### Prefix for output files
outprefix4harmonizedata=$WORKINGDIR/harmonize
outprefix4gen=$WORKINGDIR/gen
outprefix4coz=$WORKINGDIR/coz
outprefix4modavg=$WORKINGDIR/modavg
outprefix4pltnuisance=$WORKINGDIR/pltnus
outprefix4pltres=$WORKINGDIR/pltsum
outprefix4finalres=$WORKINGDIR/final


# - - - - - - - - - - - - - - - - - - - - - - - -
### Harmonize data
# - - - - - - - - - - - - - - - - - - - - - - - -
echo "Harmonize data"
$Rscriptcmd $MRCI_scriptDIR/harmonize.R \
			--scriptDIR $MRCI_scriptDIR \
			--gwas_Y1_file $gwas_Y1_file \
			--gwas_Y2_file $gwas_Y2_file \
			--prevalence_Y1 $prevalence_Y1 \
			--caseProp_Y1 $caseProp_Y1 \
			--prevalence_Y2 $prevalence_Y2 \
			--caseProp_Y2 $caseProp_Y2 \
			--name4Y1 $Y1 \
			--name4Y2 $Y2 \
			--out_prefix $outprefix4harmonizedata


# - - - - - - - - - - - - - - - - - - - - - - - -
### Run each model
# - - - - - - - - - - - - - - - - - - - - - - - -
model_list="
comp4full
comp3noh1
comp3noh2
comp3h1h2
comp2pleio
"

for model in $model_list
do
	echo ${model}

	echo "--> genetic estimation"
	$Rscriptcmd $MRCI_scriptDIR/main4gen.R \
				--harmonizedata $outprefix4harmonizedata.RData \
				--funDIR $MRCI_scriptDIR \
				--thread $threadnum \
				--model $model \
				--out_prefix $outprefix4gen

	echo "--> causal inference"
	$Rscriptcmd $MRCI_scriptDIR/main4coz.R \
				--harmonizedata $outprefix4harmonizedata.RData \
				--funDIR $MRCI_scriptDIR \
				--est_gen_file $outprefix4gen.$model.estimate.txt \
				--thread $threadnum \
				--model $model \
				--out_prefix $outprefix4coz
done


# - - - - - - - - - - - - - - - - - - - - - - - -
### Model averaging
# - - - - - - - - - - - - - - - - - - - - - - - -
$Rscriptcmd $MRCI_scriptDIR/main4modavg.R \
			--harmonizedata $outprefix4harmonizedata.RData \
			--funDIR $MRCI_scriptDIR \
			--cozest_prefix $outprefix4coz \
			--thread $threadnum \
			--out_prefix $outprefix4modavg


# - - - - - - - - - - - - - - - - - - - - - - - -
### Summarize final results
# - - - - - - - - - - - - - - - - - - - - - - - -
### Plot nuisance params
$Rscriptcmd $MRCI_scriptDIR/sumplt_nuisance.R \
			--name4Y1 $Y1 \
			--name4Y2 $Y2 \
			--outprefix4coz $outprefix4coz \
			--outprefix4modavg $outprefix4modavg \
			--outprefix4pltnuisance $outprefix4pltnuisance

### Plot major results
$Rscriptcmd $MRCI_scriptDIR/sumplt_res.R \
			--name4Y1 $Y1 \
			--name4Y2 $Y2 \
			--outprefix4harmonizedata $outprefix4harmonizedata \
			--outprefix4coz $outprefix4coz \
			--outprefix4modavg $outprefix4modavg \
			--outprefix4pltres $outprefix4pltres \
			--outprefix4finalres $outprefix4finalres


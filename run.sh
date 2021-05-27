#!/bin/bash

# <Replace with the description and/or purpose of this shell script.>
expFile=$1 # "data/DEG.AtGenExpress.signed_zstats.heat_shoots"
seedFile=$2 # "data/DEG.AtGenExpress.signed_binary.heat_shoots"
nwkFile=$3 # "data/Ath_TF_list.gene"
TFliFile=$4 # "data/Ath_TF_list.gene"
gSet=$5 # (optional) additional gene list file if the user wants
resD="result"

# EDIT: To run manually:
# expFile=single_data2/exp.2.txt
# seedFile=single_data2/bin.txt
# nwkFile=single_data2/Temp_Net.txt
# TFliFile=single_data2/TF_list.txt
# gSet=
# resD="result"

if [ -z ${expFile} ]; then
    expFile="data/DEG.AtGenExpress.signed_zstats.heat_shoots"
fi
if [ -z ${seedFile} ]; then
    seedFile="data/DEG.AtGenExpress.signed_binary.heat_shoots"
fi
if [ -z ${nwkFile} ]; then
	nwkFile="data/templateNetwork"
fi
if [ -z ${TFliFile} ]; then
	TFliFile="data/Ath_TF_list.gene"
fi

prefix="PropaNet"
n=$(head -n 1 $seedFile| awk '{print NF}')

if [ -e result/intermediate_results ]; then
    echo
else
    mkdir result/intermediate_results
fi

# Weighted template network construction
python network_weight.py \
    -nwk ${nwkFile} \
    -exp ${expFile} \
    -o ${resD}/intermediate_results/${prefix}.templateNetwork \
    || exit 1
# python network_weight.py -nwk data/templateNetwork -exp ${expFile} -o ${resD}/intermediate_results/${prefix}.templateNetwork || exit 1

# Extract optimal TF list
python TF_adding_NP_noCtrl.py \
    ${TFliFile} \
    ${resD}/intermediate_results/${prefix}.templateNetwork \
    ${expFile} \
    ${seedFile} \
    -p 10 \
    -cond ${prefix} \
    -outD ${resD}/intermediate_results \
    -c 5.9 \
    -coverNo 20000 \
    --serial \
    || exit 1

# Extract target genes of the optimal TFs
for ((i=1;i<=$[$n-1];i++)); do
    echo $i
    python Target_genes.py\
        ${resD}/intermediate_results/${prefix}.nwk.t$i\
        ${resD}/intermediate_results/${prefix}.DEGli.t$i\
        $TFliFile\
        ${resD}/intermediate_results/${prefix}.TF_rank.t$i.trim ${gSet} ${resD}/TG $i &
done; wait

# Final Result : Networks are comprised of resulting TFs/TGs
python makeTGDesc.py ${resD}

# cd ${resD}
# awk '{print}' *.trim|sort|uniq > ../propanet_${cond}_TFs
# awk '{print}' TG/*.TG|sort|uniq > ../propanet_${cond}_TGs
# cd ..
# cat propanet_${cond}_TFs propanet_${cond}_TGs|sort|uniq >propanet_${cond}_TFTGs

# # Exit with an explicit exit status.
# exit 0

#!/bin/bash

# <Replace with the description and/or purpose of this shell script.>
expFile=$1 # "data/DEG.AtGenExpress.signed_zstats.heat_shoots"
seedFile=$2 # "data/DEG.AtGenExpress.signed_binary.heat_shoots"
TFliFile=$3 # "data/Ath_TF_list.gene"
gSet=$4 # (optional) additional gene list file if the user wants
resD="result"

if [ -z ${expFile} ]; then
    expFile="data/DEG.AtGenExpress.signed_zstats.heat_shoots"
fi
if [ -z ${seedFile} ]; then
    seedFile="data/DEG.AtGenExpress.signed_binary.heat_shoots"
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
python2.7 network_weight.py -nwk data/templateNetwork -exp ${expFile} -o ${resD}/intermediate_results/${prefix}.templateNetwork || exit 1

# Extract optimal TF list
python2.7 TF_adding_NP_noCtrl.py ${TFliFile} ${resD}/intermediate_results/${prefix}.templateNetwork ${expFile} ${seedFile} -p 10 -cond ${prefix} -outD ${resD}/intermediate_results || exit 1

# Extract target genes of the optimal TFs
for ((i=1;i<=$[$n-1];i++));do
    python2.7 Target_genes.py ${resD}/intermediate_results/${prefix}.nwk.t$i ${resD}/intermediate_results/${prefix}.DEGli.t$i $TFliFile ${resD}/intermediate_results/${prefix}.TF_rank.t$i.trim ${gSet} ${resD}/TG $i&
done; wait

# Final Result : Networks are comprised of resulting TFs/TGs
python2.7 makeTGDesc.py ${resD}

# cd ${resD}
# awk '{print}' *.trim|sort|uniq > ../propanet_${cond}_TFs
# awk '{print}' TG/*.TG|sort|uniq > ../propanet_${cond}_TGs
# cd ..
# cat propanet_${cond}_TFs propanet_${cond}_TGs|sort|uniq >propanet_${cond}_TFTGs

# # Exit with an explicit exit status.
# exit 0

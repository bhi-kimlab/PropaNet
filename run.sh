#!/bin/bash

# <Replace with the description and/or purpose of this shell script.>
TFliFile="data/Ath_TF_list.gene"
data=$1 # AtGenExpress
cond=$2 # heat_shoots
gSet=$3 # (optional) additional gene list file if the user wants
resD="result/TFTG_list_per_timepoint"

seedFile="data/DEG.$data.signed_binary.$cond"
expFile="data/DEG.$data.signed_zstats.$cond"
n=$(head -n 1 $seedFile| awk '{print NF}')

# Weighted template network construction
python network_weight.py -nwk data/templateNetwork -exp ${expFile} -o ${resD}/templateNetwork.heat_shoots

# Extract optimal TF list
python TF_adding_NP_noCtrl.py ${TFliFile} ${resD}/templateNetwork.heat_shoots ${expFile} ${seedFile} -p 10 -cond ${data}.${cond} -outD ${resD}

# Extract target genes of the optimal TFs
for ((i=1;i<=$[$n-1];i++));do
	echo "===========$i=============="
    python Target_genes.py ${resD}/$data.$cond.nwk.t$i ${resD}/$data.$cond.DEGli.t$i $TFliFile ${resD}/$data.$cond.TF_rank.t$i.trim $gSet $data.$cond.t$i ${resD}/TG $i&   
done; wait

cd ${resD}
awk '{print}' *.trim|sort|uniq > ../propanet_${cond}_TFs
awk '{print}' TG/*.TG|sort|uniq > ../propanet_${cond}_TGs
cd ..
cat propanet_${cond}_TFs propanet_${cond}_TGs|sort|uniq >propanet_${cond}_TFTGs

# Exit with an explicit exit status.
exit 0

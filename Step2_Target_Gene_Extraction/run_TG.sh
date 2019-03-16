#!/bin/bash

out=$1
seedFile="exp/exp.$out.DEG.zvalue"
#out="cold.D2"
n=$(head -n 1 $seedFile| awk '{print NF}')

for ((i=1;i<=$[$n-1];i++));do
	echo "===========$i=============="
	/package/anaconda_new/bin/python2.7 Target_genes.py "network/templateNetwork.$out.corr" "exp/DEG_per_tp/all.degs.$out.t$i" "exp/res_DETF/exp.$out.DEG.zvalue.$i.trim" "$out" "$i"
done




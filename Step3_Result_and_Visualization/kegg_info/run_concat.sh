#!/bin/bash
flag=$1
cntTP=$1

TGsetli=$(ls $flag/*.AT*)

for f in $TGsetli; do
	python kegginfo_to_TGlist.py $f -o ".out"
done	

for ((i=1;i<=cntTP;i++)); do
	paste *.$i > "TG_w_symbol.$i"
done

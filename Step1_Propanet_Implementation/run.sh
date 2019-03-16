#!/bin/bash
#Given condition specific time-series data (z-value, binaryFile) from limma: exp.heat.D2.DEG.zvalue, exp.heat.D2.DEG.binary
nwk=$1
TF=$2
DEG.zval=$3
DEG.bin=$4
out=$5

nwk.out.corr=$nwk.$out.corr

#STEP1. calculate network edge score with correlation
    python network_weight.py -nwk $nwk -exp $DEG.zval -p 15 -o $nwk.out.corr 

#STEP2. Influence maximization to prioritize TF
#STEP3. Network propagation by adding TF one by one (exclue if correlation score decreases)
    python TF_adding_NP.py $TF $nwk.out.corr $DEG.zval $DEG.bin -p 10 -o $out 

##==> Retrieve master regulatory TF for each timepoint

#STEP4(optional). TGs for each TF, GO enrichment
    #python DEGnet.py

import os
import numpy as np
import argparse
import networkx as nx

def makeGeneSet(DEGFile,geneSetF,outD):
        DEGSet = set(np.genfromtxt(DEGFile, dtype=np.str))
        if geneSetF=="None" : geneSet = DEGSet
        else : geneSet = set(np.genfromtxt(geneSetF, dtype=np.str))
        return set(DEGSet)&geneSet

def Target_genes_noTF(TF, g, DEGli, TFSet) :
    visited = set()
    res=set()
    DEGset=set(DEGli)
    stack = []
    stack.append(TF)
    while len(stack)>0:
        v = stack.pop()
        if v not in visited :
            visited.add(v)
            if v not in TFSet:
                res.add(v)
            for s in g.successors(v) :
                if s not in visited:
                    stack.append(s)
    return res

def Target_genes(TF, g, DEGli, TFSet, geneSet) :
    visited = set()
    res=set()
    DEGset=set(DEGli)
    coverTGs = Target_genes_noTF(TF, g, DEGli, TFSet)
    allEdges = False if len(coverTGs&geneSet)/float(len(geneSet))>0.5 else True
    stack = []
    stack.append(TF)
    subnet = nx.DiGraph()
    while len(stack)>0:
        v = stack.pop()
        if v not in visited :
            visited.add(v)
            if v in DEGset and v not in TFSet:
                res.add(v)
            for s in g.successors(v) :
                if s not in visited and (allEdges or g[v][s]['weight']>0.8):
                    subnet.add_edge(v, s, weight=g[v][s]['weight'])
                    stack.append(s)
        # include non-DE TFs in a shortest path w/ the highest weight
        edgeSet = set()
        midTFs = set()
        for v in res :
            currPath = []
            currPathScore = -np.inf
            for p in nx.all_shortest_paths(subnet, source=TF, target=v) :
                pscore = 0
                for i in range(len(p)-1) :
                    pscore += g[p[i]][p[i+1]]['weight']
                pscore = pscore / float(len(p)-1)
                if pscore > currPathScore : 
                    currPathScore = pscore
                    currPath = p
            midTFs = midTFs | set(currPath)
            for i in range(len(currPath)-1):
                edgeSet.add((currPath[i], currPath[i+1], subnet[currPath[i]][currPath[i+1]]['weight']))
    return edgeSet, (res|midTFs)-set([TF]), res


if __name__=="__main__":
    parser = argparse.ArgumentParser(usage='python %(prog)s network/templateNetwork.txt exp/DEG_per_tp/all.degs.cold.D2.t1 TF/Ath_TF_list.gene  outDirectory/exp.cold.D2.DEG.zvalue.1.trim cold.D2.t1 outDirectory')
    parser.add_argument('TFliFile',help='TF list File')
    parser.add_argument('nwkFile',help='Network file')
    parser.add_argument('DEGliFile',help='DEG list file for n-th timepoint')
    parser.add_argument('TFtrimmedFile',help='TF result file from NP')
    parser.add_argument('-geneSet')
    parser.add_argument('-out')
    parser.add_argument('-outD')
    parser.add_argument('-n',type=int,help='nth time point')
    args=parser.parse_args()
    if not os.path.exists(args.outD) : os.mkdir(args.outD)
    exp = pd.read_csv(args.expFile, sep='\t',index_col=0)
    with open(args.TFliFile) as tfF:
        TFli=tfF.read().strip().split()
    TFset=set(TFli)
    i=args.n
    network = nx.read_edgelist(args.nwkFile,data=(('weight',float),),create_using=nx.DiGraph())
    DEGli = set(np.genfromtxt(args.geneSet, dtype=np.str)) 
    if args.geneSet != None: geneSet = set(np.genfromtxt(args.geneSet, dtype=np.str))&DEGli 
    else : geneSet = set(DEGli)
    
    #network = network.subgraph(DEGli)
    edgeSet = set()
    for idx,TF in enumerate(geneSet):
        edges, TGenes, TGs =Target_genes(TF,network,DEGli,TFSet,geneSet)
        with open('{}/{}.{}.{}.TG'.format(args.outD,args.out,idx,TF),'w') as f3:
            f3.write('\n'.join(list(TGenes)))
            edgeSet |= edges
        with open('{}/subnetwork.{}'.format(args.outD,args.n),'w') as f4:
            for (a, b, c) in edgeSet : f4.write(a+'\t'+b+'\t'+str(c)+'\n')
                    

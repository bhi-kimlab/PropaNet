import pandas as pd
import os
import numpy as np
import argparse
import networkx as nx

def makeGeneSet(DEGFile,geneSetF,n,outD):
    ##-control
    DEGSet = set(np.genfromtxt(DEGFile, dtype=np.str))
    if geneSetF==None : geneSet = DEGSet
    else : geneSet = set(np.genfromtxt(geneSetF, dtype=np.str))
    return set(DEGSet)&geneSet

def Target_genes_noTF(TF, g, TFSet) :
    visited = set()
    res=set()
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
    coverTGs = Target_genes_noTF(TF, g, TFSet)
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
    parser = argparse.ArgumentParser()
    parser.add_argument('nwkFile',help='network file')
    parser.add_argument('DEGfile')
    parser.add_argument('TFliFile')
    parser.add_argument('TFfile')
    parser.add_argument('-geneSetFile')
    parser.add_argument('out')
    parser.add_argument('outD')
    parser.add_argument('n')
    args=parser.parse_args()
    if not os.path.exists(args.outD) : os.mkdir(args.outD)
    geneSet = makeGeneSet(args.DEGfile, args.geneSetFile, args.n, args.outD)
    network = nx.read_edgelist(args.nwkFile,data=(('weight',float),),create_using=nx.DiGraph())
    TFSet = set(np.genfromtxt(args.TFliFile, dtype=np.str, delimiter='\t'))
    with open(args.DEGfile) as f, open(args.TFfile) as f2:
        DEGli=f.read().strip().split()
        TFli=f2.read().strip().split()
        print args.TFfile, TFli
    #network = network.subgraph(DEGli)
        edgeSet = set()
    for idx,TF in enumerate(TFli):
        edges, TGenes, TGs =Target_genes(TF,network,DEGli,TFSet,geneSet)
        with open('{}/{}.{}.{}.TG'.format(args.outD,args.out,idx,TF),'w') as f3:
            f3.write('\n'.join(list(TGenes)))
            edgeSet |= edges
        with open('{}/subnetwork.{}'.format(args.outD,args.n),'w') as f4:
            for (a, b, c) in edgeSet : f4.write(a+'\t'+b+'\t'+str(c)+'\n')
                    

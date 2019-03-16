import networkx as nx
import argparse
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('nwkFile')
parser.add_argument('-bin',default=None)#
parser.add_argument('-TF',default=None)#
parser.add_argument('geneSetFile')
parser.add_argument('out')
parser.add_argument('-n', help='timepoint',default=None)#
args=parser.parse_args()

network = nx.read_edgelist(args.nwkFile,data=(('weight',float),),create_using=nx.DiGraph())
geneSet=pd.read_csv(args.geneSetFile,header=None,names=['gene'])
#without considering timepoint --> edge	
if args.n == args.bin == args.TF == None:	
	network = network.subgraph(geneSet['gene'].tolist())
	#nx.write_edgelist(network,'edge.'+args.out+'.t'+args.n,data=['weight'])
	nx.write_edgelist(network,'edge.'+args.out,data=['weight'])


#considering timepoint	--> node
elif args.n != None and args.bin != None and args.TF != None:	
	with open(args.TF) as f:
		TFli=f.read().strip().split()
	geneSet=geneSet.set_index('gene')
	bins=pd.read_csv(args.bin,sep='\t').iloc[:,[0,44+int(args.n)]].set_index('gene')
	
	geneSet=pd.concat([bins,geneSet],join='inner',axis=1).reset_index()
	geneSet.columns=['gene','DEG']

	geneSet['TF']=0
	geneSet.loc[geneSet['gene'].isin(TFli),'TF']=1
	geneSet.columns=['t{}_gene'.format(args.n), 't{}_DEG'.format(args.n),'t{}_TF'.format(args.n)]
	geneSet[['t{}_gene'.format(args.n), 't{}_DEG'.format(args.n),'t{}_TF'.format(args.n)]].to_csv('node.'+args.out+'.t'+args.n,sep='\t',index=None)

else:
	print('Something went wrong!!')

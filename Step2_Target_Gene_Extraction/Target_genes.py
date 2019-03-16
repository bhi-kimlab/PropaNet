import argparse
import networkx as nx

def Target_genes(TF, g, DEGli) :
	visited = set()
	res=set()
	DEGset=set(DEGli)
	stack = []
	stack.append(TF)
	while len(stack)>0:
		v = stack.pop()
		if v not in visited :
			visited.add(v)
			if v in DEGset:
				res.add(v)
			for s in g.successors(v) :
				if s not in visited:
				#if s not in visited and g[v][s]['weight']>0.5:
					stack.append(s)
	return res

if __name__=="__main__":
	parser = argparse.ArgumentParser(usage='python %(prog)s network/templateNetwork.cold.D2.corr exp/DEG_per_tp/all.degs.cold.D2.t1 exp/DEG_TFlist/deg.TFlist.cold.D2.t1 res_DETF/cold.D2.t1')
	parser.add_argument('nwkFile',help='network file')
	parser.add_argument('DEGfile')
	parser.add_argument('TFfile')
	parser.add_argument('out')
	args=parser.parse_args()
	network = nx.read_edgelist(args.nwkFile,data=(('weight',float),),create_using=nx.DiGraph())
	with open(args.DEGfile) as f,open(args.TFfile) as f2:
		DEGli=f.read().strip().split()
		TFli=f2.read().strip().split()[:10]
	for idx,TF in enumerate(TFli):
		TGenes=Target_genes(TF,network,DEGli)
		with open('{}.{}.{}.TG'.format(args.out,idx,TF),'w') as f3:
			f3.write('\n'.join(list(TGenes)))	

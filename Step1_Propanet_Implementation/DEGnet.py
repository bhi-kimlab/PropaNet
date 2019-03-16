
flist = ['cold.D2','cold.D3','heat.D2','heat.D3']

for f in flist :
  oriNet = open('templateNetwork.'+f+'.corr', 'r').readlines()
  l = 6 if f[:4]=='cold' else 4
  for i in range(1, l+1) :
    outf = open('templateNetwork.DEG.'+f+'.'+str(i), 'w')
    DEGset = set(open('../exp/DEG_per_tp/all.degs.'+f+'.t'+str(i), 'r').read().strip().split())
    for line in oriNet :
      tp = line.split()
      if tp[0] in DEGset and tp[1] in DEGset : outf.write(line)

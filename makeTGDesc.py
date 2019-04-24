import pandas as pd
import argparse
import sys
import os

def findDETG(resPath, cond) : # cond =='cold'
  bins = pd.read_csv('exp/exp.'+cond+'.DEG.binary', sep='\t')
  for f in os.listdir(resPath) :
    if f[-4:]!='trim': continue
    tp = f.split('.')[5]
    TFSet = set(open(resPath+'/'+f, 'r').read().strip().split()[2::2])
    DEGSet = set(bins['gene'][bins.iloc[:,44+int(tp)]!=0])
    print tp, len(DEGSet&TFSet), DEGSet&TFSet
  return   

def makeDict(keggF) :
  symDict = {}
  descDict = {}
  for line in open(keggF, 'r').readlines():
    tp = line.split('\t')
    kid = tp[0].split(':')[1]
    sd = tp[1].split(';')
    sym = '' if len(sd)==1 else sd[0].split(', ')[0]
    desc = sd[0] if len(sd)==1 else sd[1]
    symDict[kid]=sym.rstrip()
    descDict[kid]=desc.rstrip()
  return symDict, descDict

def main_TG():
  parser = argparse.ArgumentParser(usage='python %(prog)s resPath')
  parser.add_argument('resPath')
  args=parser.parse_args()
  symDict, descDict = makeDict('data/kegg_ath_gene.txt')
  outF = open(args.resPath+'/TGdesc.txt', 'w')
  for f in os.listdir(args.resPath) :
    if os.path.isdir(args.resPath+'/'+f) or f[-2:]!='TG': continue
    tp = f.split('.')[2].upper()
    TF = f.split('.')[4]
    for line in open(args.resPath+'/'+f, 'r').readlines():
      TG = line.rstrip()
      TFsym = '' if TF not in symDict else symDict[TF]
      TFdesc = '' if TF not in descDict else descDict[TF]
      if TG not in symDict : outF.write('{}\t{}\t{}\t{}\t{}\t\t\n'.format(tp, TF, TFsym, TFdesc, TG))
      else : outF.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(tp, TF, TFsym, TFdesc, TG, symDict[TG], descDict[TG]))

def main():
  parser = argparse.ArgumentParser(usage='python %(prog)s resPath')
  parser.add_argument('resPath')
  args=parser.parse_args()
  symDict, descDict = makeDict('data/kegg_ath_gene.txt')
  outF = open(args.resPath+'/network_edges.txt', 'w')
  for f in os.listdir(args.resPath+'/TG') :
    if os.path.isdir(args.resPath+'/TG/'+f) or f[-2:]=='TG' or f[-3:]=='txt': continue
    print f
    tp = 'T'+f.split('.')[1]
    for line in open(args.resPath+'/TG/'+f, 'r').readlines():
      temp = line.split('\t')
      source = temp[0]
      target = temp[1]
      TFsym = '' if source not in symDict else symDict[source]
      TFdesc = '' if source not in descDict else descDict[source]
      if target not in symDict : outF.write('{}\t{}\t{}\t{}\t{}\t\t\n'.format(tp, source, TFsym, TFdesc, target))
      else : outF.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(tp, source, TFsym, TFdesc, target, symDict[target], descDict[target]))

if __name__=='__main__':
        main()

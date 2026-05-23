import pandas as pd, os
os.makedirs("/mnt/data/gold_index",exist_ok=True)
rg=pd.read_csv("/mnt/data/ref/hg38/refGene.txt",sep="\t",header=None,
 names=['bin','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score','name2','s1','s2','frames'])
rg=rg[rg.chrom.str.match(r'^chr(\d+|[XY])$')]; rg=rg[rg.cdsStart<rg.cdsEnd].copy()
rg['span']=rg.txEnd-rg.txStart
longest=rg.sort_values('span').groupby('name2').tail(1)
iv=set(); genes=set()
for r in longest.itertuples():
    ss=[int(x) for x in r.exonStarts.rstrip(',').split(',') if x]
    ee=[int(x) for x in r.exonEnds.rstrip(',').split(',') if x]
    for a,b in zip(ss,ee):
        a2,b2=max(a,r.cdsStart),min(b,r.cdsEnd)
        if b2>a2:
            iv.add((r.chrom,a2,b2)); genes.add((r.chrom,a2,b2,r.name2,r.strand))
W="/mnt/data/gold_index"
with open(W+"/cds.bed","w") as f:
    for c,a,b in sorted(iv): f.write(f"{c}\t{a}\t{b}\n")
with open(W+"/cds_genes.tsv","w") as f:
    for c,a,b,g,st in sorted(genes): f.write(f"{c}\t{a}\t{b}\t{g}\t{st}\n")
print("cds intervals",len(iv),"gene-intervals",len(genes))

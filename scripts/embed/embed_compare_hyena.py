"""Reproducibility: HyenaDNA embedding for the F1 full-vs-simple comparison (run on GPU node).
Reads /tmp/compare_sites.csv (key,src,label,seq[513bp C-oriented]); writes mean-pool and
center-token (+-16) embeddings. Used with transformers==4.39.3 (pinned; 5.x breaks HyenaDNA remote code).
Run: PYTHONPATH=/tmp/tlibs python embed_compare_hyena.py [mean|center]"""
import csv, torch, time, sys
csv.field_size_limit(10**7)
from transformers import AutoModel, AutoTokenizer
MODEL="LongSafari/hyenadna-small-32k-seqlen-hf"
MODE=sys.argv[1] if len(sys.argv)>1 else "mean"; R=16
keys=[];seqs=[]
for row in csv.DictReader(open('/tmp/compare_sites.csv')):
    keys.append(row['key']);seqs.append(row['seq'])
print(f"embedding {len(seqs)} ({MODE})",flush=True)
tok=AutoTokenizer.from_pretrained(MODEL,trust_remote_code=True)
model=AutoModel.from_pretrained(MODEL,trust_remote_code=True).to('cuda').eval()
emb={};B=64;t0=time.time()
for i in range(0,len(seqs),B):
    bk=keys[i:i+B];bs=seqs[i:i+B]
    enc=tok(bs,return_tensors='pt',padding=True,truncation=True);enc={k:v.to('cuda') for k,v in enc.items()}
    with torch.no_grad():
        out=model(**enc);h=out.last_hidden_state if hasattr(out,'last_hidden_state') else out[0]
    if MODE=="center":
        L=h.shape[1];c=L//2;e=h[:,max(0,c-R):c+R+1,:].mean(1)
    else:
        e=h.mean(1)
    e=e.detach().cpu().float()
    for k,v in zip(bk,e): emb[k]=v.clone()
out=f"/tmp/compare_hyena{'_center' if MODE=='center' else ''}.pt"
torch.save(emb,out); print(f"DONE {len(emb)} -> {out} in {time.time()-t0:.0f}s",flush=True)

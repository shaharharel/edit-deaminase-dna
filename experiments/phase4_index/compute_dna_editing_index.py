"""DNA Editing Index per gene (Levanon-style, adapted for DNA WGS).

Inputs: parsed-pileup parquets for treated (BE4_clone1) and control (Parent_WGS),
each with cols (chrom,pos,strand,motif,gene,depth,n_A,n_C,n_G,n_T,vaf_main,vaf_noise,sample).

Three sequential controls:
 1) low-VAF filter: at each position, drop reads from positions with overall VAF_main > VAF_CAP
    (germline/clonal-fixed mutations). This is per-position, not per-read.
 2) within-sample mismatch-direction noise floor: subtract noise direction.
 3) treated minus control (Parent): final differential.

Plus motif-negative control (non-TpC C's) — the signal should disappear there.

Outputs:
 - per_gene_index.parquet : per-gene DEI in treated, control, differential; per-gene N_TpC, N_npc
 - qc_summary.txt          : motif-spectrum, VAF distribution, mismatch-direction sanity
"""
import sys, pandas as pd, numpy as np

TREATED = sys.argv[1]
CONTROL = sys.argv[2]
OUT_PARQ = sys.argv[3]
OUT_QC = sys.argv[4]
VAF_CAP = float(sys.argv[5]) if len(sys.argv) > 5 else 0.40   # drop positions with VAF > 0.40
MIN_TPC_PER_GENE = int(sys.argv[6]) if len(sys.argv) > 6 else 20

T = pd.read_parquet(TREATED)
C = pd.read_parquet(CONTROL)

# tag samples
T['sample'] = 'treated'; C['sample'] = 'control'

def low_vaf(df, cap):
    """Drop positions whose overall VAF_main exceeds the cap (germline/clonal). Keep all reads at
    surviving positions. This produces a per-position pass/fail mask."""
    return df[df.vaf_main <= cap].copy()

def aggregate(df, motif):
    """Per-gene aggregate after VAF filter. Sum variant-direction reads / sum total reads."""
    sub = df[df.motif == motif]
    by_gene = sub.groupby('gene').agg(
        n_pos = ('pos','count'),
        depth_total = ('depth','sum'),
        var_main = ('vaf_main', lambda x: int(np.round((x.values * sub.loc[x.index,'depth'].values).sum()))),
        var_noise = ('vaf_noise', lambda x: int(np.round((x.values * sub.loc[x.index,'depth'].values).sum()))),
    )
    by_gene['EI_main'] = by_gene.var_main / by_gene.depth_total.clip(lower=1)
    by_gene['EI_noise'] = by_gene.var_noise / by_gene.depth_total.clip(lower=1)
    by_gene['EI_net'] = by_gene.EI_main - by_gene.EI_noise
    return by_gene

def run(df, motif, label):
    pre_n = len(df[df.motif == motif])
    df_f = low_vaf(df[df.motif == motif], VAF_CAP)
    post_n = len(df_f)
    agg = aggregate(pd.concat([df_f, df[df.motif != motif]], ignore_index=True), motif)
    return agg, pre_n, post_n

T_tpc, t_pre, t_post = run(T, 'tpc', 'treated_tpc')
T_npc, _, _         = run(T, 'npc', 'treated_npc')
C_tpc, c_pre, c_post = run(C, 'tpc', 'control_tpc')
C_npc, _, _         = run(C, 'npc', 'control_npc')

# merge per-gene
out = T_tpc[['n_pos','depth_total','EI_main','EI_noise','EI_net']].rename(columns=lambda c: f'tT_{c}')
out = out.join(C_tpc[['n_pos','depth_total','EI_main','EI_noise','EI_net']].rename(columns=lambda c: f'tC_{c}'), how='outer')
out = out.join(T_npc[['n_pos','depth_total','EI_main','EI_noise','EI_net']].rename(columns=lambda c: f'nT_{c}'), how='outer')
out = out.join(C_npc[['n_pos','depth_total','EI_main','EI_noise','EI_net']].rename(columns=lambda c: f'nC_{c}'), how='outer')
out = out.fillna({c:0 for c in out.columns})

# The DNA Editing Index (final): differential of net (within-sample mismatch-direction corrected),
# treated minus control, at TpC sites.
out['DEI_tpc'] = out.tT_EI_net - out.tC_EI_net
# motif-negative control (should be near zero)
out['DEI_npc'] = out.nT_EI_net - out.nC_EI_net
# motif-specificity ratio per gene (informative for QC; clipped to avoid division noise)
out['motif_specificity'] = out.DEI_tpc / (out.DEI_npc.abs() + 1e-6)

# coverage filter
out = out[(out.tT_n_pos >= MIN_TPC_PER_GENE) & (out.tC_n_pos >= MIN_TPC_PER_GENE)].copy()
out = out.reset_index().rename(columns={'index':'gene'})
out.to_parquet(OUT_PARQ)

# QC summary
def summ(arr): return f'n={len(arr)} mean={arr.mean():.5f} med={arr.median():.5f} q90={arr.quantile(0.9):.5f}'
lines = [
    f"=== DNA Editing Index QC ===",
    f"VAF_CAP: {VAF_CAP}   MIN_TPC_PER_GENE: {MIN_TPC_PER_GENE}",
    f"TpC positions, treated: pre-filter {t_pre:,}  post-filter {t_post:,} ({100*t_post/max(1,t_pre):.1f}%)",
    f"TpC positions, control: pre-filter {c_pre:,}  post-filter {c_post:,} ({100*c_post/max(1,c_pre):.1f}%)",
    f"",
    f"--- Per-gene DEI distributions ---",
    f"Treated EI_net at TpC (within-sample C->T minus noise): {summ(out.tT_EI_net)}",
    f"Control EI_net at TpC:                                  {summ(out.tC_EI_net)}",
    f"Treated EI_net at non-TpC (motif-negative):             {summ(out.nT_EI_net)}",
    f"",
    f"--- Differential DEI (the headline) ---",
    f"DEI at TpC (treated - control):    {summ(out.DEI_tpc)}",
    f"DEI at non-TpC (motif-negative):   {summ(out.DEI_npc)}",
    f"",
    f"--- Motif-specificity sanity ---",
    f"Fraction of genes with DEI_tpc > DEI_npc: {(out.DEI_tpc > out.DEI_npc).mean():.3f}",
    f"Mean(DEI_tpc) / Mean(DEI_npc) ratio:      {out.DEI_tpc.mean()/max(out.DEI_npc.mean(),1e-9):.2f}",
    f"",
    f"PASS criteria (Levanon-style, adapted):",
    f"  [a] DEI_tpc > 0 and significantly above DEI_npc (motif specificity)",
    f"  [b] DEI_tpc > DEI_npc in >70% of genes",
    f"  [c] No huge mass at VAF=1.0 left after filter (germline removed)",
]
with open(OUT_QC,'w') as fh: fh.write('\n'.join(lines))
print('\n'.join(lines))

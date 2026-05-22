"""F3: per-gene intrinsic rAPOBEC1 (BE4) susceptibility ranking + cancer-driver tiering.
The safety-gate prototype output: which cancer-driver genes are intrinsically most vulnerable
to this editor's chemistry -> prioritize for off-target screening. HONEST: intrinsic (sequence/motif)
prior, accessibility-agnostic, UNCALIBRATED (no continuous DNA gold; calibrated burden is >=1Mb only)."""
import numpy as np, pandas as pd, math
CANCER=set("""ABL1 AKT1 ALK AMER1 APC AR ARID1A ASXL1 ATM ATRX AXIN1 BAP1 BCL2 BCOR BRAF BRCA1 BRCA2 CARD11 CBL CDC73 CDH1 CDKN2A CEBPA CIC CREBBP CTNNB1 DAXX DNMT3A EGFR EP300 ERBB2 EZH2 FBXW7 FGFR2 FGFR3 FLT3 FOXL2 GATA3 GNA11 GNAQ GNAS HNF1A HRAS IDH1 IDH2 JAK1 JAK2 JAK3 KDM5C KDM6A KIT KMT2D KRAS MAP2K1 MAP3K1 MED12 MEN1 MET MLH1 MPL MSH2 MSH6 MYD88 NF1 NF2 NFE2L2 NOTCH1 NOTCH2 NPM1 NRAS PAX5 PBRM1 PDGFRA PIK3CA PIK3R1 PPP2R1A PRDM1 PTCH1 PTEN PTPN11 RB1 RET RNF43 SETD2 SF3B1 SMAD2 SMAD4 SMARCA4 SMARCB1 SMO SOCS1 SPOP STAG2 STK11 TET2 TNFAIP3 TP53 TSC1 U2AF1 VHL WT1 AKT2 ARID2 ATR B2M BARD1 BCL6 BRD4 BRIP1 BTK CALR CASP8 CBFB CCND1 CCND2 CCND3 CCNE1 CD274 CD79A CD79B CDK12 CDK4 CDK6 CDKN1B CDKN2C CHEK2 CRLF2 CSF1R CSF3R CTCF CXCR4 DDR2 ERBB3 ERBB4 ERG ESR1 ETV6 FANCA FANCC FGFR1 FGFR4 FLCN FUBP1 GATA1 GATA2 HGF IKZF1 IRF4 JUN KDM5A KDR KEAP1 KMT2A KMT2C MAP2K2 MAP2K4 MAPK1 MDM2 MDM4 MITF MTOR MUTYH MYC MYCL MYCN NKX2-1 NSD2 NSD3 NTRK1 NTRK2 NTRK3 PALB2 PDGFRB PHF6 PIM1 PRKAR1A RAD21 RAF1 RARA ROS1 RUNX1 SDHA SDHB SDHC SDHD SOX2 SPEN SRC SRSF2 STAT3 SUFU SYK TGFBR2 TMPRSS2 TNFRSF14 TSC2 XPO1 AKT3 ARAF ARID1B AURKA AURKB AXL BCL10 BCORL1 BCR BIRC3 BLM BTG1 CDK8 CDKN2B CHEK1 CRKL CYLD DOT1L EED EIF4A2 EPHA3 EPHB1 ERCC4 ETV1 FAS FGF19 FGF3 FGF4 FH FLT1 FLT4 FOXO1 FOXP1 GRIN2A GSK3B ID3 IGF1R IKBKE IL7R INPP4B IRS2 KLF4 LMO1 MALT1 MAP3K13 MCL1 MEF2B MRE11 MSH3 MSI2 NBN NCOR1 NFKBIA NSD1 NT5C2 P2RY8 PDCD1 PIK3CB PMS2 POLD1 POLE POT1""".split())
df=pd.read_parquet('data/processed/pred_gene_burden_v1.parquet')
df['cancer']=df.gene.isin(CANCER)
df=df.sort_values('pred_burden',ascending=False).reset_index(drop=True)
df['rank']=np.arange(1,len(df)+1); df['pct']=df['rank']/len(df)
print(f"genes ranked: {len(df)} | cancer-driver genes present: {df.cancer.sum()}")
# enrichment of cancer genes in top decile
N=len(df); K=df.cancer.sum(); n=int(0.1*N); ov=int(df.head(n).cancer.sum()); exp=K*n/N
hp=sum(math.comb(K,i)*math.comb(N-K,n-i) for i in range(ov,min(K,n)+1))/math.comb(N,n) if ov>0 else 1
print(f"cancer genes in top-10% susceptibility: {ov} (expected {exp:.1f}, {ov/max(exp,1):.2f}x, HG p={hp:.2e})")
print(f"median susceptibility: cancer={df[df.cancer].pred_burden.median():.3f} vs other={df[~df.cancer].pred_burden.median():.3f}")
print("\n=== TOP 20 cancer-driver genes by intrinsic rAPOBEC1(BE4) susceptibility (screening priority) ===")
top=df[df.cancer].head(20)[['rank','gene','pred_burden','n_C_scored','gc_mean']]
print(top.to_string(index=False))
df.to_parquet('data/processed/risk_ranking_v1.parquet',index=False)
print("\nwrote data/processed/risk_ranking_v1.parquet")
print("[honest] intrinsic sequence/motif prior, accessibility-agnostic, UNCALIBRATED.")
print("Use = prioritize these genes for experimental off-target screening for this editor.")

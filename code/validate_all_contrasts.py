#!/usr/bin/env python3
"""Full end-to-end validation: re-run GSEA from all 37 per-contrast DE files
and compare NES against the delivered nes_matrix.csv. Writes a per-contrast
validation table."""
import pandas as pd, numpy as np, gseapy, warnings, os
warnings.filterwarnings('ignore')

nes_ref = pd.read_csv('/mnt/user-uploads/nes_matrix.csv', index_col=0).dropna(how='all')
nes_ref = nes_ref[nes_ref.index.notna()]

gmt = {}
with open('/mnt/user-uploads/curated_gene_sets.gmt') as f:
    for line in f:
        parts = line.strip().split('\t')
        gmt[parts[0]] = [g for g in parts[2:] if g]

# explicit DE file stem -> nes_matrix column mapping (all 37, verified 1:1)
STEM_TO_COL = {
    "NYGC_Cervical_Cord": "NYGC_Cervical",
    "NYGC_Lumbar_Cord": "NYGC_Lumbar",
    "NYGC_Motor_Cortex": "NYGC_MotorCortex",
    "OSD863_noMicro_FvG": "OSD863_healthy_noMicro_FvG",
    "OSD863_withMicro_FvG": "OSD863_healthy_withMicro_FvG",
    "OSD863_Flight_MicroVsNo": "OSD863_healthy_Flight_MicroVsNo",
    "OSD863_Ground_MicroVsNo": "OSD863_healthy_Ground_MicroVsNo",
    "OSD871_noMicro_FvG": "OSD871_healthy_noMicro_FvG",
    "OSD871_withMicro_FvG": "OSD871_healthy_withMicro_FvG",
    "OSD871_Flight_MicroVsNo": "OSD871_healthy_Flight_MicroVsNo",
    "OSD871_Ground_MicroVsNo": "OSD871_healthy_Ground_MicroVsNo",
    "OSD999_simGCR_50cGy": "OSD-999",
    "OSD202_IR_1mon": "OSD-202_IR_only_1mon_vs_Control_1mon",
    "OSD202_IR_4mon": "OSD-202_IR_only_4mon_vs_Control_4mon",
    "OSD202_HLU_1mon": "OSD-202_HLU_only_1mon_vs_Control_1mon",
    "OSD202_HLU_4mon": "OSD-202_HLU_only_4mon_vs_Control_4mon",
    "OSD202_HLU_IR_1mon": "OSD-202_HLU+IR_1mon_vs_Control_1mon",
    "OSD202_HLU_IR_4mon": "OSD-202_HLU+IR_4mon_vs_Control_4mon",
    "GSE226395_clinostat": "GSE226395_MN_SMG",
    "RRRM2_FL_Young_ISST_FvG": "RRRM2_FL_Young_FLT_vs_GC",
    "RRRM2_FL_Mature_ISST_FvG": "RRRM2_FL_Old_FLT_vs_GC",
    "RRRM2_FL_Age_GC_Old_vs_Young": "RRRM2_FL_Age_GC_Old_vs_Young",
    "RRRM2_FL_Age_FLT_Old_vs_Young": "RRRM2_FL_Age_FLT_Old_vs_Young",
    "RRRM2_CB_Young_ISST_FvG": "RRRM2_CB_Young_ISST_FLT_vs_GC",
    "RRRM2_CB_Young_LAR_FvG": "RRRM2_CB_Young_LAR_FLT_vs_GC",
    "RRRM2_CB_Mature_ISST_FvG": "RRRM2_CB_Old_FLT_vs_GC",
    "RRRM2_CB_Age_GC_Old_vs_Young": "RRRM2_CB_Age_GC_Old_vs_Young",
    "RRRM2_CB_Age_FLT_Old_vs_Young": "RRRM2_CB_Age_FLT_Old_vs_Young",
    "RRRM2_CB_Young_ISST_vs_LAR": "RRRM2_CB_Young_ISST_vs_LAR",
    "RRRM2_HPC_Young_ISST_FvG": "RRRM2_HPC_Young_ISST_FLT_vs_GC",
    "RRRM2_HPC_Young_LAR_FvG": "RRRM2_HPC_Young_LAR_FLT_vs_GC",
    "RRRM2_HPC_Mature_ISST_FvG": "RRRM2_HPC_Old_ISST_FLT_vs_GC",
    "RRRM2_HPC_Mature_LAR_FvG": "RRRM2_HPC_Old_LAR_FLT_vs_GC",
    "RRRM2_HPC_Age_GC_Old_vs_Young": "RRRM2_HPC_Age_GC_Old_vs_Young",
    "RRRM2_HPC_Age_ISST_FLT_Old_vs_Young": "RRRM2_HPC_Age_ISST_FLT_Old_vs_Young",
    "RRRM2_HPC_Mature_ISST_vs_LAR": "RRRM2_HPC_Old_ISST_vs_LAR",
    "RRRM2_HPC_Young_ISST_vs_LAR": "RRRM2_HPC_Young_ISST_vs_LAR",
}

# map DE file stems -> nes_matrix columns
DE_DIR = '/mnt/results/de_raw_files/per_contrast_DE'
files = sorted(f[:-4] for f in os.listdir(DE_DIR) if f.endswith('.csv'))
assert set(files) == set(STEM_TO_COL), (
    f"stem mismatch: {set(files) ^ set(STEM_TO_COL)}")
assert len(set(STEM_TO_COL.values())) == 37
rows = []
for stem in files:
    df = pd.read_csv(f'{DE_DIR}/{stem}.csv')
    rnk = df.dropna(subset=['stat']).groupby('gene')['stat'].mean().reset_index()
    rnk = rnk.sort_values('stat', ascending=False)
    try:
        r = gseapy.prerank(rnk=rnk, gene_sets=gmt, nperm=10000, min_size=1,
                           max_size=500, weight_score=1.0, threads=4,
                           outdir=None, verbose=False)
        got = r.res2d.set_index('Term')['NES'].astype(float)
    except Exception as e:
        rows.append({'contrast_file': stem, 'error': str(e)})
        continue
    # explicit mapping (no fuzzy matching)
    refcol = STEM_TO_COL[stem]
    common = got.index.intersection(nes_ref.index)
    ref = nes_ref.loc[common, refcol]
    g = got.loc[common]
    mask = ref.notna() & g.notna()
    corr = np.corrcoef(ref[mask], g[mask])[0, 1] if mask.sum() > 2 else np.nan
    kd = abs(g.get('T1_KEGG_ALS', np.nan) - nes_ref.loc['T1_KEGG_ALS', refcol]) \
        if 'T1_KEGG_ALS' in got.index else np.nan
    rows.append({'contrast_file': stem, 'ref_column': refcol,
                 'n_sets_compared': int(mask.sum()), 'nes_corr': round(corr, 4),
                 'median_abs_dNES': round(float(np.median(abs(ref[mask] - g[mask]))), 4),
                 'KEGG_ALS_abs_diff': round(float(kd), 4)})

out = pd.DataFrame(rows)
out.to_csv('/mnt/results/canonical_outputs/gsea_revalidation_all37.csv', index=False)
print(out.to_string(index=False))
print(f"\nMedian correlation across contrasts: {out['nes_corr'].median():.5f}")

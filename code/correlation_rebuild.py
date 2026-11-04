"""Rebuild gene-level Wald-statistic correlations (methods doc Section 7).

Spearman rho of DESeq2 Wald statistics between each ALS reference and each
human spaceflight contrast, on shared genes. Permutation p-values with 1,000
permutations (minimum reportable p = 0.001). NOTE: the other machine's run
excluded 806 cell-type marker genes; the T4 marker-set gene membership was
not included in the handoff, so this rebuild runs WITHOUT the exclusion.
Differences from the published values gauge how much the exclusion matters.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata

DE = "/mnt/results/de_raw_files/per_contrast_DE"
REFS = ["NYGC_Cervical_Cord", "NYGC_Lumbar_Cord", "NYGC_Motor_Cortex"]
HUMAN_SF = ["OSD863_noMicro_FvG", "OSD863_withMicro_FvG",
            "OSD863_Flight_MicroVsNo", "OSD863_Ground_MicroVsNo",
            "OSD871_noMicro_FvG", "OSD871_withMicro_FvG",
            "OSD871_Flight_MicroVsNo", "OSD871_Ground_MicroVsNo",
            "OSD999_simGCR_50cGy", "GSE226395_clinostat"]
N_PERM = 1000
rng = np.random.default_rng(42)


def load_stat(stem):
    df = pd.read_csv(f"{DE}/{stem}.csv")
    s = df.dropna(subset=["stat"]).groupby("gene")["stat"].mean()
    return s[s.index.notna()]


stats = {s: load_stat(s) for s in REFS + HUMAN_SF}

rows = []
for ref in REFS:
    r = stats[ref]
    for sf in HUMAN_SF:
        x = stats[sf]
        shared = r.index.intersection(x.index)
        a = rankdata(r.loc[shared].values)
        b = rankdata(x.loc[shared].values)
        rho = np.corrcoef(a, b)[0, 1]
        # permutation: shuffle b's gene labels via index permutation
        n = len(shared)
        count = 0
        for _ in range(N_PERM):
            bp = b[rng.permutation(n)]
            if abs(np.corrcoef(a, bp)[0, 1]) >= abs(rho):
                count += 1
        pval = max((count + 1) / (N_PERM + 1), 0.001)
        rows.append({"als_reference": ref, "spaceflight_contrast": sf,
                     "n_shared_genes": n, "rho": round(rho, 3),
                     "R2": round(rho**2, 3), "p_perm": pval})

out = pd.DataFrame(rows)
out.to_csv("/mnt/results/canonical_outputs/gene_level_correlations.csv",
           index=False)
print(out.to_string(index=False))
print(f"\nrho range: {out['rho'].min():.2f} to {out['rho'].max():.2f} | "
      f"max R2: {out['R2'].max():.3f} | "
      f"pairs with R2<0.03: {(out['R2']<0.03).sum()}/{len(out)}")

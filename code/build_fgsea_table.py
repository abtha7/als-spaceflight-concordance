"""Build the full fgsea output table (155 gene sets x 37 contrasts) from the
NES and padj matrices, with gene set category and contrast metadata attached.

Outputs: /mnt/results/fgsea_full_results.csv
"""
import pandas as pd
import numpy as np

# ---------------------------------------------------------------- load
nes = pd.read_csv("/mnt/user-uploads/nes_matrix.csv").dropna(subset=["pathway"])
padj = pd.read_csv("/mnt/user-uploads/padj_matrix.csv").dropna(subset=["pathway"])

assert nes.shape == (155, 38) and padj.shape == (155, 38), "unexpected matrix shape"
assert list(nes.columns) == list(padj.columns), "column mismatch"
assert (nes["pathway"].tolist() == padj["pathway"].tolist()), "pathway order mismatch"

contrasts = list(nes.columns[1:])
assert len(contrasts) == 37

# --------------------------------------------- gene set category mapping
# T1 sets: disease signatures (sub-groups per t1_pathway_guide.csv)
def t1_subgroup(name: str) -> str:
    if "GSE77702" in name or "GSE77707" in name:
        return "RBP knockdown signature"
    if "GSE52946" in name or "GSE52202" in name:
        return "ALS patient signature"
    return "Disease pathway"

# T2 sets: mechanism categories encoded as T2_<Category>__<set name>
def t2_category(name: str) -> str:
    body = name.split("__", 1)[0].replace("T2_", "", 1)
    return body.replace("_", " ")

rows = []
for _, r in nes.iterrows():
    p = r["pathway"]
    if p.startswith("T1_"):
        tier, cat = "Tier 1", t1_subgroup(p)
    else:
        tier, cat = "Tier 2", t2_category(p)
    rows.append({"pathway": p, "tier": tier, "category": cat})

sets = pd.DataFrame(rows)

# ------------------------------------------------------ contrast metadata
meta = {
    "NYGC_MotorCortex": ("NYGC Motor Cortex (GSE124439)", "Human", "Post-mortem motor cortex", "ALS vs control", "ALS reference"),
    "NYGC_Cervical": ("NYGC Cervical Spinal Cord (Zenodo 6385747)", "Human", "Post-mortem cervical spinal cord", "ALS vs control", "ALS reference"),
    "NYGC_Lumbar": ("NYGC Lumbar Spinal Cord (Zenodo 6385747)", "Human", "Post-mortem lumbar spinal cord", "ALS vs control", "ALS reference"),
    "OSD863_healthy_noMicro_FvG": ("OSD-863 cortical organoids (GLDS-716)", "Human", "Cortical organoids, no microglia", "ISS flight vs ground", "Stressor"),
    "OSD863_healthy_withMicro_FvG": ("OSD-863 cortical organoids (GLDS-716)", "Human", "Cortical organoids, with microglia", "ISS flight vs ground", "Stressor"),
    "OSD863_healthy_Flight_MicroVsNo": ("OSD-863 cortical organoids (GLDS-716)", "Human", "Cortical organoids, in flight", "Microglia vs no microglia", "Modifier"),
    "OSD863_healthy_Ground_MicroVsNo": ("OSD-863 cortical organoids (GLDS-716)", "Human", "Cortical organoids, on ground", "Microglia vs no microglia", "Modifier"),
    "OSD871_healthy_noMicro_FvG": ("OSD-871 dopaminergic organoids (GLDS-720)", "Human", "Dopaminergic organoids, no microglia", "ISS flight vs ground", "Stressor"),
    "OSD871_healthy_withMicro_FvG": ("OSD-871 dopaminergic organoids (GLDS-720)", "Human", "Dopaminergic organoids, with microglia", "ISS flight vs ground", "Stressor"),
    "OSD871_healthy_Flight_MicroVsNo": ("OSD-871 dopaminergic organoids (GLDS-720)", "Human", "Dopaminergic organoids, in flight", "Microglia vs no microglia", "Modifier"),
    "OSD871_healthy_Ground_MicroVsNo": ("OSD-871 dopaminergic organoids (GLDS-720)", "Human", "Dopaminergic organoids, on ground", "Microglia vs no microglia", "Modifier"),
    "OSD-999": ("OSD-999 cortical organoids (GSE309755)", "Human", "Cortical brain organoids", "simGCR 50 cGy 5-ion vs sham", "Stressor"),
    "OSD-202_IR_only_1mon_vs_Control_1mon": ("OSD-202 mouse brain (GLDS-202)", "Mouse", "Whole brain, 1 month", "Ionizing radiation only vs control", "Stressor"),
    "OSD-202_IR_only_4mon_vs_Control_4mon": ("OSD-202 mouse brain (GLDS-202)", "Mouse", "Whole brain, 4 months", "Ionizing radiation only vs control", "Stressor"),
    "OSD-202_HLU_only_1mon_vs_Control_1mon": ("OSD-202 mouse brain (GLDS-202)", "Mouse", "Whole brain, 1 month", "Hindlimb unloading only vs control", "Stressor"),
    "OSD-202_HLU_only_4mon_vs_Control_4mon": ("OSD-202 mouse brain (GLDS-202)", "Mouse", "Whole brain, 4 months", "Hindlimb unloading only vs control", "Stressor"),
    "OSD-202_HLU+IR_1mon_vs_Control_1mon": ("OSD-202 mouse brain (GLDS-202)", "Mouse", "Whole brain, 1 month", "HLU + radiation vs control", "Stressor"),
    "OSD-202_HLU+IR_4mon_vs_Control_4mon": ("OSD-202 mouse brain (GLDS-202)", "Mouse", "Whole brain, 4 months", "HLU + radiation vs control", "Stressor"),
    "GSE226395_MN_SMG": ("GSE226395 motor neurons", "Human", "iPSC spinal motor neurons", "2D clinostat simulated microgravity vs normal gravity", "Stressor"),
}

# RRRM-2 contrasts (18): region, age, condition
rrrm2 = {
    "RRRM2_FL_Young_FLT_vs_GC": ("Frontal lobe", "Young (12 wk)", "ISS-terminal flight vs ground control"),
    "RRRM2_FL_Old_FLT_vs_GC": ("Frontal lobe", "Mature (29 wk)", "ISS-terminal flight vs ground control"),
    "RRRM2_FL_Age_GC_Old_vs_Young": ("Frontal lobe", "Age contrast (GC)", "Ground control old vs young"),
    "RRRM2_FL_Age_FLT_Old_vs_Young": ("Frontal lobe", "Age contrast (FLT)", "Flight old vs young"),
    "RRRM2_CB_Young_ISST_FLT_vs_GC": ("Cerebellum", "Young (12 wk)", "ISS-terminal flight vs ground control"),
    "RRRM2_CB_Young_LAR_FLT_vs_GC": ("Cerebellum", "Young (12 wk)", "LAR flight vs ground control"),
    "RRRM2_CB_Old_FLT_vs_GC": ("Cerebellum", "Mature (29 wk)", "ISS-terminal flight vs ground control"),
    "RRRM2_CB_Age_GC_Old_vs_Young": ("Cerebellum", "Age contrast (GC)", "Ground control old vs young"),
    "RRRM2_CB_Age_FLT_Old_vs_Young": ("Cerebellum", "Age contrast (FLT)", "Flight old vs young"),
    "RRRM2_CB_Young_ISST_vs_LAR": ("Cerebellum", "Young (12 wk)", "ISS-terminal vs LAR"),
    "RRRM2_HPC_Young_ISST_FLT_vs_GC": ("Hippocampus", "Young (12 wk)", "ISS-terminal flight vs ground control"),
    "RRRM2_HPC_Young_LAR_FLT_vs_GC": ("Hippocampus", "Young (12 wk)", "LAR flight vs ground control"),
    "RRRM2_HPC_Old_ISST_FLT_vs_GC": ("Hippocampus", "Mature (29 wk)", "ISS-terminal flight vs ground control"),
    "RRRM2_HPC_Old_LAR_FLT_vs_GC": ("Hippocampus", "Mature (29 wk)", "LAR flight vs ground control"),
    "RRRM2_HPC_Age_GC_Old_vs_Young": ("Hippocampus", "Age contrast (GC)", "Ground control old vs young"),
    "RRRM2_HPC_Age_ISST_FLT_Old_vs_Young": ("Hippocampus", "Age contrast (FLT)", "ISS-terminal flight old vs young"),
    "RRRM2_HPC_Old_ISST_vs_LAR": ("Hippocampus", "Mature (29 wk)", "ISS-terminal vs LAR"),
    "RRRM2_HPC_Young_ISST_vs_LAR": ("Hippocampus", "Young (12 wk)", "ISS-terminal vs LAR"),
}
for c, (region, age, cond) in rrrm2.items():
    meta[c] = (f"RRRM-2 {region}", "Mouse", f"{region}, {age}", cond,
               "Modifier" if ("Age" in age or "vs LAR" in cond) else "Stressor")

assert set(meta) == set(contrasts), f"metadata missing: {set(contrasts) - set(meta)}"

# ------------------------------------------------------------- melt/join
nes_long = nes.melt(id_vars="pathway", var_name="contrast", value_name="NES")
padj_long = padj.melt(id_vars="pathway", var_name="contrast", value_name="padj")
full = nes_long.merge(padj_long, on=["pathway", "contrast"], validate="one_to_one")
full = full.merge(sets, on="pathway", validate="many_to_one")

meta_df = pd.DataFrame(
    [{"contrast": c, "dataset": v[0], "species": v[1], "tissue_model": v[2],
      "comparison": v[3], "contrast_type": v[4]} for c, v in meta.items()])
full = full.merge(meta_df, on="contrast", validate="many_to_one")

# ------------------------------------------------------------- annotate
full["significant_FDR0.05"] = full["padj"] < 0.05
full["NES_direction"] = np.sign(full["NES"]).map({1.0: "positive", -1.0: "negative"})

# ALS-spinal-cord concordance direction per category (empirical anchoring):
# positive NES in ALS spinal cord = concordant direction for that category
als_ref = full[(full["contrast_type"] == "ALS reference") &
               (full["contrast"] != "NYGC_MotorCortex")]
anchor = (als_ref.groupby("category")["NES"].mean()
          .apply(lambda v: "up" if v > 0 else ("down" if v < 0 else pd.NA))
          .rename("ALS_concordant_direction"))
full = full.merge(anchor, on="category", how="left")
full["relation_to_ALS_spinal_cord"] = full.apply(
    lambda r: ("concordant" if r["NES_direction"] == r["ALS_concordant_direction"]
               else "opposite") if pd.notna(r["ALS_concordant_direction"]) else "",
    axis=1)

# ------------------------------------------------------------- finalize
col_order = ["pathway", "tier", "category", "contrast", "dataset", "species",
             "tissue_model", "comparison", "contrast_type", "NES", "padj",
             "significant_FDR0.05", "NES_direction", "relation_to_ALS_spinal_cord"]
full = full[col_order].sort_values(["tier", "category", "pathway", "contrast"],
                                   key=lambda s: s.str.lower() if s.dtype == object else s)
full.to_csv("/mnt/results/fgsea_full_results.csv", index=False)

# ------------------------------------------------------------- report
print(f"Rows: {len(full)} (expect 155 x 37 = {155*37})")
print(f"Pathways: {full['pathway'].nunique()} | Contrasts: {full['contrast'].nunique()}")
print(f"Significant pairs (FDR<0.05): {full['significant_FDR0.05'].sum()} "
      f"({100*full['significant_FDR0.05'].mean():.1f}%)")
print()
print("Gene sets per category:")
print(full.groupby(["tier", "category"])["pathway"].nunique().to_string())
print()
print("Sanity check vs Table 2 (KEGG ALS, first 6 contrasts):")
kegg = full[full["pathway"] == "T1_KEGG_ALS"][["contrast", "NES", "padj"]]
print(kegg.head(6).to_string(index=False))

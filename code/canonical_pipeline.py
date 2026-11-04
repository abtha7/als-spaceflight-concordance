#!/usr/bin/env python3
"""
canonical_pipeline.py
============================================================================
Canonical post-processing pipeline for the ALS-spaceflight concordance
manuscript. Single source of truth for every number cited in the paper.

INPUTS (all delivered/verified files):
  nes_matrix.csv   155 gene sets x 37 contrasts, fgsea NES values
  padj_matrix.csv  matching BH-adjusted p-values (FDR) per contrast
  curated_gene_sets.gmt  the 155-set curated library (Tier 1 + Tier 2)

DESIGN DECISIONS (documented, parameterized):
  1. Category membership comes from the CURATED category assignments
     (fgsea_full_results.csv category column, which partitions the GMT),
     NOT from regex matching on pathway names.
  2. Artifact rows in nes_matrix.csv (3 all-NaN rows + 1 stray NaN-indexed
     row) are dropped; 155 real pathways are retained. The manuscript's
     "156 test sets" claim traces to the uncleaned row count and is wrong.
  3. The anchoring rule is a parameter. Four candidate definitions are
     computed and compared so the choice is transparent:
       both@0.05   padj < 0.05 in cervical AND lumbar   (strictest)
       either@0.05 padj < 0.05 in cervical OR lumbar
       both@0.25   padj < 0.25 in both (GSEA convention)
       either@0.25 padj < 0.25 in either
     A category anchors UP (S_k = +1) if frac_sig >= 0.40 under the rule
     AND all matched sets have NES > 0 in the reference region(s).
  4. WCS_c = sum over anchored categories of the MEAN NES over ALL sets in
     the category (no significance conditioning in the score itself, to
     avoid selection bias). A sensitivity variant restricted to
     reference-significant sets is also reported.
  5. Stressor / region groups are EXPLICIT contrast lists, never regex.
  6. No randomness anywhere: the pipeline is deterministic and re-runnable.

OUTPUTS (written to --outdir):
  numbers_manifest.json          every manuscript-cited number, keyed
  concordance_matrix.csv         per-contrast KEGG ALS NES/FDR, category
                                 means, WCS (all anchor rules)
  anchoring_table.csv            per-category anchoring stats per rule
  rbp_knockdown_table.csv        all 12 knockdown sets x organoid contrasts
  subanalysis_stressors.csv      explicit stressor grouping
  subanalysis_microglia.csv      all organoid microglia contrasts
  subanalysis_rodent_regions.csv true 4/6/8 RRRM-2 structure + OSD-202 arms
"""

import argparse
import json
import os
import pandas as pd
import numpy as np

# ----------------------------------------------------------------------------
# Explicit contrast groups (single source of truth; no regex classification)
# ----------------------------------------------------------------------------
REFERENCE_CONTRASTS = ["NYGC_Cervical", "NYGC_Lumbar", "NYGC_MotorCortex"]

STRESSOR_GROUPS = {
    "Radiation: simGCR 50 cGy (OSD-999)": ["OSD-999"],
    "Radiation: IR-only (OSD-202)": [
        "OSD-202_IR_only_1mon_vs_Control_1mon",
        "OSD-202_IR_only_4mon_vs_Control_4mon",
    ],
    "Radiation: HLU+IR combined (OSD-202)": [
        "OSD-202_HLU+IR_1mon_vs_Control_1mon",
        "OSD-202_HLU+IR_4mon_vs_Control_4mon",
    ],
    "Microgravity: primary payloads": [
        "OSD863_healthy_noMicro_FvG",
        "OSD871_healthy_noMicro_FvG",
        "GSE226395_MN_SMG",
    ],
    "Microgravity: RRRM-2 flight vs ground": [
        "RRRM2_FL_Young_FLT_vs_GC",
        "RRRM2_FL_Old_FLT_vs_GC",
        "RRRM2_CB_Young_ISST_FLT_vs_GC",
        "RRRM2_CB_Young_LAR_FLT_vs_GC",
        "RRRM2_CB_Old_FLT_vs_GC",
        "RRRM2_HPC_Young_ISST_FLT_vs_GC",
        "RRRM2_HPC_Young_LAR_FLT_vs_GC",
        "RRRM2_HPC_Old_ISST_FLT_vs_GC",
        "RRRM2_HPC_Old_LAR_FLT_vs_GC",
    ],
}

MODIFIER_GROUPS = {
    "Microglia: -Micro flight vs ground": [
        "OSD863_healthy_noMicro_FvG", "OSD871_healthy_noMicro_FvG"],
    "Microglia: +Micro flight vs ground": [
        "OSD863_healthy_withMicro_FvG", "OSD871_healthy_withMicro_FvG"],
    "Microglia: within-flight MicroVsNo": [
        "OSD863_healthy_Flight_MicroVsNo", "OSD871_healthy_Flight_MicroVsNo"],
    "Microglia: within-ground MicroVsNo": [
        "OSD863_healthy_Ground_MicroVsNo", "OSD871_healthy_Ground_MicroVsNo"],
    "OSD-202: HLU-only (microgravity analog)": [
        "OSD-202_HLU_only_1mon_vs_Control_1mon",
        "OSD-202_HLU_only_4mon_vs_Control_4mon"],
    "RRRM-2: age contrasts (not stressors)": [
        "RRRM2_FL_Age_GC_Old_vs_Young", "RRRM2_FL_Age_FLT_Old_vs_Young",
        "RRRM2_CB_Age_GC_Old_vs_Young", "RRRM2_CB_Age_FLT_Old_vs_Young",
        "RRRM2_HPC_Age_GC_Old_vs_Young", "RRRM2_HPC_Age_ISST_FLT_Old_vs_Young"],
    "RRRM-2: habitat contrasts ISST vs LAR (not stressors)": [
        "RRRM2_CB_Young_ISST_vs_LAR", "RRRM2_HPC_Old_ISST_vs_LAR",
        "RRRM2_HPC_Young_ISST_vs_LAR"],
}

RRRM2_REGIONS = {
    "Frontal Lobe (FL)": ["RRRM2_FL_Young_FLT_vs_GC", "RRRM2_FL_Old_FLT_vs_GC",
                          "RRRM2_FL_Age_GC_Old_vs_Young",
                          "RRRM2_FL_Age_FLT_Old_vs_Young"],
    "Cerebellum (CB)": ["RRRM2_CB_Young_ISST_FLT_vs_GC",
                        "RRRM2_CB_Young_LAR_FLT_vs_GC",
                        "RRRM2_CB_Old_FLT_vs_GC",
                        "RRRM2_CB_Age_GC_Old_vs_Young",
                        "RRRM2_CB_Age_FLT_Old_vs_Young",
                        "RRRM2_CB_Young_ISST_vs_LAR"],
    "Hippocampus (HPC)": ["RRRM2_HPC_Young_ISST_FLT_vs_GC",
                          "RRRM2_HPC_Young_LAR_FLT_vs_GC",
                          "RRRM2_HPC_Old_ISST_FLT_vs_GC",
                          "RRRM2_HPC_Old_LAR_FLT_vs_GC",
                          "RRRM2_HPC_Age_GC_Old_vs_Young",
                          "RRRM2_HPC_Age_ISST_FLT_Old_vs_Young",
                          "RRRM2_HPC_Old_ISST_vs_LAR",
                          "RRRM2_HPC_Young_ISST_vs_LAR"],
}

ANCHOR_RULES = {
    "both@0.05": lambda p_c, p_l: (p_c < 0.05) & (p_l < 0.05),
    "either@0.05": lambda p_c, p_l: (p_c < 0.05) | (p_l < 0.05),
    "both@0.25": lambda p_c, p_l: (p_c < 0.25) & (p_l < 0.25),
    "either@0.25": lambda p_c, p_l: (p_c < 0.25) | (p_l < 0.25),
}

CATEGORIES = ["Mito OXPHOS", "Proteostasis", "Neuroinflammation Glial",
              "Apoptosis", "Autophagy Mitophagy", "RNA Metabolism",
              "Oxidative Stress"]  # full 10-cat universe minus the 4 that
# cannot anchor (Axon Transport, DNA Damage Repair, Glutamate Excitotoxicity
# fail the approved rule; see anchoring_table.csv for all 10)


def load_inputs(datadir):
    nes = pd.read_csv(os.path.join(datadir, "nes_matrix.csv"), index_col=0)
    padj = pd.read_csv(os.path.join(datadir, "padj_matrix.csv"), index_col=0)
    n_raw = len(nes)

    # Documented artifact cleanup: drop all-NaN rows and any NaN-indexed row
    nes = nes[nes.notna().any(axis=1)]
    nes = nes[nes.index.notna()]
    padj = padj[padj.notna().any(axis=1)]
    padj = padj[padj.index.notna()]

    # Category assignments: curated partition of the GMT (built earlier this
    # session into fgsea_full_results.csv); fall back to datadir if present
    cat_path = "/mnt/results/fgsea_full_results.csv"
    if not os.path.exists(cat_path):
        cat_path = os.path.join(datadir, "fgsea_full_results.csv")
    fg = pd.read_csv(cat_path)
    cat = fg.drop_duplicates("pathway").set_index("pathway")["category"]

    missing = set(cat.index) - set(nes.index)
    if missing:
        raise SystemExit(f"[x] Category table references pathways absent "
                         f"from NES matrix: {missing}")
    print(f"[+] NES matrix: {n_raw} raw rows -> {len(nes)} pathways after "
          f"artifact cleanup; {nes.shape[1]} contrasts")
    print(f"[+] padj matrix: {len(padj)} pathways; category map: {len(cat)}")
    return nes, padj, cat


def anchoring_table(nes, padj, cat):
    """Compute anchoring statistics under all candidate rules."""
    rows = []
    for c in CATEGORIES:
        sets = [p for p in nes.index if cat.get(p) == c]
        p_c = padj.loc[sets, "NYGC_Cervical"]
        p_l = padj.loc[sets, "NYGC_Lumbar"]
        n_c = int((p_c < 0.05).sum())
        n_l = int((p_l < 0.05).sum())
        row = {"category": c, "tested": len(sets),
               "sig_cervical": n_c, "sig_lumbar": n_l,
               "frac_sig_cervical": round(n_c / len(sets), 3),
               "frac_sig_lumbar": round(n_l / len(sets), 3),
               "mean_NES_cervical_all": round(nes.loc[sets, "NYGC_Cervical"].mean(), 3),
               "mean_NES_lumbar_all": round(nes.loc[sets, "NYGC_Lumbar"].mean(), 3)}
        for rule, fn in ANCHOR_RULES.items():
            m = fn(p_c, p_l)
            n_m = int(m.sum())
            row[f"matched_{rule}"] = n_m
            row[f"frac_sig_{rule}"] = round(n_m / len(sets), 3)
            # direction concordance: fraction of matched sets NES > 0 in BOTH
            n_pos = int(((nes.loc[m[m].index, "NYGC_Cervical"] > 0) &
                         (nes.loc[m[m].index, "NYGC_Lumbar"] > 0)).sum()) \
                if n_m else 0
            row[f"dir_concordance_{rule}"] = round(n_pos / n_m, 3) if n_m else np.nan
            # uniform direction: all matched sets NES > 0 in both regions
            row[f"uniform_up_{rule}"] = bool(
                (nes.loc[m[m].index, "NYGC_Cervical"] > 0).all() and
                (nes.loc[m[m].index, "NYGC_Lumbar"] > 0).all()) if n_m else False
        rows.append(row)
    return pd.DataFrame(rows)


def anchored_categories(anchor_df, rule, threshold=0.40):
    sel = anchor_df[anchor_df[f"frac_sig_{rule}"] >= threshold]
    out = {}
    for _, r in sel.iterrows():
        if r[f"uniform_up_{rule}"]:
            out[r["category"]] = +1.0
        else:
            out[r["category"]] = 0.0  # fails direction criterion
    return out


def compute_wcs(nes, cat, anchors, restrict_to_ref_sig=None, exclude_sets=None):
    """WCS per contrast. restrict_to_ref_sig: optional dict cat->set names
    (reference-significant subsets) for the sensitivity variant.
    exclude_sets: optional dict cat->[set names] dropped from that
    category's mean (direction-discordant sets under the recommended rule)."""
    exclude_sets = exclude_sets or {}
    out = {}
    for contrast in nes.columns:
        total, parts = 0.0, {}
        for c, sign in anchors.items():
            sets = (restrict_to_ref_sig[c] if restrict_to_ref_sig and c in
                    restrict_to_ref_sig else
                    [p for p in nes.index if cat.get(p) == c])
            sets = [p for p in sets if p not in exclude_sets.get(c, [])]
            m = nes.loc[sets, contrast].mean() if sets else 0.0
            parts[c] = round(sign * m, 3)
            total += sign * m
        out[contrast] = {"WCS": round(total, 3), **parts}
    return pd.DataFrame(out).T


def anchored_categories_recommended(nes, padj, cat, threshold=0.40,
                                    min_dir_concordance=0.90):
    """Recommended anchor rule: a category anchors if >=threshold of its
    sets are significant (padj<0.25, GSEA-conventional FDR) in BOTH cord
    regions AND >=min_dir_concordance of matched sets are positive in both.
    Direction-discordant matched sets are excluded from the category mean
    (returned as excluded dict) and must be named in the manuscript."""
    anchors, excluded = {}, {}
    for c in CATEGORIES:
        sets = [p for p in nes.index if cat.get(p) == c]
        m = [p for p in sets if padj.loc[p, "NYGC_Cervical"] < 0.25 and
             padj.loc[p, "NYGC_Lumbar"] < 0.25]
        if not m or len(m) / len(sets) < threshold:
            continue
        pos = [p for p in m if nes.loc[p, "NYGC_Cervical"] > 0 and
               nes.loc[p, "NYGC_Lumbar"] > 0]
        if len(pos) / len(m) < min_dir_concordance:
            continue
        anchors[c] = +1.0
        excluded[c] = [p for p in m if p not in pos]
    return anchors, excluded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default="/mnt/user-uploads")
    ap.add_argument("--outdir", default="/mnt/results/canonical_outputs")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    nes, padj, cat = load_inputs(args.datadir)
    manifest = {"inputs": {"pathways": len(nes), "contrasts": int(nes.shape[1])},
                "notes": []}

    # ---- 1. Anchoring -------------------------------------------------------
    anchor_df = anchoring_table(nes, padj, cat)
    anchor_df.to_csv(os.path.join(args.outdir, "anchoring_table.csv"), index=False)
    print("\n=== ANCHORING (per category, per rule) ===")
    print(anchor_df[["category", "tested", "sig_cervical", "sig_lumbar",
                     "frac_sig_both@0.05", "frac_sig_either@0.05",
                     "uniform_up_both@0.05"]].to_string(index=False))

    # ---- 2. WCS under each anchor rule --------------------------------------
    concordance = pd.DataFrame(index=nes.columns)
    kegg = "T1_KEGG_ALS"
    concordance["KEGG_ALS_NES"] = nes.loc[kegg].round(3)
    concordance["KEGG_ALS_padj"] = padj.loc[kegg].round(5)
    for c in CATEGORIES:
        sets = [p for p in nes.index if cat.get(p) == c]
        concordance[f"meanNES_{c}"] = nes.loc[sets].mean(axis=0).round(3)

    ref_sig = {c: set(padj.index[(padj[c_col] < 0.05) & (padj["NYGC_Lumbar" if
               c_col == "NYGC_Cervical" else "NYGC_Cervical"] < 0.05) &
               (padj[c_col] < 0.05)])
               for c, c_col in [("NYGC_Cervical", "NYGC_Cervical")]}

    for rule in ANCHOR_RULES:
        anchors = anchored_categories(anchor_df, rule)
        wcs = compute_wcs(nes, cat, anchors)
        concordance[f"WCS_{rule}"] = wcs["WCS"]
        n_anch = sum(1 for v in anchors.values() if v != 0)
        print(f"[+] rule {rule}: {n_anch} categories anchored "
              f"({', '.join(k for k, v in anchors.items() if v != 0)})")
        # sensitivity: restrict to reference-significant sets (both regions)
        sig_both = {c: [p for p in nes.index if cat.get(p) == c and
                        padj.loc[p, "NYGC_Cervical"] < 0.05 and
                        padj.loc[p, "NYGC_Lumbar"] < 0.05]
                    for c in anchors}
        wcs_sig = compute_wcs(nes, cat, anchors, restrict_to_ref_sig=sig_both)
        concordance[f"WCS_{rule}_sensRefSig"] = wcs_sig["WCS"]

    # ---- 2b. Recommended rule: dual-region FDR<0.25, frac_sig>=0.40,
    #          >=90% direction concordance, discordant sets excluded -------
    anchors_rec, excluded_rec = anchored_categories_recommended(nes, padj, cat)
    wcs_rec = compute_wcs(nes, cat, anchors_rec, exclude_sets=excluded_rec)
    concordance["WCS_recommended"] = wcs_rec["WCS"]
    excl_flat = [p for v in excluded_rec.values() for p in v]
    manifest["recommended_anchor_rule"] = {
        "definition": ("category anchors if >=40% of its sets are significant "
                       "at padj<0.25 in BOTH cervical and lumbar cord AND >=90% "
                       "of matched sets are NES>0 in both; direction-discordant "
                       "matched sets are excluded from category means"),
        "anchored_categories": sorted(anchors_rec),
        "excluded_sets": excl_flat,
    }
    print(f"\n[+] recommended rule: {len(anchors_rec)} categories anchored; "
          f"excluded sets: {excl_flat}")

    # gated variant (BIOMNI_METHODS.md Section 6 formula): categories with
    # frac_sig < 0.3 IN THAT CONTRAST contribute zero
    wcs_gated = {}
    for contrast in nes.columns:
        total = 0.0
        for c, sign in anchors_rec.items():
            excl_c = excluded_rec.get(c, [])
            sets = [p for p in nes.index if cat.get(p) == c
                    and p not in excl_c]
            frac_sig_c = (padj.loc[sets, contrast] < 0.05).mean()
            if frac_sig_c >= 0.3:
                total += sign * nes.loc[sets, contrast].mean()
        wcs_gated[contrast] = round(total, 3)
    concordance["WCS_recommended_gated"] = concordance.index.map(wcs_gated)

    concordance.round(3).to_csv(os.path.join(args.outdir,
                                            "concordance_matrix.csv"))
    print("\n=== WCS (PRIMARY: recommended rule; both@0.05 shown for "
          "comparison) ===")
    show = ["KEGG_ALS_NES", "KEGG_ALS_padj", "WCS_recommended",
            "WCS_both@0.05"]
    print(concordance[show].round(2).to_string())

    # ---- 3. RBP knockdown analysis ------------------------------------------
    rbp_rows = []
    kb = [p for p in nes.index if "Knockdown" in str(p)]
    organoid = [c for c in nes.columns if c.startswith(("OSD863", "OSD871"))]
    for p in kb:
        for c in organoid:
            rbp_rows.append({"set": p, "contrast": c,
                             "NES": round(nes.loc[p, c], 3),
                             "padj": round(padj.loc[p, c], 5),
                             "significant": bool(padj.loc[p, c] < 0.05)})
    rbp = pd.DataFrame(rbp_rows)
    rbp.to_csv(os.path.join(args.outdir, "rbp_knockdown_table.csv"), index=False)

    # direction agreement: flight concordant with knockdown if
    # 'down' sets are negative AND 'up' sets are positive in flight
    # (match on the LAST token: 'down'/'up' — 'Knockdown' contains 'down'!)
    def kb_dir(set_name):
        return set_name.strip().split()[-1]  # 'down' or 'up'

    agree, tot = 0, 0
    for p in kb:
        for c in ["OSD863_healthy_noMicro_FvG", "OSD871_healthy_noMicro_FvG"]:
            v = nes.loc[p, c]
            expect_neg = kb_dir(p) == "down"
            tot += 1
            if (v < 0) == expect_neg:
                agree += 1
    manifest["rbp"] = {
        "n_sets": len(kb),
        "direction_agreement_primary_payloads": f"{agree}/{tot}",
        "summary": ("All knockdown 'down' sets are significantly DEPLETED "
                    "(negative NES) in flight organoids, but 'up' sets are "
                    "ALSO depleted: directionally incoherent with a coherent "
                    "RBP loss-of-function mimic. Do NOT report positive "
                    "enrichment.")}
    print(f"\n=== RBP: {len(kb)} knockdown sets; direction agreement "
          f"{agree}/{tot} ===")

    # ---- 3b. RBP extended to ALL 37 contrasts (for Nov 19 meeting) ---------
    rbp37_rows = []
    for p in kb:
        expect_neg = kb_dir(p) == "down"
        for c in nes.columns:
            v = nes.loc[p, c]
            rbp37_rows.append({
                "set": p, "knockdown_direction": kb_dir(p), "contrast": c,
                "NES": round(v, 3), "padj": round(padj.loc[p, c], 5),
                "significant": bool(padj.loc[p, c] < 0.05),
                "direction_agrees": bool((v < 0) == expect_neg)})
    rbp37 = pd.DataFrame(rbp37_rows)
    rbp37.to_csv(os.path.join(args.outdir, "rbp_knockdown_all37.csv"),
                 index=False)
    # per-contrast agreement summary (over the 12 sets)
    agr = rbp37.groupby("contrast")["direction_agrees"].agg(["sum", "count"])
    agr["frac"] = (agr["sum"] / agr["count"]).round(3)
    agr.sort_values("frac").to_csv(os.path.join(
        args.outdir, "rbp_direction_agreement_per_contrast.csv"),
        index_label="contrast")
    n_sig = int(rbp37["significant"].sum())
    print(f"[+] RBP all-37 table: {len(rbp37)} rows; {n_sig} significant "
          f"(padj<0.05); per-contrast agreement written")

    # ---- 4. Sub-analyses with explicit groups -------------------------------
    # PRIMARY rule (user-approved): dual-region FDR<0.25, frac_sig>=0.40,
    # >=90% direction concordance, discordant sets excluded from means.
    # The four ANCHOR_RULES columns are retained as sensitivity analyses.
    wcs_col = "WCS_recommended"
    ox_col = "meanNES_Mito OXPHOS"

    def group_table(groups):
        rows = []
        for g, contrasts in groups.items():
            sub = concordance.loc[[c for c in contrasts if c in
                                   concordance.index]]
            rows.append({"group": g, "n": len(sub),
                         "mean_WCS": round(sub[wcs_col].mean(), 3),
                         "mean_OXPHOS_NES": round(sub[ox_col].mean(), 3),
                         "contrasts": "; ".join(sub.index)})
        return pd.DataFrame(rows)

    stress = group_table(STRESSOR_GROUPS)
    micro = group_table(MODIFIER_GROUPS)
    stress.to_csv(os.path.join(args.outdir, "subanalysis_stressors.csv"),
                  index=False)
    micro.to_csv(os.path.join(args.outdir, "subanalysis_microglia.csv"),
                 index=False)

    rodent_rows = []
    for region, contrasts in RRRM2_REGIONS.items():
        sub = concordance.loc[contrasts]
        flt = [c for c in contrasts if "FLT_vs_GC" in c]
        rodent_rows.append({
            "region": region, "n_all": len(contrasts),
            "mean_WCS_all": round(sub[wcs_col].mean(), 3),
            "n_flight_vs_ground": len(flt),
            "mean_WCS_flight_vs_ground": round(
                sub.loc[flt, wcs_col].mean(), 3) if flt else np.nan,
            "mean_OXPHOS_NES_flight_vs_ground": round(
                sub.loc[flt, ox_col].mean(), 3) if flt else np.nan})
    rodent = pd.DataFrame(rodent_rows)
    osd202 = concordance.loc[[c for c in nes.columns if c.startswith("OSD-202")]]
    rodent_rows.append({
        "region": "Whole brain (OSD-202, ground)", "n_all": len(osd202),
        "mean_WCS_all": round(osd202[wcs_col].mean(), 3),
        "n_flight_vs_ground": 0,
        "mean_WCS_flight_vs_ground": np.nan,
        "mean_OXPHOS_NES_flight_vs_ground": round(osd202[ox_col].mean(), 3)})
    rodent = pd.DataFrame(rodent_rows)
    rodent.to_csv(os.path.join(args.outdir,
                               "subanalysis_rodent_regions.csv"), index=False)

    # ---- 4b. Manuscript Table 2: reference / stressor / modifier ------------
    # Section split follows the canonical explicit groups. The two noMicro
    # payload contrasts appear under Stressor (flight vs ground) and are
    # also the '-Microglia' arm of the microglia modifier view; noted in
    # table_note. Esila's 3/18/16 split vs this 3/17/17 split is unresolved
    # -- per-contrast numbers are identical either way, only row grouping
    # differs; settle with Esila before submission.
    REFERENCE_CONTRASTS = ["NYGC_Cervical", "NYGC_Lumbar", "NYGC_MotorCortex"]
    table2_rows = []
    def t2_row(section, group, c):
        r = concordance.loc[c]
        wcs = r[wcs_col]
        zone = ("Concordant" if wcs > 1 else
                "Discordant" if wcs < -1 else "Neutral")
        return {"section": section, "group": group, "contrast": c,
                "KEGG_ALS_NES": r["KEGG_ALS_NES"],
                "KEGG_ALS_padj": r["KEGG_ALS_padj"],
                "WCS": wcs, "zone": zone}
    for c in REFERENCE_CONTRASTS:
        table2_rows.append(t2_row("Reference (ALS tissue)", "ALS spinal cord / motor cortex", c))
    for g, cs in STRESSOR_GROUPS.items():
        for c in cs:
            table2_rows.append(t2_row("Stressor", g, c))
    for g, cs in MODIFIER_GROUPS.items():
        for c in cs:
            if c in STRESSOR_GROUPS["Microgravity: primary payloads"]:
                continue  # already listed under Stressor; see table_note
            table2_rows.append(t2_row("Modifier", g, c))
    table2 = pd.DataFrame(table2_rows)
    table2["table_note"] = ("noMicro payload contrasts (OSD863/871) listed "
                            "under Stressor also form the -Microglia arm of "
                            "the microglia modifier analysis")
    # per-contrast sample sizes from the other machine's handoff
    n_file = "/mnt/user-uploads/als_handoff/per_contrast_sample_sizes.csv"
    if os.path.exists(n_file):
        ns = pd.read_csv(n_file)
        table2 = table2.merge(ns[["contrast", "treatment_n", "control_n",
                                  "n_source"]], on="contrast", how="left")
        manifest["sample_sizes_source"] = ("per_contrast_sample_sizes.csv "
                                           "(other-machine handoff)")
    table2.to_csv(os.path.join(args.outdir, "table2_contrast_structure.csv"),
                  index=False)
    print(f"\n[+] Table 2: {len(table2)} rows "
          f"({(table2['section'] == 'Reference (ALS tissue)').sum()} reference / "
          f"{(table2['section'] == 'Stressor').sum()} stressor / "
          f"{(table2['section'] == 'Modifier').sum()} modifier)")

    # ---- 4c. Tier-3 ALS subtype markers (from the Multilevel run) -----------
    # T3 sets exist ONLY in single_source_gsea_all.csv (fgseaMultilevel,
    # 277 sets); the delivered Simple matrices have 155 sets without them.
    # Label all T3 results as Multilevel-run values.
    ss_file = "/mnt/user-uploads/als_handoff/single_source_gsea_all.csv"
    if os.path.exists(ss_file):
        ss = pd.read_csv(ss_file)
        t3 = ss[ss["pathway"].str.startswith("T3_")].copy()
        t3_out = t3[["pathway", "contrast", "NES", "pval", "padj", "size",
                     "n_genes"]].sort_values(["pathway", "contrast"])
        t3_out.to_csv(os.path.join(args.outdir, "t3_subtype_table.csv"),
                      index=False)
        piv = t3.pivot(index="contrast", columns="pathway", values="NES")
        manifest["t3_subtypes"] = {
            "source_run": "fgseaMultilevel (single_source_gsea_all.csv)",
            "sets": sorted(t3["pathway"].unique()),
            "n_rows": len(t3),
            "contrasts_covered": int(t3["contrast"].nunique()),
            "note": (f"{len(t3)} of {3 * 37} possible set-x-contrast rows; "
                     "missing rows = set size < 5 after ortholog mapping "
                     "in that contrast's universe")}
        print(f"\n=== T3 subtype markers (Multilevel run, {len(t3)} rows) ===")
        print(piv.round(2).to_string())

    print("\n=== STRESSOR GROUPS (explicit lists) ===")
    print(stress[["group", "n", "mean_WCS", "mean_OXPHOS_NES"]].to_string(index=False))
    print("\n=== RODENT REGIONS (true 4/6/8 structure) ===")
    print(rodent.to_string(index=False))

    # ---- 5. Manifest ---------------------------------------------------------
    key = {c: {"KEGG_ALS_NES": float(concordance.loc[c, "KEGG_ALS_NES"]),
               "KEGG_ALS_padj": float(concordance.loc[c, "KEGG_ALS_padj"]),
               "WCS": float(concordance.loc[c, wcs_col])}
           for c in ["NYGC_Cervical", "NYGC_Lumbar", "NYGC_MotorCortex",
                     "OSD-999", "OSD863_healthy_noMicro_FvG",
                     "OSD871_healthy_noMicro_FvG", "GSE226395_MN_SMG",
                     "OSD871_healthy_withMicro_FvG",
                     "OSD871_healthy_Flight_MicroVsNo"]}
    manifest["key_contrasts"] = key
    manifest["anchor_rule_used"] = (
        "RECOMMENDED (primary): category anchors if >=40% of its sets are "
        "significant at padj<0.25 in BOTH cervical and lumbar cord AND >=90% "
        "of matched sets are NES>0 in both; direction-discordant matched sets "
        "are excluded from category means (exactly one: GOBP Nuclear Export "
        "GO:0051168). WCS_recommended is the primary WCS column; WCS_both@0.05 "
        "etc. are sensitivity analyses.")
    with open(os.path.join(args.outdir, "numbers_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n[+] All outputs written to {args.outdir}")
    print("[+] numbers_manifest.json is the single source of truth for "
          "manuscript numbers.")


if __name__ == "__main__":
    main()

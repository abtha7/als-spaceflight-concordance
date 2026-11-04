# Cross-Dataset Concordance and Divergence of ALS and Spaceflight Transcriptomic Signatures — Data and Code Deposit

Data and results tables: Zenodo - DOI 10.5281/zenodo.23950969.

This archive contains all data, code, and results needed to reproduce the
analyses in the manuscript. Every number reported in the manuscript traces to
`results/canonical_outputs/numbers_manifest.json` or the tables in
`results/canonical_outputs/`.

## Layout

    code/
      canonical_pipeline.py        # primary analysis: anchoring, WCS, Table 2,
                                   #   RBP knockdown, T3 subtypes, subanalyses
      build_fgsea_table.py         # rebuilds the 155-set fgsea results table
                                   #   from the delivered NES/padj matrices
      validate_all_contrasts.py    # independent re-ranking + gseapy prerank
                                   #   validation of all 37 contrasts
      correlation_rebuild.py       # gene-level Spearman correlations with
                                   #   permutation p-values
    data/
      de_raw_files/                # per-contrast DESeq2 output (37 contrasts)
        per_contrast_DE/*.csv      #   gene, log2FC, pvalue, padj, Wald stat
      delivered_matrices/          # fgseaSimple NES and BH-adjusted p matrices
        nes_matrix.csv             #   155 gene sets x 37 contrasts
        padj_matrix.csv
      gene_sets/
        curated_gene_sets.gmt      # 155 curated gene sets (5-498 genes)
      handoff/                     # fgseaMultilevel run (seed 42) + metadata
        single_source_gsea_all.csv       # 277 sets x 37 contrasts, incl.
                                         #   T3 subtypes, T4 cell-type markers,
                                         #   KEGG-ALS-unique, 100 controls
        kegg_als_unique_and_random_fgsea.csv
        combined_category_map.csv        # category assignment for all 277 sets
        per_contrast_sample_sizes.csv    # treatment/control n per contrast
    results/
      fgsea_full_results.csv       # 155-set results table with categories
      canonical_outputs/           # ALL manuscript numbers (see below)

## Reproduction chain

1. Differential expression: DESeq2 (R/Bioconductor) per dataset — see
   manuscript Methods for designs and covariates. RRRM-2 uses NASA GeneLab
   pre-computed DESeq2 results. Per-contrast outputs are in
   `data/de_raw_files/per_contrast_DE/`.
2. Gene sets: `data/gene_sets/curated_gene_sets.gmt` (155 sets; the full
   curation also contains KEGG-ALS-unique, T3/T4 and control sets, included in
   the Multilevel run below).
3. GSEA (primary, fgseaSimple): rankings = Wald statistic (mean across
   duplicate genes), descending. fgsea, 10,000 permutations, size 5-500,
   BH correction within contrast, seed 42. Outputs:
   `data/delivered_matrices/`. Independent Python re-ranking + gseapy
   prerank (nperm = 10,000) reproduces all 37 contrasts at r >= 0.9999
   (`code/validate_all_contrasts.py`).
4. GSEA (Multilevel run, 277 sets incl. controls): fgseaMultilevel, seed 42.
   Output: `data/handoff/single_source_gsea_all.csv`.
5. Primary analysis: `python code/canonical_pipeline.py` — writes
   `results/canonical_outputs/` (anchoring, concordance/WCS, Table 2, RBP,
   T3, subanalyses, numbers manifest).
6. Gene-level correlations: `python code/correlation_rebuild.py`.

## Key manuscript numbers (all in results/canonical_outputs/)

- Primary concordance score: `concordance_matrix.csv`, column
  `WCS_recommended` (7 anchored categories; anchor rule defined in
  `numbers_manifest.json` -> `recommended_anchor_rule`).
- Negative controls: random size-matched sets 25/1,850 significant (1.4%);
  random gene removal 721/1,850 significant with 0/721 direction reversals;
  targeted OXPHOS/proteasome removal 11/37 significant with 10/11 reversed.
- Gene-level correlations: `canonical_outputs/gene_level_correlations.csv`.

## License

Data: CC-BY 4.0. Code: MIT.

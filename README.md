# Moorea_North_Shore_Sediment

Lean analysis path for Moorea north-shore sediment **Fine Fraction**
(Cooks Bay + Opunohu Bay vs out-of-bay; within-bay outlet distance).

## Paper figures (script 04)

Run from this folder:

```bash
Rscript scripts/04_paper_figures.R
```

| Figure | File | Content |
|------|------|---------|
| Fig 1 | `outputs/04_fig1_kriging_fine_fraction.png` | Kriging map of Fine Fraction |
| Fig 2 | `outputs/04_fig2_in_vs_out_bay.png` | In-bay (C+O) vs out-of-bay Welch t-test |
| Fig 3 | `outputs/04_fig3_inbay_lm.png` | In-bay-only LM: Fine Fraction ~ outlet distance |

Supporting tables/summary: `outputs/04_fig2_ttest.csv`, `outputs/04_fig3_lm_results.csv`, `outputs/04_summary.md`.

### Definitions
- **Response:** spreadsheet `Fine Fraction`
- **In-bay:** Site prefixes **C** (Cooks) and **O** (Opunohu)
- **Out-of-bay:** all other sites
- Fig 3 uses `log1p(Fine Fraction) ~ dist_outlet_m` (no grain-size terms). Shore can be re-enabled later via `fig3_include_shore` in Settings.

### Dependencies
R packages: `dplyr`, `tidyr`, `readr`, `stringr`, `ggplot2`, `sf`, `sp`, `gstat`  
Basemap cache: `outputs/02_osm_moorea.rds` (included so Fig 1 runs without re-downloading OSM).

### Note
This repository is the **paper-scoped** path only. Broader exploratory analyses live in the local `Sediment_Project_2026` working folder and are not mirrored here.

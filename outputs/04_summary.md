# Paper figures summary (script 04)

Run: 2026-10-17 12:01:06 -10
Response: Fine Fraction
In-bay prefixes: C, O (Cook's + Opunohu); all other sites = out of bay.
n analysis-ready: 50
n in-bay / out-of-bay: 23 / 27

## Fig 1
- Ordinary kriging of Fine Fraction: 04_fig1_kriging_fine_fraction.png

## Fig 2
- Welch t-test on log1p(Fine Fraction), in vs out bay: p = 0.00127
- Mean Fine Fraction in-bay = 21.9; out-of-bay = 7.62
- Plot: 04_fig2_in_vs_out_bay.png; table: 04_fig2_ttest.csv

## Fig 3 (in-bay only; no grain size)
- Main model: log1p(Fine Fraction) ~ outlet (fig3_include_shore = FALSE)
- log1p(Fine Fraction) ~ outlet: adj R² = 0.401, outlet p = 0.000707, AIC = 64.18
- Outlet + shore model not run (set fig3_include_shore <- TRUE to enable).
- Plot: 04_fig3_inbay_lm.png
- Coefficients: 04_fig3_lm_results.csv

## Note
- Script 03 exploratory analyses are retained as archive context; script 04 is the lean paper figure path.

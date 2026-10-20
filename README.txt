HAB mixed-model analysis (Supplemental Fig. 8 and Supplemental Table 2)
=======================================================================

These files reproduce the negative-binomial mixed models reported in Supplemental Fig. 8 and
Supplemental Table 2. They start from HAB counts per larva and time point. They do not include
the HAB-detection code that converts EthoVision activity into HABs.

data/
  fig7_long.csv        No-treatment experiment (Fig. 7)
  fig8_pre_long.csv    Indomethacin experiment, before treatment (as plotted in Fig. 8A; includes the Fig. 7 larvae)
  fig8_post_long.csv   Indomethacin experiment, after treatment
  fig9_pre_long.csv    Retigabine experiment, before treatment
  fig9_post_long.csv   Retigabine experiment, after treatment
  fig10_pre_long.csv   MCC950 experiment, before treatment
  fig10_post_long.csv  MCC950 experiment, after treatment

  Columns: experiment, segment, larva (ID unique within a file), group, tbi (1 = injured),
  drug (1 = drug, 0 = vehicle, empty before treatment), plate, session (plate x time point),
  tp (time point, h post-injury), tp_order, count (HABs on the recording's own bin scale),
  bins (number of one-minute bins in the recording: 120, 180 for 3-h recordings, 103 for the
  shortened Fig. 8 recording).

supp_fig8_models.R
  R version using glmmTMB (Laplace approximation). Run from this folder:  Rscript supp_fig8_models.R
  Writes results_R.csv and results_R_tests.csv and prints the original values for comparison.

python/fit_from_csv.py, python/glmm.py
  The original analysis (maximum likelihood with adaptive Gauss-Hermite quadrature, NumPy/SciPy).
  Run from this folder:  python python/fit_from_csv.py data results_python.csv

results_python.csv, results_python_tests.csv
  Output of the original analysis: rate ratio, 95% CI, unadjusted and Holm-adjusted p for each
  time point and comparison, and the likelihood-ratio tests.

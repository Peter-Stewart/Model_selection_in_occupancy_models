## Overview ##
This folder contains subfolders with the results dataframes for the main text (simulations with K=2) and supplementary material (simulations with K=40). The list of dataframes in each subfolder, and the variable descriptions for each dataframe, are:

## File list ##
* **resultsdf_m1_occ.Rdata**  - model 1 results for scenario with M-bias in occupancy process only
* **resultsdf_m2_occ.Rdata**  - model 2 results for scenario with M-bias in occupancy process only
* **resultsdf_m1_det.Rdata**  - model 1 results for scenario with M-bias in detection process only
* **resultsdf_m2_det.Rdata** - model 2 results for scenario with M-bias in detection process only
* **resultsdf_m1_occ_det.Rdata** - model 1 results for scenario with M-bias in occupancy and detection processes
* **resultsdf_m2_occ_det.Rdata** - model 2 results for scenario with M-bias in occupancy and detection processes
* **resultsdf_m3_occ_det.Rdata** - model 3 results for scenario with M-bias in occupancy and detection processes
* **resultsdf_m4_occ_det.Rdata** - model 4 results for scenario with M-bias in occupancy and detection processes
* **predictions_all_occ.Rdata** - predictive accuracy metrics for all models in scenario with M-bias in occupancy process only
* **retrodictions_all_occ.Rdata** - retrodictive accuracy metrics for all models in scenario with M-bias in occupancy process only
* **predictions_all_det.Rdata** - predictive accuracy metrics for all models in scenario with M-bias in detection process only
* **retrodictions_all_det.Rdata** - retrodictive accuracy metrics for all models in scenario with M-bias in detection process only
* **predictions_all_occ_det.Rdata** - predictive accuracy metrics for all models in scenario with M-bias in occupancy and detection processes
* **retrodictions_all_occ_det.Rdata** - retrodictive accuracy metrics for all models in scenario with M-bias in occupancy and detection processes
* **Akaike_all_occ.Rdata** - AIC and BIC values and weights for all models in scenario with M-bias in occupancy process only
* **Akaike_all_det.Rdata** - AIC and BIC values and weights for all models in scenario with M-bias in detection process only
* **Akaike_all_occ_det.Rdata** - AIC and BIC values and weights for all models in scenario with M-bias in occupancy and detection processes


## Variable descriptions ##
**resultsdf_m1_occ**:
* **column_label** - iteration of simulation
* **Intercept.est** - estimate of occupancy intercept
* **X.est** - estimate of effect of covariate X on occupancy (psi)
* **Intercept.se** - standard error of estimate of occupancy intercept
* **X.se** - standard error of estimate of effect of covariate X on occupancy (psi)
* **Intercept.zscore** - z-score for estimate of occupancy intercept
* **X.zscore** - z-score for estimate of effect of covariate X on occupancy (psi)
* **Intercept.pval** - p-value for estimate of occupancy intercept
* **X.pval** - pvalue for estimate of effect of covariate X on occupancy (psi)
* **det.est** - estimate of detection intercept
* **det.se** - standard error of estimate of detection intercept
* **det.z** - z-score for estimate of detection intercept
* **det.pval** - p-value for estimate of detection intercept
* **B1** - true effect of X on occupancy (psi)
* **B2** - true effect of A on X
* **B3** - true effect of A on D
* **B4** - true effect of C on D
* **B5** - true effect of C on occupancy (psi)
* **alpha** - true occupancy intercept
* **det_prob** - true detection probability
* **logit_det_prob** - logit of true detection probability
* **X.ci.up** - 95% confidence interval of effect of covariate X on occupancy (psi), upper value
* **X.ci.low** - 95% confidence interval of effect of covariate X on occupancy (psi), lower value
* **X.in95ci** - 1 = true effect of X (B1) is contained in 95% confidence interval, 0 = true effect of X (B1) is not contained in 95% confidence interval
* **X.bias** - bias of estimate of effect of covariate X on occupancy (psi): bias = X.est - B1

**resultsdf_m1_occ**:
* **column_label** - iteration of simulation
* **Intercept.est** - estimate of occupancy intercept
* **X.est** - estimate of effect of covariate X on occupancy (psi)
* **D.est** - estimate of effect of covariate D on occupancy (psi)
* **Intercept.se** - standard error of estimate of occupancy intercept
* **X.se** - standard error of estimate of effect of covariate X on occupancy (psi)
* **D.se** - standard error of estimate of effect of covariate D on occupancy (psi)
* **Intercept.zscore** - z-score for estimate of occupancy intercept
* **X.zscore** - z-score for estimate of effect of covariate X on occupancy (psi)
* **D.zscore** - z-score for estimate of effect of covariate D on occupancy (psi)
* **Intercept.pval** - p-value for estimate of occupancy intercept
* **X.pval** - pvalue for estimate of effect of covariate X on occupancy (psi)
* **D.pval** - pvalue for estimate of effect of covariate D on occupancy (psi)
* **det.est** - estimate of detection intercept
* **det.se** - standard error of estimate of detection intercept
* **det.z** - z-score for estimate of detection intercept
* **det.pval** - p-value for estimate of detection intercept
* **B1** - true effect of X on occupancy (psi)
* **B2** - true effect of A on X
* **B3** - true effect of A on D
* **B4** - true effect of C on D
* **B5** - true effect of C on occupancy (psi)
* **alpha** - true occupancy intercept
* **det_prob** - true detection probability
* **logit_det_prob** - logit of true detection probability
* **X.ci.up** - 95% confidence interval of effect of covariate X on occupancy (psi), upper value
* **X.ci.low** - 95% confidence interval of effect of covariate X on occupancy (psi), lower value
* **X.in95ci** - 1 = true effect of X (B1) is contained in 95% confidence interval, 0 = true effect of X (B1) is not contained in 95% confidence interval
* **X.bias** - bias of estimate of effect of covariate X on occupancy (psi): bias = X.est - B1

**resultsdf_m1_det**:
* **column_label** - iteration of simulation
* **psi.intercept.est** - estimate of occupancy intercept
* **X.est** - estimate of effect of covariate X on occupancy (psi)
* **psi.intercept.se** - standard error of estimate of occupancy intercept
* **X.se** - standard error of estimate of effect of covariate X on occupancy (psi)
* **psi.intercept.zscore** - z-score for estimate of occupancy intercept
* **X.zscore** - z-score for estimate of effect of covariate X on occupancy (psi)
* **psi.intercept.pval** - p-value for estimate of occupancy intercept
* **X.pval** - pvalue for estimate of effect of covariate X on occupancy (psi)
* **det.intercept.est** - estimate of detection intercept
* **det.U.est** - estimate of effect of covariate U on detection
* **det.intercept.se** - standard error of estimate of detection intercept
* **det.U.se** - standard error of estimate of effect of covariate U on detection
* **det.intercept.z** - z-score for estimate of detection intercept
* **det.U.z** - z-score for estimate of effect of covariate U on detection
* **det.intercept.pval** - p-value for estimate of detection intercept
* **det.U.pval** - p-value for estimate of effect of covariate U on detection
* **B1** - true effect of X on occupancy (psi)
* **B6** - true effect of U on detection
* **B7** - true effect of Q on U
* **B8** - true effect of Q on R
* **B9** - true effect of V on R
* **B10** - true effect of V on detection
* **alpha** - true occupancy intercept
* **alpha2** - true detection intercept
* **X.ci.up** - 95% confidence interval of effect of covariate X on occupancy (psi), upper value
* **X.ci.low** - 95% confidence interval of effect of covariate X on occupancy (psi), lower value
* **X.in95ci** - 1 = true effect of X (B1) is contained in 95% confidence interval, 0 = true effect of X (B1) is not contained in 95% confidence interval
* **X.bias** - bias of estimate of effect of covariate X on occupancy (psi): bias = X.est - B1

**resultsdf_m2_det**:
* **column_label** - iteration of simulation
* **psi.intercept.est** - estimate of occupancy intercept
* **X.est** - estimate of effect of covariate X on occupancy (psi)
* **psi.intercept.se** - standard error of estimate of occupancy intercept
* **X.se** - standard error of estimate of effect of covariate X on occupancy (psi)
* **psi.intercept.zscore** - z-score for estimate of occupancy intercept
* **X.zscore** - z-score for estimate of effect of covariate X on occupancy (psi)
* **psi.intercept.pval** - p-value for estimate of occupancy intercept
* **X.pval** - pvalue for estimate of effect of covariate X on occupancy (psi)
* **det.intercept.est** - estimate of detection intercept
* **det.U.est** - estimate of effect of covariate U on detection
* **det.R.est** - estimate of effect of covariate R on detection
* **det.intercept.se** - standard error of estimate of detection intercept
* **det.U.se** - standard error of estimate of effect of covariate U on detection
* **det.R.se** - standard error of estimate of effect of covariate R on detection
* **det.intercept.z** - z-score for estimate of detection intercept
* **det.U.z** - z-score for estimate of effect of covariate U on detection
* **det.R.z** - z-score for estimate of effect of covariate R on detection
* **det.intercept.pval** - p-value for estimate of detection intercept
* **det.U.pval** - p-value for estimate of effect of covariate U on detection
* **det.R.pval** - p-value for estimate of effect of covariate U on detection
* **B1** - true effect of X on occupancy (psi)
* **B6** - true effect of U on detection
* **B7** - true effect of Q on U
* **B8** - true effect of Q on R
* **B9** - true effect of V on R
* **B10** - true effect of V on detection
* **alpha** - true occupancy intercept
* **alpha2** - true detection intercept
* **X.ci.up** - 95% confidence interval of effect of covariate X on occupancy (psi), upper value
* **X.ci.low** - 95% confidence interval of effect of covariate X on occupancy (psi), lower value
* **X.in95ci** - 1 = true effect of X (B1) is contained in 95% confidence interval, 0 = true effect of X (B1) is not contained in 95% confidence interval
* **X.bias** - bias of estimate of effect of covariate X on occupancy (psi): bias = X.est - B1

**resultsdf_m1_occ_det**:
* **column_label** - iteration of simulation
* **psi.intercept.est** - estimate of occupancy intercept
* **X.est** - estimate of effect of covariate X on occupancy (psi)
* **psi.intercept.se** - standard error of estimate of occupancy intercept
* **X.se** - standard error of estimate of effect of covariate X on occupancy (psi)
* **Intercept.zscore** - z-score for estimate of occupancy intercept
* **X.zscore** - z-score for estimate of effect of covariate X on occupancy (psi)
* **Intercept.pval** - p-value for estimate of occupancy intercept
* **X.pval** - pvalue for estimate of effect of covariate X on occupancy (psi)
* **det.intercept.est** - estimate of detection intercept
* **det.U.est** - estimate of effect of covariate U on detection
* **det.intercept.se** - standard error of estimate of detection intercept
* **det.U.se** - standard error of estimate of effect of covariate U on detection
* **det.intercept.z** - z-score for estimate of detection intercept
* **det.U.z** - z-score for estimate of effect of covariate U on detection
* **det.intercept.pval** - p-value for estimate of detection intercept
* **det.U.pval** - p-value for estimate of effect of covariate U on detection
* **B1** - true effect of X on occupancy (psi)
* **B2** - true effect of A on X
* **B3** - true effect of A on D
* **B4** - true effect of C on D
* **B5** - true effect of C on occupancy (psi)
* **B6** - true effect of U on detection
* **B7** - true effect of Q on U
* **B8** - true effect of Q on R
* **B9** - true effect of V on R
* **B10** - true effect of V on detection
* **alpha** - true occupancy intercept
* **alpha2** - true detection intercept
* **X.ci.up** - 95% confidence interval of effect of covariate X on occupancy (psi), upper value
* **X.ci.low** - 95% confidence interval of effect of covariate X on occupancy (psi), lower value
* **X.in95ci** - 1 = true effect of X (B1) is contained in 95% confidence interval, 0 = true effect of X (B1) is not contained in 95% confidence interval
* **X.bias** - bias of estimate of effect of covariate X on occupancy (psi): bias = X.est - B1

**resultsdf_m2_occ_det**:
* **column_label** - iteration of simulation
* **psi.intercept.est** - estimate of occupancy intercept
* **X.est** - estimate of effect of covariate X on occupancy (psi)
* **D.est** - estimate of effect of covariate D on occupancy (psi)
* **psi.intercept.se** - standard error of estimate of occupancy intercept
* **X.se** - standard error of estimate of effect of covariate X on occupancy (psi)
* **D.se** - standard error of estimate of effect of covariate D on occupancy (psi)
* **Intercept.zscore** - z-score for estimate of occupancy intercept
* **X.zscore** - z-score for estimate of effect of covariate X on occupancy (psi)
* **D.zscore** - z-score for estimate of effect of covariate D on occupancy (psi)
* **Intercept.pval** - p-value for estimate of occupancy intercept
* **X.pval** - pvalue for estimate of effect of covariate X on occupancy (psi)
* **D.pval** - pvalue for estimate of effect of covariate D on occupancy (psi)
* **det.intercept.est** - estimate of detection intercept
* **det.U.est** - estimate of effect of covariate U on detection
* **det.intercept.se** - standard error of estimate of detection intercept
* **det.U.se** - standard error of estimate of effect of covariate U on detection
* **det.intercept.z** - z-score for estimate of detection intercept
* **det.U.z** - z-score for estimate of effect of covariate U on detection
* **det.intercept.pval** - p-value for estimate of detection intercept
* **det.U.pval** - p-value for estimate of effect of covariate U on detection
* **B1** - true effect of X on occupancy (psi)
* **B2** - true effect of A on X
* **B3** - true effect of A on D
* **B4** - true effect of C on D
* **B5** - true effect of C on occupancy (psi)
* **B6** - true effect of U on detection
* **B7** - true effect of Q on U
* **B8** - true effect of Q on R
* **B9** - true effect of V on R
* **B10** - true effect of V on detection
* **alpha** - true occupancy intercept
* **alpha2** - true detection intercept
* **X.ci.up** - 95% confidence interval of effect of covariate X on occupancy (psi), upper value
* **X.ci.low** - 95% confidence interval of effect of covariate X on occupancy (psi), lower value
* **X.in95ci** - 1 = true effect of X (B1) is contained in 95% confidence interval, 0 = true effect of X (B1) is not contained in 95% confidence interval
* **X.bias** - bias of estimate of effect of covariate X on occupancy (psi): bias = X.est - B1

**resultsdf_m3_occ_det**:
* **column_label** - iteration of simulation
* **psi.intercept.est** - estimate of occupancy intercept
* **X.est** - estimate of effect of covariate X on occupancy (psi)
* **psi.intercept.se** - standard error of estimate of occupancy intercept
* **X.se** - standard error of estimate of effect of covariate X on occupancy (psi)
* **Intercept.zscore** - z-score for estimate of occupancy intercept
* **X.zscore** - z-score for estimate of effect of covariate X on occupancy (psi)
* **Intercept.pval** - p-value for estimate of occupancy intercept
* **X.pval** - pvalue for estimate of effect of covariate X on occupancy (psi)
* **det.intercept.est** - estimate of detection intercept
* **det.U.est** - estimate of effect of covariate U on detection
* **det.R.est** - estimate of effect of covariate R on detection
* **det.intercept.se** - standard error of estimate of detection intercept
* **det.U.se** - standard error of estimate of effect of covariate U on detection
* **det.R.se** - standard error of estimate of effect of covariate R on detection
* **det.intercept.z** - z-score for estimate of detection intercept
* **det.U.z** - z-score for estimate of effect of covariate U on detection
* **det.R.z** - z-score for estimate of effect of covariate R on detection
* **det.intercept.pval** - p-value for estimate of detection intercept
* **det.U.pval** - p-value for estimate of effect of covariate U on detection
* **det.R.pval** - p-value for estimate of effect of covariate U on detection
* **B1** - true effect of X on occupancy (psi)
* **B2** - true effect of A on X
* **B3** - true effect of A on D
* **B4** - true effect of C on D
* **B5** - true effect of C on occupancy (psi)
* **B6** - true effect of U on detection
* **B7** - true effect of Q on U
* **B8** - true effect of Q on R
* **B9** - true effect of V on R
* **B10** - true effect of V on detection
* **alpha** - true occupancy intercept
* **alpha2** - true detection intercept
* **X.ci.up** - 95% confidence interval of effect of covariate X on occupancy (psi), upper value
* **X.ci.low** - 95% confidence interval of effect of covariate X on occupancy (psi), lower value
* **X.in95ci** - 1 = true effect of X (B1) is contained in 95% confidence interval, 0 = true effect of X (B1) is not contained in 95% confidence interval
* **X.bias** - bias of estimate of effect of covariate X on occupancy (psi): bias = X.est - B1

**resultsdf_m4_occ_det**:
* **column_label** - iteration of simulation
* **psi.intercept.est** - estimate of occupancy intercept
* **X.est** - estimate of effect of covariate X on occupancy (psi)
* **D.est** - estimate of effect of covariate D on occupancy (psi)
* **psi.intercept.se** - standard error of estimate of occupancy intercept
* **X.se** - standard error of estimate of effect of covariate X on occupancy (psi)
* **D.se** - standard error of estimate of effect of covariate D on occupancy (psi)
* **Intercept.zscore** - z-score for estimate of occupancy intercept
* **X.zscore** - z-score for estimate of effect of covariate X on occupancy (psi)
* **D.zscore** - z-score for estimate of effect of covariate D on occupancy (psi)
* **Intercept.pval** - p-value for estimate of occupancy intercept
* **X.pval** - pvalue for estimate of effect of covariate X on occupancy (psi)
* **D.pval** - pvalue for estimate of effect of covariate D on occupancy (psi)
* **det.intercept.est** - estimate of detection intercept
* **det.U.est** - estimate of effect of covariate U on detection
* **det.R.est** - estimate of effect of covariate R on detection
* **det.intercept.se** - standard error of estimate of detection intercept
* **det.U.se** - standard error of estimate of effect of covariate U on detection
* **det.R.se** - standard error of estimate of effect of covariate R on detection
* **det.intercept.z** - z-score for estimate of detection intercept
* **det.U.z** - z-score for estimate of effect of covariate U on detection
* **det.R.z** - z-score for estimate of effect of covariate R on detection
* **det.intercept.pval** - p-value for estimate of detection intercept
* **det.U.pval** - p-value for estimate of effect of covariate U on detection
* **det.R.pval** - p-value for estimate of effect of covariate U on detection
* **B1** - true effect of X on occupancy (psi)
* **B2** - true effect of A on X
* **B3** - true effect of A on D
* **B4** - true effect of C on D
* **B5** - true effect of C on occupancy (psi)
* **B6** - true effect of U on detection
* **B7** - true effect of Q on U
* **B8** - true effect of Q on R
* **B9** - true effect of V on R
* **B10** - true effect of V on detection
* **alpha** - true occupancy intercept
* **alpha2** - true detection intercept
* **X.ci.up** - 95% confidence interval of effect of covariate X on occupancy (psi), upper value
* **X.ci.low** - 95% confidence interval of effect of covariate X on occupancy (psi), lower value
* **X.in95ci** - 1 = true effect of X (B1) is contained in 95% confidence interval, 0 = true effect of X (B1) is not contained in 95% confidence interval
* **X.bias** - bias of estimate of effect of covariate X on occupancy (psi): bias = X.est - B1

**predictions_all_occ**, **retrodictions_all_occ**, **predictions_all_det**, and **retrodictions_all_det**:
* **column_label** - iteration of simulation
* **m1_mean_error** - mean error (mean difference between predicted/retrodicted and true psi across sites) for model 1
* **m1_mean_abs_error** - mean absolute error (mean absolute difference between predicted/retrodicted and true psi across sites) for model 1
* **m1_prop_in95ci** - proportion of sites where true occupancy probability is contained within 95% confidence interval of estimate for model 1
* **m2_mean_error** - mean error (mean difference between predicted/retrodicted and true psi across sites) for model 2
* **m2_mean_abs_error** - mean absolute error (mean absolute difference between predicted/retrodicted and true psi across sites) for model 2
* **m2_prop_in95ci** - proportion of sites where true occupancy probability is contained within 95% confidence interval of estimate for model 2

**predictions_all_occ_det** and **retrodictions_all_occ_det**:
* **column_label** - iteration of simulation
* **m1_mean_error** - mean error (mean difference between predicted/retrodicted and true psi across sites) for model 1
* **m1_mean_abs_error** - mean absolute error (mean absolute difference between predicted/retrodicted and true psi across sites) for model 1
* **m1_prop_in95ci** - proportion of sites where true occupancy probability is contained within 95% confidence interval of estimate for model 1
* **m2_mean_error** - mean error (mean difference between predicted/retrodicted and true psi across sites) for model 2
* **m2_mean_abs_error** - mean absolute error (mean absolute difference between predicted/retrodicted and true psi across sites) for model 2
* **m2_prop_in95ci** - proportion of sites where true occupancy probability is contained within 95% confidence interval of estimate for model 2
* * **m3_mean_error** - mean error (mean difference between predicted/retrodicted and true psi across sites) for model 3
* **m3_mean_abs_error** - mean absolute error (mean absolute difference between predicted/retrodicted and true psi across sites) for model 3
* **m3_prop_in95ci** - proportion of sites where true occupancy probability is contained within 95% confidence interval of estimate for model 3
* * **m4_mean_error** - mean error (mean difference between predicted/retrodicted and true psi across sites) for model 4
* **m4_mean_abs_error** - mean absolute error (mean absolute difference between predicted/retrodicted and true psi across sites) for model 4
* **m4_prop_in95ci** - proportion of sites where true occupancy probability is contained within 95% confidence interval of estimate for model 4

**Akaike_all_occ** and **Akaike_all_det**:
* **column_label** - iteration of simulation
* **m1_AIC** - Akaike's Information Criterion for model 1
* **m2_AIC** - Akaike's Information Criterion for model 2
* **m1_AIC_delta** - delta AIC (difference between model's AIC and best AIC in candidate set) for model 1
* **m2_AIC_delta** - delta AIC (difference between model's AIC and best AIC in candidate set) for model 2
* **m1_weight** - Akaike weight for model 1
* **m2_weight** - Akaike weight for model 2
* **m1_BIC** - Schwarz criterion (aka Bayesian Information Criterion) for model 1
* **m2_BIC** - Schwarz criterion (aka Bayesian Information Criterion) for model 2
* **m1_BIC_delta** - delta BIC (difference between model's BIC and best BIC in candidate set) for model 1
* **m2_BIC_delta** - delta BIC (difference between model's BIC and best BIC in candidate set) for model 2
* **m1_BIC_weight** - BIC weight for model 1
* **m2_BIC_weight** - BIC weight for model 2

**Akaike_all_occ_det**:
* **column_label** - iteration of simulation
* **m1_AIC** - Akaike's Information Criterion for model 1
* **m2_AIC** - Akaike's Information Criterion for model 2
* **m3_AIC** - Akaike's Information Criterion for model 3
* **m4_AIC** - Akaike's Information Criterion for model 4
* **m1_AIC_delta** - delta AIC (difference between model's AIC and best AIC in candidate set) for model 1
* **m2_AIC_delta** - delta AIC (difference between model's AIC and best AIC in candidate set) for model 2
* **m3_AIC_delta** - delta AIC (difference between model's AIC and best AIC in candidate set) for model 3
* **m4_AIC_delta** - delta AIC (difference between model's AIC and best AIC in candidate set) for model 4
* **m1_weight** - Akaike weight for model 1
* **m2_weight** - Akaike weight for model 2
* **m3_weight** - Akaike weight for model 3
* **m4_weight** - Akaike weight for model 4
* **m1_BIC** - Schwarz criterion (aka Bayesian Information Criterion) for model 1
* **m2_BIC** - Schwarz criterion (aka Bayesian Information Criterion) for model 2
* **m3_BIC** - Schwarz criterion (aka Bayesian Information Criterion) for model 3
* **m4_BIC** - Schwarz criterion (aka Bayesian Information Criterion) for model 4
* **m1_BIC_delta** - delta BIC (difference between model's BIC and best BIC in candidate set) for model 1
* **m2_BIC_delta** - delta BIC (difference between model's BIC and best BIC in candidate set) for model 2
* **m3_BIC_delta** - delta BIC (difference between model's BIC and best BIC in candidate set) for model 3
* **m4_BIC_delta** - delta BIC (difference between model's BIC and best BIC in candidate set) for model 4
* **m1_BIC_weight** - BIC weight for model 1
* **m2_BIC_weight** - BIC weight for model 2
* **m3_BIC_weight** - BIC weight for model 3
* **m4_BIC_weight** - BIC weight for model 4




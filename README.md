# Global SDM Model Generalization

Study II: independent multi-region extension of [Study I](https://github.com/ys30/Ecological-SDM-XGBoost/tree/a5231e0766239b539edb638b8150cf73d1498f9e). Study I is read-only; this snapshot is the observed repository head, not yet verified as the exact submitted version.

## Questions
1. Do nonlinear specifications rival algorithm family across regions?
2. How do spatial validation, prevalence, sample size and environmental novelty modify paired performance differences?
3. How much pooled AUC reflects between-species rather than within-species comparisons?

## Data
Pinned disdat 1.1-0, commit 71b878e85ea02764386cfcdbff97816a8f18a8df; Elith et al. (2020), doi:10.17161/bi.v15i2.13384. Six regions (AWT, CAN, NSW, NZ, SA, SWI), 226 candidate species. Audit all, including original taxa; no outcome-based exclusions. PA and PO/background experiments remain separate. Do not redistribute raw data without checking terms.

## Gates
M00 source audit; M01-02 data inventory; M03 predictor and spatial audit; M04 exact original registry; M05 NSW reproduction; M06 pilot; M07 full benchmark; M08 pooled AUC decomposition; M09 hierarchical analysis; M10 independent external dataset; M11 manuscript.

Primary endpoint: paired within-species ROC AUC on identical held-out observations, macro average by species and region. Secondary: AP, Brier and calibration for PA probability models; PO scores are not occurrence probabilities. Nested spatial tuning; transformations fitted only on training data. Group coincident sites; retain failed folds and exclusions. Leave-one-region-out evaluates selection rules, not transfer between unrelated taxa. Region-level uncertainty is limited by six regions. No results yet.

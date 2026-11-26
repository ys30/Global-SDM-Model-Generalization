# Study I replication source audit

Read-only audit of ys30/Ecological-SDM-XGBoost, public main at a5231e0766239b539edb638b8150cf73d1498f9e. This is not a certified submission snapshot.
Public README reports v6.26m, 1167 plots, 8 taxa, 9336 plot-taxon records, 747 detections; 14 configurations (12 common across 23 schemes). GLM quad-4 wins 16/23; quadratic regression wins 19/23. The manuscript itself is private. README identifies research commit 6d54748421421048f8d2be13d862861ea6c120ce and archive DOI 10.5281/zenodo.24451053; archive contents and license not verified here.

## Source inventory
- README.md: manuscript summary and historical pipeline distinctions.
- results/v625/table_tuned_summary.csv: author-provided aggregate results, not independent CI replay.
- results/v625/README.md: explicitly states full package not publicly uploaded.
- Public scripts 50–63 have different names from manuscript-package scripts 50–63; do not equate the numbering.
- Original historical PA MaxEnt uses surveyed nondetections as background approximation; Study II PO MaxEnt uses independent presence/background training. These are different tasks.
- Study I dataset must be matched by actual site IDs, species IDs, predictor definitions, and outcomes, not simply region NSW. No exact match established.

## Replication gates
1. Recover exact submitted manuscript/SI and author archive or licensed equivalent.
2. Inventory all 14 named configurations, predictor sets N4/RT/categorical, penalty strengths, software versions, pooling structures and hyperparameter search.
3. Verify exact 1167 sites and 8 taxa, response matrix and predictor values.
4. Recover all 23 fold assignments, buffers, eligibility and weighting definitions.
5. Replay and compare reported AUC and winner counts; document tolerances before running.
6. Apply exact configurations to eligible Study II groups separately from expanded benchmark.
Until gates pass, call Study II a conceptual extension, not exact replication.

## Verified public summary numbers (author provided, not replayed)
| Design | Quad GLM | tuned per-taxon XGB | tuned pooled XGB | MaxEnt |
|---|---:|---:|---:|---:|
| grid4_s2026_b0 | 0.8188012010601716 | 0.7983896915279673 | 0.798065909190339 | 0.7866038871121852 |
| km15_b0 | 0.8201231202045075 | 0.7936687436791966 | 0.7895594757761486 | 0.810168060551469 |
| grid4_s2026_b10 | 0.8146043018424263 | 0.7801463260094813 | 0.7887777100551758 | 0.7812856688543242 |

Study II's current exploratory findings cannot be interpreted as a direct contradiction of Study I until these gates pass.

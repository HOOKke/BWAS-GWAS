## HCP validations

- `01_compute_hcp_idp_behavior_asso.R`

  Computes HCP structural IDP-behavior association vectors. For each selected
  behavioral phenotype, the script tests associations between structural IDPs
  and phenotype values while adjusting for age, sex/gender, education,
  handedness, and race.

- `02_hcp_ukb_bwas_pattern_similarity.py`

  Aligns the HCP IDP-behavior association vectors with UK Biobank CCA brain
  loading vectors for the same selected structural IDPs.


- `03_prs_hcp_validation.R`

  Runs the HCP PRS validation workflow in two stages. First, for each PRS mode,
  the script tests associations between the mode-specific PRS scores and the
  four targeted HCP phenotypes. PRS scores and phenotype values are residualized
  for age, sex/gender, handedness, race, and genetic PCs 1-10, and adjusted
  PRS-phenotype correlations are computed across PRS P-value thresholds. Second,
  after the batch association analysis, the strongest PRS P-value threshold is
  used for four-quantile visualization.

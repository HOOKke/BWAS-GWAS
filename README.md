# BWAS-GWAS
## Introduction
This repository contains code supporting the manuscript **“Genetic mechanisms underlying the associations between brain structure and mental health.”**
## Data
Our analysis pipeline relies on individual-level **UK Biobank** data (phenotypic, neuroimaging, and genotype), which are provided under controlled access agreements and cannot be redistributed or uploaded to third-party platforms. Access to these data is available to bona fide researchers through application to the UK Biobank website (https://www.ukbiobank.ac.uk). Additional information about registration for access to the data is available at http://www.ukbiobank.ac.uk/register-apply/.
## Analysis
Scripts are organized into three subfolders within the analysis directory, corresponding to the order of analyses reported in the manuscript.
### 01_CCA:
Canonical Correlation Analysis (CCA) was performed to identify multivariate associations between 7 mental health measures and 277 brain structural IDPs in UK Biobank individuals. Dimensionality reduction via principal component analysis (PCA) was applied to brain structural IDPs before running CCA, yielding top 70 principal components that explained 72.6% of the variance.
### 02_GWAS:
For each significant CCA mode, we carried out subsequent genetic analyses of the brain scores (U). Prior to genetic analysis, we estimated the impact of brain scores on survival using Cox proportional hazard regression models. Then, genome-wide association studies (GWAS) was performed on brain scores for each significant CCA mode using PLINK 1.9 with additional adjusted for the top 40 genetic principal components. 
### 03_Post-GWAS:
1. **SNP_heritability**: Using the summary statistics, we estimated the SNP-based heritability of the brain scores using LDSC (https://github.com/bulik/ldsc). LD scores were sourced from the European 1000 Genomes sample.
2. **Genetic_correlations**: We used LDSC to estimate the genetic correlations between the brain scores and 15 other complex traits, including reaction time, educational attainment (overall, cognitive component, and non-cognitive component), cognitive performance, intelligence, extremely high intelligence, alcohol dependence, anxiety disorder, cannabis use disorder, ADHD, PTSD, SCZ, ASD, and MDD (Supplementary Table 4).
3. **PRS_analysis**: We selected independent white British subjects for which did not overlap with the GWAS sample and which both genetic and behavioral data were available. The polygenic risk scores (PRSs) for the brain scores of these individuals were calculated using PLINK based on our GWAS results. Then, the associations between PRS and mental health traits were estimated using Pearson’s correlation analysis.

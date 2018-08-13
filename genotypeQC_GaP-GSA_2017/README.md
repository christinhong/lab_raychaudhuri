# genotypeQC_GaP-GSA_2017
Project: Immunomics

*This is a private repo for Christin to develop and share code in the Raychaudhuri Lab.  It can be forked into the lab GitHub account (https://github.com/immunogenomics) without counting against the lab's private repo quota.*

## Project Notes

#### Data source: GaP Long Island Registry (http://www.gapregistry.org/)

#### Platform: Illumina Global Screening Array (genotype array)

#### Script purpose(s):

1. 1_genotypeQC_gap2017.sh
	1. Process and QC GaP genotype data, do PCA.
	2. Merge with cleaned genotypes from 1000 Genomes, do PCA.

2. 2_pcaAncestryGenderAge_gap2017.R
	1. Plot PCA results (adapted from Dr. Yang Luo).
	2. Determine European ancestry of GaP participants by calculating Mahalanobis distances from PC1 and PC2 for Europeans in 1000 Genomes.  Threshold: 2 SD.
	3. Analyze and plot gender, age, and number of European GaP individuals per category.
	4. Select candidates for Immunomics project.


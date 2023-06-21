# DACompare
# Language: R
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 2.0, R 4.0.0
# Dependencies:  [1] doParallel_1.0.17       iterators_1.0.14        foreach_1.5.2
# [4] dimRed_0.2.6            DRR_0.0.4               CVST_0.2-3
# [7] Matrix_1.3-2            penalizedLDA_1.1        kernlab_0.9-32
# [10] probFDA_1.0.1           clusterGeneration_1.3.7 plotly_4.10.1
# [13] mixOmics_6.14.1         ggplot2_3.4.0           lattice_0.20-45
# [16] MASS_7.3-58.1           Rcpp_1.0.9

PluMA plugin that compares multiple DA methods (PLSDA, PCA, PFDA, PLDA) across synthetic data clusters with noise.

Input is a tab-delimited file of keyword-value pairs:

nclust  number of clusters
ncomp   number of components
sepFrom Separation value start
sepTo   Separating value end
sepInc  Separation value increment (will test every value start to end, incrementing by increment)
sampFrom        number of sample start
sampTo  number of sample end
sampInc number of sample increment (will test number of samples from start to end, incrementing by increment
nRepetitions    number of repetitions
nSignals        number of signals
nNoise  amount of noise

Output will be four CSV files, using the user-specified prefix (one for each algorithm).


# Wavelet characterization of spatial pattern in allele frequency
A set of methods for using wavelet transforms to characterize spatial pattern in allele frequency

Currently incomplete but in progress.

Preprint here
https://www.biorxiv.org/content/10.1101/2022.03.21.485229v3




## Genome-wide wavelet dissimilarity
[This file "ExampleGenomicWavDissim2.R",](https://github.com/jesserlasky/WaveletSpatialGenetic/blob/main/ExampleGenomicWavDissim2.R), contains a script with example SNP and geographic coordinate datasets (included in this repository).
This shows how to use my code to calculate genome-wide wavelet dissimilarity, using the SNP matrix and the coordinates of samples. It also shows how to plot the mean dissimilarity among sites compared to a null, as well as a the local dissimilarity at each sampled site compared to a null, across a range of spatial scales.

## Scale-specific genetic variance test
SLiM simulations were conducted with linear selective gradient and a single distinct habitat patch. "SpaVarSel_linear_sigma1_sigK1_JRL.1.5.slim" has the representative linear gradient simulation. "SpaVarSel_patchy_sigma2_sigK0.5_JRL.1.5.slim" has the representative patch simulation.

An example dataset is included, resulting from a run of SpaVarSel_patchy_sigma2_sigK0.5_JRL.1.5.slim using a scale of mating and dispersal sigma = 0.5, and a strength of selection 0.125.

The script ScalSpGenVarTest_Example_1_S.R shows how to implement the scale-specific genetic variance test on this example dataset.

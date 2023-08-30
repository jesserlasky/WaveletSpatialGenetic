# WaveletSpatialGenetic
A set of methods for using wavelet transforms to characterize spatial pattern in allele frequency

Currently incomplete but in progress.

Preprint here
https://www.biorxiv.org/content/10.1101/2022.03.21.485229v1




## Genome-wide wavelet dissimilarity
[This file](https://github.com/jesserlasky/WaveletSpatialGenetic/blob/main/ExampleGenomicWavDissim2.R), "ExampleGenomicWavDissim2.R", contains a script with example SNP and geographic coordinate datasets (included in this repository).
This shows how to use my code to calculate genome-wide wavelet dissimilarity, using the SNP matrix and the coordinates of samples. It also shows how to plot the mean dissimilarity among sites compared to a null, as well as a the local dissimilarity at each sampled site compared to a null, across a range of spatial scales.

## Scale-specific genetic variation test
SLiM simulations were conducted with linear selective gradient and a single distinct habitat patch. "SpaVarSel_linear_sigma1_sigK1_JRL.1.5.slim" has the representative linear gradient simulation. "SpaVarSel_patchy_sigma2_sigK0.5_JRL.1.5.slim" has the representative patch simulation.


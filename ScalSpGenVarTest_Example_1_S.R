


source('JL_OutflankFuns5_Fst.R')
source("WavFuns5.R")
source("Likelihood functions for OutFLANK_JRL.R")

scales <- exp(seq(1, 2.5, length.out = 6))

#this is already LD thinned, and filtered for MAF < 0.1, and signed so major allele is 1 and minor allele is 0.
SNPsF2 <- read.csv('model_output_patchySelecGrad_1.5_sigma0.5_sigmaK20.125_K25_tdelta500_gen10000_1337254.csv')

#SNP information table
snz <- read.csv('~/Dropbox/jesse/WaveletQstFst/Results/SLiM_sims/SNPinfo_model_output_patchySelecGrad_1.5_sigma0.5_sigmaK20.125_K25_tdelta500_gen10000_1337254.csv')


tmp2samp <- SNPsF2[,1:2] #x - y spatial coordinates

SNPsF <- as.matrix(SNPsF2[,-(1:2)])

#distance matrix
DistMat <- as.matrix(dist(tmp2samp[,c('x', 'y')]))

locSD <- apply(SNPsF, 2, sd)# locus standard deviation
locSD <- matrix(locSD, byrow = T, nrow = nrow(SNPsF), ncol = ncol(SNPsF))

#empty matrix for scaled-specific genetic variances
scalWavVar <- matrix(NA, nrow = length(scales), ncol = ncol(SNPsF))

outFtail <- 0.025

jl_out <- list()


for(i in 1:length(scales)){ #run test for all scales

		tmpW <- twavf_allOmega_b(SNPsF, OmegaDist = DistMat, s = scales[i]) #SNPs should be matrix

		scalWavVar[i, ] <- apply(tmpW/locSD, 2, var)# locus variance in wavelet

		jl_out[[i]] <- jl_iterator(scalWavVar[i, ], 25/scales[i], quantile(scalWavVar[i, ], outFtail), quantile(scalWavVar[i, ], 1 - outFtail))

	}



jl_out[[3]][[3]][snz$SNPeff != 0] #selected SNP p-values for 3rd scale

rank(jl_out[[3]][[3]])[snz$SNPeff != 0] #selected SNP ranks for 3rd scale (s = 4.95, approximately half the habitat patch size)

jl_out[[3]][[4]][snz$SNPeff != 0] #selected SNP q-values for 3rd scale









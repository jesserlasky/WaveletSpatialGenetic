
#this is modified from Whitlock & Lotterhos https://www.jstor.org/stable/10.1086/682949
#make the iterator function
#takes as input the wavlet transformed values, starting degrees of freedom

jl_iterator <- function(FSTvec, Nsamp_q, SmallestFstLimit, LargestLimit){
	#give Nsamp_q as landscape size / scale
	#FSTvec is the vector of wavelet variances, scaled
count <- 0

	keepGoing <- T

  while(keepGoing){

    count <- count+1
    if(count>19) {
      keepGoing <- FALSE  
      writeLines("Exceeded iteration maximum.") ###Try with increased maximum value for count two lines above.
    }

if(count == 1){ tmpFstbar <- mean(FSTvec)

	}else{


		tmpFstbar <- mean(FSTvec[FSTvec > SmallestFstLimit & qz > 0.05])

	}


	propEff <- EffectiveNumberSamplesMLE(FstVect = FSTvec, Fstbar = tmpFstbar, NumberOfSamples = Nsamp_q, SmallestFstInTrimmedList = SmallestFstLimit, LargestFstInTrimmedList = LargestLimit)

		#get qvalues
		pz <- pchisq((propEff - 1) * (FSTvec / tmpFstbar), df = propEff - 1, lower.tail = F)

		qz <- p.adjust(pz, method = 'fdr')

		if(count == 1) qz_old <- qz

		if(count > 1 & all.equal(qz > 0.05, qz_old > 0.05) == T){

			keepGoing <- F

		}else{
			qz_old <- qz
		}


	}

	list(propEff, tmpFstbar, pz, qz) #items are: inferred number of populations, mean wavelet transform for trimmed distribution, p values, q values 

}












#############		big function for simulation verification to take slim output file tables and generate wavelet fst assessment


simAssWavQ <- function(slimfile, scales, nsampz, outFtail, isNeut){
	#slimfile is the full path to the slim output file
	#scales are the scales to analyse
	#nsampz are number of sampled individuals
	#outFtail is what prroportion of the tail to remove
	#isNeut is logical to indicate whether it is neutral, or if there are loci truly under selection

library('sp')
library('slimr')
library('raster')

source("~/Dropbox/jesse/WaveletQstFst/Code/full_output.R")

source("~/Dropbox/jesse/WaveletQstFst/Code/MySLiM_Rfuns.R")

source("~/Dropbox/jesse/WaveletQstFst/Code/FstQst_WavFuns2.R")

source("~/Dropbox/jesse/WaveletQstFst/Code/Likelihood functions for OutFLANK.R")


tmp <- read_slim(slimfile)

if(isNeut == F) tmp2 <- read_individuals(tmp, colnamez =  c("indiv_id",    "sex", "genome1_id", "genome2_id", "x", "y", "phenotype", "age")) 
if(isNeut == T) tmp2 <- read_individuals(tmp, colnamez =  c("indiv_id",    "sex", "genome1_id", "genome2_id", "x", "y", "age")) # some are neutral


cat('slim file read', '\n')

nindiv <- nrow(tmp2)

nsamp <- nsampz


#apply grid to sample individuals
tmp2$xG <- round(tmp2$x / sqrt((25^2)/nsamp))
tmp2$yG <- round(tmp2$y / sqrt((25^2)/nsamp))
tmp2$Gr <- paste(tmp2$xG, tmp2$yG)


sampind <- c()

for(i in 1:length(unique(tmp2$Gr))){
	
	tmpI <- sample(tmp2$indiv_id[tmp2$Gr == unique(tmp2$Gr)[i]], 1)
	
	sampind <- c(sampind, tmpI)
	
	}

if(length(sampind) > nsamp) sampind <- sample(sampind, nsamp)

if(length(sampind) < nsamp){

	poz <- tmp2$indiv_id[! tmp2$indiv_id %in% sampind]
	
	sampind <- c(sampind, sample(poz, nsamp - length(sampind)))
	}

tmp2samp <- tmp2[match(sampind, tmp2$indiv_id) , ]

SNPs <- jl_getSLiM_SNPs(slimsimobj = tmp, sampind = sampind, ind.colnamez = c("indiv_id",    "sex", "genome1_id", "genome2_id", "x", "y", "phenotype", "age"))#, minMAF = 0.1) #200 individuals

cat('SNPs read', '\n')

mutz <- read_mutations(tmp)

#filter for minor allele frequency
SNPsF <- SNPs[,(colSums(SNPs) >= (nsamp * 0.1 * 2)) & (colSums(SNPs) <= (nsamp * 0.9 * 2))]

rm('SNPs')
gc()


#calculate minor allele freqncy
SNPsMAF <- colSums(SNPsF)/(nsamp * 2)

SNPsMAF[SNPsMAF > 0.5] <- 1 - SNPsMAF[SNPsMAF > 0.5]


SNPeff <- mutz$s[match(as.numeric(colnames(SNPsF)), unlist(mutz$mut_id))]



SNPpos <- mutz$pos[match(as.numeric(colnames(SNPsF)), unlist(mutz$mut_id))]


###For each SNP get wavelet diff at each sample location
#rsnps <- sample(1:ncol(SNPsF), 1e3) #for subsetting
rsnps <- 1:ncol(SNPsF)

locSD <- apply(SNPsF, 2, sd)# locus standard deviation

#position of non-neutral snps
causalpos <- SNPpos[SNPeff != 0]

scalWavVar <- matrix(NA, nrow = length(scales), ncol = ncol(SNPsF))
jl_out <- list()
t1z <- c()
t1z.q <- c()
pz <- c()
qz <- c()
rankz <- c()
top25pos <- c()




for(i in 1:length(scales)){

		tmpW <- apply(SNPsF[,rsnps], 2, twavf_allOmega, Omega = as.matrix(tmp2samp[,c('x', 'y')]), s = scales[i]) #wavelet difference

		scalWavVar[i, ] <- apply(tmpW/locSD, 2, var)# locus variance in wavelet

		jl_out[[i]] <- jl_iterator(scalWavVar[i, ], 25/scales[i], quantile(scalWavVar[i, ], outFtail), quantile(scalWavVar[i, ], 1 - outFtail))



#type 1 error
		t1z <- c(t1z, sum(jl_out[[i]][[3]][SNPeff == 0] < 0.05) / sum(SNPeff == 0))
		t1z.q <- c(t1z.q, sum(jl_out[[i]][[4]][SNPeff == 0] < 0.05) / sum(SNPeff == 0))

		top25pos <- rbind(top25pos, SNPpos[order(rank(jl_out[[i]][[3]]))][1:25])
##need to order ##FIX##

		if(sum(SNPeff != 0) > 1){ #multiple causal SNPs
				pz <- rbind(pz, jl_out[[i]][[3]][SNPeff != 0])
				qz <- rbind(qz, jl_out[[i]][[4]][SNPeff != 0])
				rankz <- rbind(rankz, rank(jl_out[[i]][[3]])[SNPeff != 0])

			}else{
				if(sum(SNPeff != 0) == 1){ #only 1 causal SNP
				 pz <- c(pz, jl_out[[i]][[3]][SNPeff != 0])
				 qz <- c(qz, jl_out[[i]][[4]][SNPeff != 0])
				 rankz <- c(rankz, rank(jl_out[[i]][[3]])[SNPeff != 0])

				}
			}

			cat('scale = ', round(scales[i], 2), '\n')

		}

phenCorr <- NA
if(isNeut == F) phenCorr <- cor(tmp2$y, tmp2$phenotype)

list(slimfile, 't1 p<0.05' = t1z, 't1 q<0.05' = t1z.q, 'causal_p' = pz, 'causal_q' = qz, 'causal_rank' = rankz, 'N_causal_SNPs' = sum(SNPeff != 0), 'phen_corr_y' = phenCorr, 'top25pos' = top25pos, 'causalpos' = causalpos)



}








#This is a set of functions for implementing a Difference of Gaussians wavelet on population genomic data.




kgauss_b <- function(OmegaDist, ab, s){

exp(-((OmegaDist[ab,]/s)^2)/2) #where ab is the focal row of omega

	}



 



etaV_b <- function(OmegaDist, ab, s){ #etaV is for vector of results one for each row in Omega

	tmp <- kgauss_b(OmegaDist, ab, s)

	tmp / sum(tmp)

	}




psiV_b <- function(OmegaDist, ab, s, beta){ #psiV is for vector of results one for each row in Omega
	etaV_b(OmegaDist, ab, s) - etaV_b(OmegaDist, ab, s*beta)
}





#this function gets the wavelet filter  for location row ab scale s
twavFilterOnly <- function(ab, OmegaDist, s){ 

	psis <- psiV_b(OmegaDist, ab, s, beta = 1.87) ########do this only once (instead of for every SNP)
	
	psis
	
}



twavFilterOnly_h <- function(psis) sqrt(sum(psis^2)) 






#this is applied to 'f_uvMat' a matrix of allele frequencies where rows are samples and columns are SNPs, a distance matrix of geographic distance among samples 'OmegaDist', and the spatial scale of analysis 's'
#this will give the wavelet transform for each SNP at each location in Omega for scale s.
twavf_allOmega_b <- function(f_uvMat, OmegaDist, s){ 	

	filtOnly <- sapply(1:nrow(OmegaDist), twavFilterOnly, OmegaDist = OmegaDist, s = s) #rows in Omega become columns in the output filtOnly
	norms_h <- apply(filtOnly, 2, twavFilterOnly_h)

	# f_uvMat has rows for samples, columns are SNPs.
	rawWavM <- (t(filtOnly) %*% (f_uvMat)) #/ norms_h

 	sweep(rawWavM, 1, norms_h, `/`) #sum are NA because norms_h is zero - i.e. there is no variation observable at the scale locally

	} 








gWavDiss <- function(f_uvMat, OmegaDist, scales, nperm){ #new function to be created


	emp_wavGdis <- c()
	null_wavGdis <- list()

	for(daS in 1:length(scales)){

		tmpW <- twavf_allOmega_b(f_uvMat, OmegaDist = OmegaDist, s = scales[daS]) #wavelet difference


		wavGdis <- sqrt(rowSums(tmpW^2)) #updated - added sqrt() - for v5

		emp_wavGdis <- cbind(emp_wavGdis, wavGdis) 

		tmp_null_wavGdis <- c()

		cat('scale', daS, '\n')

		for(j in 1:nperm){

			rindiv <- sample(1:nrow(f_uvMat)) #for permuting the distance matrix

			tmpWn <- twavf_allOmega_b(f_uvMat, OmegaDist = OmegaDist[rindiv, rindiv], s = scales[daS])

			wavGdisn <- sqrt(rowSums(tmpWn^2)) #updated - added sqrt() - for v5

			tmp_null_wavGdis <- cbind(tmp_null_wavGdis, wavGdisn) 

			cat(j)

		}

		null_wavGdis[[daS]] <- tmp_null_wavGdis
	}


	emp_wavGdis[!is.finite(emp_wavGdis)] <- NA

	CIz <- matrix(nrow = length(scales), ncol = 2)

	for(i in 1:length(scales)){
		null_wavGdis[[i]][!is.finite(null_wavGdis[[i]])] <- NA

		CIz[i,] <- quantile(apply(null_wavGdis[[i]], 2, mean, na.rm = T), probs = c(0.025, 0.975), na.rm = T)
		}

	list('obs_wavelets' = emp_wavGdis, 'null_wavelets' = null_wavGdis, 'null_95bounds' = CIz)	

	}








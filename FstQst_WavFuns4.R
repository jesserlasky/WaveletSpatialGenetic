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
	#h <- sqrt(sum(psis^2)) #normalization - also do once ####
	
	#sum(psis * f_uv) / h
	
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







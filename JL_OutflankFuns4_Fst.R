
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













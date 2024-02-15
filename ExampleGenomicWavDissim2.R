
#this script will show an example of getting genome wide dissimilarity

library('RColorBrewer')

source("WavFuns5.R")

tmp2samp <- read.csv('ExampleSampleLocs.1.csv')

nsamp <- nrow(tmp2samp)

SNPsF <- read.csv('ExampleSampleSNPs.1.csv', header = F)
#I have not yet tested this with missing allele calls. Currently I recommend using it with imputed data.

#the SNP matrix is just a matrix with alternate allele counts for each individual. Rows are samples, columns are SNPs

SNPsF <-as.matrix(SNPsF)


#convert so major allele is always 1 
SNPsF[,colSums(SNPsF) < (nsamp * 0.5 * 2)] <- 2 - SNPsF[,colSums(SNPsF) < (nsamp * 0.5 * 2)]####
######


#calculate a distance matrix between sampling locations
DistMat <- as.matrix(dist(tmp2samp[,c('x', 'y')]))




# this is choosing the minimum scale for analysis, but there is a bit of an art to this. if you have zeroes in the off diagonals of the distance matrix you wouldn't want to choose the minimum scale in this manner, because it may choose zero, which doesn't work.
#I recommend aiming for a minimum scale that is not so small as to yield very few pairs of locations that close together, but not so large as to completely miss interesting things that may be going on.
minscal <- quantile(DistMat[upper.tri(DistMat)], 0.001)

maxscal <- quantile(DistMat[upper.tri(DistMat)], 0.99) / 2

sclz <- exp(seq(log(minscal), log(maxscal), length.out = 15)) #these are the scales of analysis to implement. I recommend starting somewhere above the minimum distance between pairs, and stopping at a scale maybe half the landscape extent. These are still a bit arbitrary and merit further investigation.


set.seed(1) #for reproducibility


#calculate the genome-wide wavelet dissimliarity for the selected scales, using permutations to generate a null ('nperm' gives the number of permutations)
gWavDiss_out <- gWavDiss(SNPsF, OmegaDist = DistMat, scales = sclz, nperm = 100) 


#get y-limits for plot
tmpylim <- c(min(apply(gWavDiss_out[['obs_wavelets']], 2, mean, na.rm = T)), max(apply(gWavDiss_out[['obs_wavelets']], 2, mean, na.rm = T)))



#plot the change in mean genomic wavelet dissimilarity across scales, compared to null

plot.new()
plot.window(xlim = range(sclz), ylim = tmpylim, log = 'xy')

polygon(c(sclz, rev(sclz)), c(gWavDiss_out[['null_95bounds']][,1], rev(gWavDiss_out[['null_95bounds']][,2])), col = gray(0.75), border = NA)
text(5, 270, 'Null', pos = 4, col = gray(0.75))

lines(sclz, apply(gWavDiss_out[['obs_wavelets']], 2, mean, na.rm = T))

text(2.1, 1e3, 'Observed', pos = 4)

points(sclz, apply(gWavDiss_out[['obs_wavelets']], 2, mean, na.rm = T))


axis(1)
axis(2)

title(xlab = 's (spatial scale)')
title(ylab = expression(Mean~of~wavelet~dissimilarities~D[list(a,b)]^{wav} ~(s)), line = 2)






######## now looking at individual locations deviation from null ###########

#these are some functions for getting color ramps scaled to a range of values. It is inefficient and inelegant but it works.
mepal <- colorRampPalette(rev(brewer.pal(11, 'RdBu')[c(2:4, 8:11)]))

scalefun_col2 <- function(x, Ncol, Max = NULL, Min = NULL){
	if(is.null(Max)){
			tmp <- (x - min(x, na.rm = T)) / max(x - min(x, na.rm = T), na.rm = T)
			}else{
			tmp <- (x - Min) / (Max - Min)

			}
	round(tmp * (Ncol - 1)) + 1
	}


#need to fix bug when i make 4 one the scales::
scalecolz <- c(5, 6, 10, 15) #subset of scales to plot



par(mfrow = c(2,2)) #this needs adjusted depending on how many scales you want to plot

for(scl in 1:length(scalecolz)){

	scalecol <- scalecolz[scl]


	plot(tmp2samp[,c('x', 'y')])

	par(xpd = T)
	lines(c(1, 1 + sclz[scalecol]), c(27,27)) #this plots the scale as a line, but it has to be in the same units as the map coordinates. If your scales are in km but the map is lat/lon, this will not work, and instead you need to plot these lines in units of longitude degrees.
	text(1, 29, 'scale')
	par(xpd = F)



	#hypothesis test for individual locations
	sclCI <- apply(gWavDiss_out[['null_wavelets']][[scalecol]], 1, quantile, probs = c(0.025, 0.975), na.rm = T)

	redz <- scalefun_col2(gWavDiss_out[['obs_wavelets']][,scalecol], Ncol = 40, Max = max(gWavDiss_out[['obs_wavelets']][,scalecol][gWavDiss_out[['obs_wavelets']][,scalecol] > sclCI[2,]], na.rm = T), Min = min(gWavDiss_out[['obs_wavelets']][,scalecol][gWavDiss_out[['obs_wavelets']][,scalecol] > sclCI[2,]], na.rm = T)) + 60 #the signficantly high dissimilarity locations
	redz[redz < 0] <- NA

	points(tmp2samp[,c('x', 'y')], cex = c(0, 1.5)[1 + (gWavDiss_out[['obs_wavelets']][,scalecol] > sclCI[2,])], lwd = 1.5, 
		col = mepal(100)[redz], 
		pch = 19) #this will throw a warning if there are no locations where wavelet dissimilarity is more than expected


	points(tmp2samp[,c('x', 'y')], cex = c(0, 1.5)[1 + (gWavDiss_out[['obs_wavelets']][,scalecol] < sclCI[1,])], lwd = 1.5, col = mepal(100)[scalefun_col2(gWavDiss_out[['obs_wavelets']][,scalecol], Ncol = 40, Max = max(gWavDiss_out[['obs_wavelets']][,scalecol][gWavDiss_out[['obs_wavelets']][,scalecol] < sclCI[1,]], na.rm = T), Min = min(gWavDiss_out[['obs_wavelets']][,scalecol][gWavDiss_out[['obs_wavelets']][,scalecol] < sclCI[1,]], na.rm = T)) + 2], pch = 19) #this will throw a warning if there are no locations where wavelet dissimilarity is less than expected

	}
###










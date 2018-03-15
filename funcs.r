library(coin)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


pML <- function(chrList, IDs){
	## population allele freqs, simple ML estimate
	outlst <- list()
	for(i in 1:length(chrList)){
		G <- as.matrix(chrList[[i]])
		L <- dim(G)[1]
		N <- length(unique(IDs[,2]))
		p <- matrix(NA, nrow= L, ncol= N)
		for(j in 1:L){
		    p[j,] <- tapply(X= G[j,], INDEX= IDs[,2], mean)/2
		}
	colnames(p) <- unique(IDs[,2])
	outlst[[i]] <- p
	}
	return(outlst)
}


getAFdiff <- function(p){
	# allele frequency diffs of a single pML df (use loop for pML list)
	dat <- p
	pDiff <- as.data.frame(matrix(NA, nrow= dim(p)[1], ncol= 3))
	colnames(pDiff) <- c('cgal_mros', 'pach_mros', 'cgal_pach')
	pDiff[,1] <- dat[,1]-dat[,2]
	pDiff[,2] <- dat[,3]-dat[,2]
	pDiff[,3] <- dat[,1]-dat[,3]
	return(pDiff)
}



getPquant <- function(pDiff, quantil = 0.98){
	# get the indices for outliers given xth quantile 
	dat <- abs(pDiff)
	quantPls <- list()
	for(i in 1:dim(pDiff)[2]){
		quant <- quantile(dat[,i], prob= quantil, type= 8)
		quantPls[[i]] <- which(dat[,i] > quant)
	}
	names(quantPls) <- c('cgal_mros', 'pach_mros', 'cgal_pach')
	return(quantPls)
}

getFstGLquant <- function(pFst, quantil = 0.98){
	# get the indices for outliers given xth quantile 
	#dat <- abs(pFst)
	quantPls <- list()
	for(i in 1:length(pFst)){
		quant <- quantile(pFst[[i]], prob= quantil, type= 8, na.rm= T)
		quantPls[[i]] <- which(pFst[[i]] > quant)
	}
	names(quantPls) <- c('cgal_mros', 'pach_mros', 'cgal_pach')
	return(quantPls)
}



getFquant <- function(fst, quantil = 0.98){
	quantFls <- list()
	for(i in 1:length(fst)){
			dat <- fst[[i]]$WEIR_AND_COCKERHAM
			quant <- quantile(dat, prob= quantil, na.rm= T, type= 8)
			quantFls[[i]] <- which(dat >= quant)
	}
	names(quantFls) <- c('cgal_mros', 'pach_mros', 'cgal_pach')
	return(quantFls)
}


#subFstQls <- function(fQls, n= 8000){
#	for(i in 1:length(fstQls)){




is.integer0 <- function(x){
	is.integer(x) && length(x) == 0L
}

getVarSites <- function(ps, siteIX, specX1, specX2){
	## This function takes as input the genotype data.frame of loci by individuals, a vector of indices of target loci within said df (or NULL, if not supplying those), and the two respective indices for individuals.
	## OUTPUT is a list of 2 lists, with 1) the 2 subset dfs with all variable loci (sd != 0), respectively, and 2) the indices of those remaining loci (meant to be used in the getQr* functions)

	if(!is.null(siteIX)){
		ps <- ps[siteIX,]
	}
	invar1 <- which(apply(ps[,specX1], 1, sd) == 0)
	invar2 <- which(apply(ps[,specX2], 1, sd) == 0)
	if(is.integer0(invar1) && is.integer0(invar2)){
		pps1 <- ps[,specX1]
		pps2 <- ps[,specX2]
		ixSites1 <- siteIX
		ixSites2 <- siteIX
	}
	else if(is.integer0(invar1)){
		pps1 <- ps[, specX1]
		pps2 <- ps[-invar2,specX2]
		ixSites1 <- siteIX
		ixSites2 <- siteIX[!(siteIX %in% invar2)]
	}
	else if(is.integer0(invar2)){
		pps1 <- ps[-invar1,specX1]
		pps2 <- ps[, specX2]
		ixSites1 <- siteIX[!(siteIX %in% invar1)]
		ixSites2 <- siteIX
	}
	else{
		pps1 <- ps[-invar1, specX1]
		pps2 <- ps[-invar2, specX2]
		ixSites1 <- siteIX[!(siteIX %in% invar1)]
		ixSites2 <- siteIX[!(siteIX %in% invar2)]
	}
	return(list(list(pps1, pps2), list(ixSites1, ixSites2)))
}


getQr_allChrSub <- function(chrls, pQls, pops, tChr = 2, subN= 5000){
	## get r^2 values of mean genotypes for outliers (s) and neutral sites (n) between the two given taxa, for:
	## 1) s-s, n-n, and s-n within chromosome (tChr)
	## 2) s-s, n-n, and s-n between tChr and all other chromosomes
	## output is a list of lists with length(list) == length(pops) and the resp. comparisons within each of these two lists (length = 6) 
	spec1 <- pops[1]
	spec2 <- pops[2]
	specs <- c(rep('cgal', 10), rep('mros', 10), rep('pach', 10))
	indXs <- c(1:10, 11:20, 21:30)
	specX1 <- which(specs == spec1)
	specX2 <- which(specs == spec2)
	chromnum <- c(2, 7, 10, 18, 21)
	chromnames <- paste('chr', chromnum, sep= '')
	chrX <- which(chromnum == tChr)
	gChr <- chrls[[chrX]]
	# s 
	selX <- pQls[[chrX]]							#getPquant(pDiff, qtt)#[[1]]
	ppairs <- names(selX)
	specPairX <- which(ppairs == paste(pops[1], pops[2], sep= '_'))
	selXX <- selX[[specPairX]]
	neuXX <- seq(1, dim(gChr)[1])[-selXX]
	#
	pps <- getVarSites(gChr, selXX, specX1, specX2)
	if(dim(pps[[1]][[1]])[1] > subN){
		ixS1 <- sample(nrow(pps[[1]][[1]]), subN, replace= F)
		ixS2 <- sample(nrow(pps[[1]][[2]]), subN, replace= F)
		pps1 <- pps[[1]][[1]][ixS1,]
		pps2 <- pps[[1]][[2]][ixS2,]
	}
	else{ 
		ixS1 <- pps[[2]][[1]]
		ixS2 <- pps[[2]][[2]]
		pps1 <- pps[[1]][[1]]
		pps2 <- pps[[1]][[2]]
	}

	LDmat1 <- cor(t(pps1))^2
	LDmat2 <- cor(t(pps2))^2

	# n
	nps1 <- dim(pps1)[1]	#dim(LDmat1)[1]
	nps2 <- dim(pps2)[1]	#dim(LDmat2)[1]
	#pn <- gChr[neuXX,]
	ppn <- getVarSites(gChr, neuXX, specX1, specX2)
	ppn1 <- ppn[[1]][[1]]
	ppn2 <- ppn[[1]][[2]]
	ixN1 <- ppn[[2]][[1]]
	ixN2 <- ppn[[2]][[2]]
	
	#
	snx1 <- sample(nrow(ppn1), nps1, replace= F)
	snx2 <- sample(nrow(ppn2), nps2, replace= F)
	pN1 <- ppn1[snx1,]
	pN2 <- ppn2[snx2,]
	ixN1 <- ixN1[snx1]
	ixN2 <- ixN2[snx2]
	LDmatN1 <- cor(t(pN1))^2
	LDmatN2 <- cor(t(pN2))^2

	### s - n 
	mSN1 <- matrix(nrow= nps1, ncol= nps1)
	mSN2 <- matrix(nrow= nps2, ncol= nps2)
	for(i in 1:nps1){
		selloc1 <- pps1[i,]
		mSN1[i,] <- cor(t(pN1), as.numeric(selloc1))^2
	}
	for(i in 1:nps2){
	  	selloc2 <- pps2[i,]
		mSN2[i,] <- cor(t(pN2), as.numeric(selloc2))^2
	}
	mSN1 <- as.numeric(mSN1)[seq(2, length(mSN1)-nps1, 2)]
	mSN2 <- as.numeric(mSN2)[seq(2, length(mSN2)-nps2, 2)]

	#### diff chroms (s-s)
	gChrOthers <- chrl[-chrX]
	pQchrXs <- pQls[-chrX]		# indices for all "other" chroms, with all 3 specPairs
	pQchrOthers <- lapply(pQchrXs, '[[', specPairX)		# list of indices for afDiff outliers per Chr, for the resp. pop-pair
	names(pQchrOthers) <- names(chrl[-chrX])
	t.Seleck <- list()
	t.Neuter <- list()
	for(i in 1:length(pQchrOthers)){
		t.Seleck[[i]] <- gChrOthers[[i]][pQchrOthers[[i]],]
		t.Neuter[[i]] <- gChrOthers[[i]][-pQchrOthers[[i]],]
	}
	allSothers <- do.call("rbind", t.Seleck)
	ppsOthers <- getVarSites(allSothers, siteIX= NULL, specX1, specX2)
	ppsOth1 <- ppsOthers[[1]][[1]]
	ppsOth2 <- ppsOthers[[1]][[2]]
	pSothX1 <- sample(nrow(ppsOth1), nps1, replace= F)
	pSothX2 <- sample(nrow(ppsOth2), nps2, replace= F)

	SSothers1 <- matrix(nrow= nps1, ncol= nps1)
	SSothers2 <- matrix(nrow= nps2, ncol= nps2)
	for(i in 1:nps1){
		selloc1 <- pps1[i,]
		SSothers1[i,] <- cor(t(ppsOth1[pSothX1,]), as.numeric(selloc1))^2
	}
	for(i in 1:nps2){
		selloc2 <- pps2[i,]
		SSothers2[i,] <- cor(t(ppsOth2[pSothX2,]), as.numeric(selloc2))^2
	}
	SSothers1 <- as.numeric(SSothers1)[seq(2, length(SSothers1)-nps1, 2)]
	SSothers2 <- as.numeric(SSothers2)[seq(2, length(SSothers2)-nps2, 2)]

	####    (s-n)
	allNothers <- do.call("rbind", t.Neuter)
	ppNothers <- getVarSites(allNothers, siteIX= NULL, specX1, specX1)
	ppNoth1 <- ppNothers[[1]][[1]]
	ppNoth2 <- ppNothers[[1]][[2]]
	pNothX1 <- sample(nrow(ppNoth1), nps1, replace= F)
	pNothX2 <- sample(nrow(ppNoth2), nps2, replace= F)
	#
	SNothers1 <- matrix(nrow= nps1, ncol= nps1)
	SNothers2 <- matrix(nrow= nps2, ncol= nps2)
	for(i in 1:nps1){
		selloc1 <- pps1[i,]
		SNothers1[i,] <- cor(t(ppNoth1[pNothX1,]), as.numeric(selloc1))^2
	}
	for(i in 1:nps2){
		selloc2 <- pps2[i,]
		SNothers2[i,] <- cor(t(ppNoth2[pNothX2,]), as.numeric(selloc2))^2
	}
	SNothers1 <- as.numeric(SNothers1)[seq(2, length(SNothers1)-nps1, 2)]
	SNothers2 <- as.numeric(SNothers2)[seq(2, length(SNothers2)-nps2, 2)]
	####    (n-n)  
	NNothers1 <- matrix(nrow= nps1, ncol= nps1)
	NNothers2 <- matrix(nrow= nps2, ncol= nps2)
	for(i in 1:nps1){
		neutloc1 <- pN1[i,]
		NNothers1[i,] <- cor(t(ppNoth1[pNothX1,]), as.numeric(neutloc1))^2
	}
	for(i in 1:nps2){
		neutloc2 <- pN2[i,]
		NNothers2[i,] <- cor(t(ppNoth2[pNothX2,]), as.numeric(neutloc2))^2
	}
	NNothers1 <- as.numeric(NNothers1)[seq(2, length(NNothers1)-nps1, 2)]
	NNothers2 <- as.numeric(NNothers2)[seq(2, length(NNothers2)-nps2, 2)]
	####

	spec1list <- list(LDmat1[upper.tri(LDmat1)], LDmatN1[upper.tri(LDmatN1)], mSN1, SSothers1, SNothers1, NNothers1)
	names(spec1list) <- c("LDmat", "LDmatN", "mSN", "SSothers", "SNothers", "NNothers")
	spec2list <- list(LDmat2[upper.tri(LDmat2)], LDmatN2[upper.tri(LDmatN2)], mSN2, SSothers2, SNothers2, NNothers2)
	names(spec2list) <- names(spec1list)
	LDlist <- list(spec1list, spec2list)
	names(LDlist) <- pops
	return(LDlist)
}

getQrDistIXsn <- function(chrls, pQls, pops, bPosls, tChr = 2, subN= 5000){
	## get r^2 values of mean genotypes  and their respective distances for outliers (s) and neutral sites (n) between the two given taxa
	## NOTE: only for the s, n, and s-n within the given chromosome (tChr)
	## 1) s-s, and n-n within chromosome (tChr) and distance matrices for the respective loci 
	## output is a list of lists with length(list) == length(pops) and the resp. comparisons within each of these two lists (length = 4) 
	spec1 <- pops[1]
	spec2 <- pops[2]
	specs <- c(rep('cgal', 10), rep('mros', 10), rep('pach', 10))
	indXs <- c(1:10, 11:20, 21:30)
	specX1 <- which(specs == spec1)
	specX2 <- which(specs == spec2)
	chromnum <- c(2, 7, 10, 18, 21)
	chromnames <- paste('chr', chromnum, sep= '')
	chrX <- which(chromnum == tChr)
	gChr <- chrls[[chrX]]
	# s 
	selX <- pQls[[chrX]]							#getPquant(pDiff, qtt)#[[1]]
	ppairs <- names(selX)
	specPairX <- which(ppairs == paste(pops[1], pops[2], sep= '_'))
	selXX <- selX[[specPairX]]
	neuXX <- seq(1, dim(gChr)[1])[-selXX]
	#
	pps <- getVarSites(gChr, selXX, specX1, specX2)
	
	# selXX, ixS1, ixS2, ixN1, ixN2
	
	
	if(dim(pps[[1]][[1]])[1] > subN){
		ixS1 <- sample(nrow(pps[[1]][[1]]), subN, replace= F)
		ixS2 <- sample(nrow(pps[[1]][[2]]), subN, replace= F)
		pps1 <- pps[[1]][[1]][ixS1,]
		pps2 <- pps[[1]][[2]][ixS2,]
	}
	else{ 
		ixS1 <- pps[[2]][[1]]
		ixS2 <- pps[[2]][[2]]
		pps1 <- pps[[1]][[1]]
		pps2 <- pps[[1]][[2]]
	}

	LDmat1 <- cor(t(pps1))^2
	LDmat2 <- cor(t(pps2))^2
	# n
	nps1 <- dim(pps1)[1]	#dim(LDmat1)[1]
	nps2 <- dim(pps2)[1]	#dim(LDmat2)[1]
	#pn <- gChr[neuXX,]
	ppn <- getVarSites(gChr, neuXX, specX1, specX2)
	ppn1 <- ppn[[1]][[1]]
	ppn2 <- ppn[[1]][[2]]
	ixN1 <- ppn[[2]][[1]]
	ixN2 <- ppn[[2]][[2]]
	
	#
	snx1 <- sample(nrow(ppn1), nps1, replace= F)
	snx2 <- sample(nrow(ppn2), nps2, replace= F)
	pN1 <- ppn1[snx1,]
	pN2 <- ppn2[snx2,]
	ixN1 <- ixN1[snx1]
	ixN2 <- ixN2[snx2]
	LDmatN1 <- cor(t(pN1))^2
	LDmatN2 <- cor(t(pN2))^2
#	
#	### s - n 
#	mSN1 <- matrix(nrow= nps1, ncol= nps1)
#	mSN2 <- matrix(nrow= nps2, ncol= nps2)
#	for(i in 1:nps1){
#		selloc1 <- pps1[i,]
#		mSN1[i,] <- cor(t(pN1), as.numeric(selloc1))^2
#	}
#	for(i in 1:nps2){
#	  	selloc2 <- pps2[i,]
#		mSN2[i,] <- cor(t(pN2), as.numeric(selloc2))^2
#	}
#	mSN1 <- as.numeric(mSN1)[seq(2, length(mSN1)-nps1, 2)]
#	mSN2 <- as.numeric(mSN2)[seq(2, length(mSN2)-nps2, 2)]
	bPos <- bPosls[[chrX]]
	bpS1 <- bPos[ixS1,2]
	bpS2 <- bPos[ixS2,2]
	bpN1 <- bPos[ixN1,2]
	bpN2 <- bPos[ixN2,2]
	distS1 <- dist(bpS1)
	distS2 <- dist(bpS2)
	distN1 <- dist(bpN1)
	distN2 <- dist(bpN2)


	return(list(list(LDmat1[upper.tri(LDmat1)], LDmatN1[upper.tri(LDmatN1)], distS1, distN1), list(LDmat2[upper.tri(LDmat2)], LDmatN2[upper.tri(LDmatN2)], distS2, distN2)))				#ixS1, ixN1), list(ixS2, ixN2))
}


getSubQld <- function(chrls, pQls, pops, tChr = 2, subN= 5000){
	## get r^2 values of mean genotypes for outliers (s) and neutral sites (n) between the two given taxa, for:
	## NOTE: only for s, n, and s-n within the given chromosome (tChr)
	## s-s, n-n, and s-n within chromosome (tChr)
	## output is a list of lists with length(list) == length(pops) and the resp. comparisons within each of these two lists (length = 3) 
	spec1 <- pops[1]
	spec2 <- pops[2]
	specs <- c(rep('cgal', 10), rep('mros', 10), rep('pach', 10))
	indXs <- c(1:10, 11:20, 21:30)
	specX1 <- which(specs == spec1)
	specX2 <- which(specs == spec2)
	chromnum <- c(2, 7, 10, 18, 21)
	chromnames <- paste('chr', chromnum, sep= '')
	chrX <- which(chromnum == tChr)
	gChr <- chrls[[chrX]]
	# s 
	selX <- pQls[[chrX]]							#getPquant(pDiff, qtt)#[[1]]
	ppairs <- names(selX)
	specPairX <- which(ppairs == paste(pops[1], pops[2], sep= '_'))
	selXX <- selX[[specPairX]]
	neuXX <- seq(1, dim(gChr)[1])[-selXX]
	#
	pps <- getVarSites(gChr, selXX, specX1, specX2)
	if(dim(pps[[1]][[1]])[1] > subN){
		ixS1 <- sample(nrow(pps[[1]][[1]]), subN, replace= F)
		ixS2 <- sample(nrow(pps[[1]][[2]]), subN, replace= F)
		pps1 <- pps[[1]][[1]][ixS1,]
		pps2 <- pps[[1]][[2]][ixS2,]
	}
	else{ 
		ixS1 <- pps[[2]][[1]]
		ixS2 <- pps[[2]][[2]]
		pps1 <- pps[[1]][[1]]
		pps2 <- pps[[1]][[2]]
	}
	#	
	LDmat1 <- cor(t(pps1))^2
	LDmat2 <- cor(t(pps2))^2
	# n
	nps1 <- dim(pps1)[1]	#dim(LDmat1)[1]
	nps2 <- dim(pps2)[1]	#dim(LDmat2)[1]
	#pn <- gChr[neuXX,]
	ppn <- getVarSites(gChr, neuXX, specX1, specX2)
	ppn1 <- ppn[[1]][[1]]
	ppn2 <- ppn[[1]][[2]]
	ixN1 <- ppn[[2]][[1]]
	ixN2 <- ppn[[2]][[2]]
	
	#
	snx1 <- sample(nrow(ppn1), nps1, replace= F)
	snx2 <- sample(nrow(ppn2), nps2, replace= F)
	pN1 <- ppn1[snx1,]
	pN2 <- ppn2[snx2,]
	ixN1 <- ixN1[snx1]
	ixN2 <- ixN2[snx2]
	LDmatN1 <- cor(t(pN1))^2
	LDmatN2 <- cor(t(pN2))^2

	### s - n 
	mSN1 <- matrix(nrow= nps1, ncol= nps1)
	mSN2 <- matrix(nrow= nps2, ncol= nps2)
	for(i in 1:nps1){
		selloc1 <- pps1[i,]
		mSN1[i,] <- cor(t(pN1), as.numeric(selloc1))^2
	}
	for(i in 1:nps2){
	  	selloc2 <- pps2[i,]
		mSN2[i,] <- cor(t(pN2), as.numeric(selloc2))^2
	}
	mSN1 <- as.numeric(mSN1)[seq(2, length(mSN1)-nps1, 2)]
	mSN2 <- as.numeric(mSN2)[seq(2, length(mSN2)-nps2, 2)]

	spec1list <- list(LDmat1[upper.tri(LDmat1)], LDmatN1[upper.tri(LDmatN1)], mSN1)
	names(spec1list) <- c("LDmat", "LDmatN", "mSN")
	spec2list <- list(LDmat2[upper.tri(LDmat2)], LDmatN2[upper.tri(LDmatN2)], mSN2)
	names(spec2list) <- names(spec1list)
	LDlist <- list(spec1list, spec2list)
	names(LDlist) <- pops
	return(LDlist)
}



plotFstPerSpec <- function(x, tchr){
	chromnum <- c(2, 7, 10, 18, 21)
	#specs <- c(rep('c. galanthus', 10), rep('m. rosina', 10), rep('pachinus', 10))
	chromnames <- paste('chrom', chromnum, sep= '')
	chromIX <- which(chromnum == tchr)
	chromN <- chromnames[chromIX]
	par(mfrow= c(3,1), oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	for(i in 1:length(x)){
		specsX <- names(x)
		dat <- x[[i]][[chromIX]]
		plot(dat$WEIGHTED_FST, type= 'l', ylab= expression(F[ST]), xlab= paste(specsX[i], chromN, sep= '   -   '), cex.lab= 1.4, cex= 1.3, ylim= c(0,1))
	}
}

plotVioChr <- function(tajDobj, tchrom){
	specNames <- c('c. galanthus', 'm. rosina', 'pachinus')
	chromnum <- c(2, 7, 10, 18, 21)
	chromnames <- paste('chrom', chromnum, sep= '')
	chromIX <- which(chromnum == tchrom)
	chromN <- chromnames[chromIX]
	#
	vioplot(na.omit(tajDobj$cgal[[chromIX]]$TajimaD), na.omit(tajDobj$mros[[chromIX]]$TajimaD), na.omit(tajDobj$pach[[chromIX]]$TajimaD), names= specNames, col= brewer.pal(3, 'Set1')[2])	#'white') 
}

plotVioSpec <- function(tajDobj, spec){
	specNames <- c('c. galanthus', 'm. rosina', 'pachinus')
	chromnum <- c(2, 7, 10, 18, 21)
	chromnames <- paste('chrom', chromnum, sep= '')
	specs <- c('cgal', 'mros', 'pach')
	specIX <- which(specs == spec)
	Obj <- tajDobj[[specIX]]
	#
	vioplot(na.omit(Obj$chrom2$TajimaD), na.omit(Obj$chrom7$TajimaD), na.omit(Obj$chrom10$TajimaD), na.omit(Obj$chrom18$TajimaD), na.omit(Obj$chrom21$TajimaD), names= chromnames, col= brewer.pal(3, 'Set2')[1])
	abline(h= 0, lty= 3)
}

# tajD for each chromosome 
plotTajimaD <- function(tajDobj, spec = 'cgal'){
	par(mfrow= c(5,1), oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	for(i in 1:5){
		specs <- c('cgal', 'mros', 'pach')
		specNames <- c('c. galanthus', 'm. rosina', 'pachinus')
		specIX <- which(specs == spec)
		Obj <- tajDobj[[specIX]][[i]]
		nBins <- dim(Obj)[1]
		if(i == 1){
			plot(seq(1, nBins, 1), Obj$TajimaD, type= 'l', xlab= chromnames[i], ylab= "Tajima's D", cex= 1.3, cex.axis= 1.4, main= specNames[specIX], cex.main= 1.5) 	
			abline(h= 0, lty= 3)
		}
		else if(i == 4){
			plot(seq(1, nBins, 1), Obj$TajimaD, type= 'l', xlab= chromnames[i], ylab= "Tajima's D", cex= 1.3, cex.axis= 1.4, ) 
			abline(h= 0, lty= 3)
			abline(v= Obj$BIN_START[which(Obj$BIN_START == 7e+05)], col= 'green')
		}
		else{
			plot(seq(1, nBins, 1), Obj$TajimaD, type= 'l', xlab= chromnames[i], ylab= "Tajima's D", cex= 1.3, cex.axis= 1.4, ) 
			abline(h= 0, lty= 3)
		}
	}
}

printQ <- function(tajDat, thresh= 2){
	dat <- tajDat$TajimaD
	lowerq = quantile(dat, na.rm= T)[2]
	upperq = quantile(dat, na.rm= T)[4]
	iqr = upperq - lowerq #Or use IQR(data)
	#Compute the bounds for an outlier, given thresh (1.5, 3,...)
	upper = (iqr * thresh) + upperq
	lower = lowerq - (iqr * thresh)
	return(tajDat[which(dat < lower | dat > upper),])
}

plotDxy <- function(x, tchr){
	chromnum <- c(2, 7, 10, 18, 21)
	#specs <- c(rep('c. galanthus', 10), rep('m. rosina', 10), rep('pachinus', 10))
	chromnames <- paste('chrom', chromnum, sep= '')
	chromIX <- which(chromnum == tchr)
	chromN <- chromnames[chromIX]
	par(mfrow= c(3,1), oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	#for(i in 1:length(x)){
		#specsX <- names(x)
	dat <- x[[chromIX]]
	plot(dat$dxy_cgal_mros, type= 'l', ylab= expression(D[XY]), xlab= 'cgal_mros', cex.lab= 1.4, cex= 1.3, ylim= c(0,1), main= chromN, cex.main= 1.6)
	#abline(h= 0, lty= 3)
	plot(dat$dxy_mros_pach, type= 'l', ylab= expression(D[XY]), xlab= 'pach_mros', cex.lab= 1.4, cex= 1.3, ylim= c(0,1))
	#abline(h= 0, lty= 3)
	plot(dat$dxy_cgal_pach, type= 'l', ylab= expression(D[XY]), xlab= 'cgal_pach', cex.lab= 1.4, cex= 1.3, ylim= c(0,1))
	#abline(h= 0, lty= 3)
}

plotFst <- function(x, tchr){
	chromnum <- c(2, 7, 10, 18, 21)
	#specs <- c(rep('c. galanthus', 10), rep('m. rosina', 10), rep('pachinus', 10))
	chromnames <- paste('chrom', chromnum, sep= '')
	chromIX <- which(chromnum == tchr)
	chromN <- chromnames[chromIX]
	par(mfrow= c(3,1), oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	#for(i in 1:length(x)){
		#specsX <- names(x)
	dat <- x[[chromIX]]
	plot(dat$Fst_cgal_mros, type= 'l', ylab= expression(F[ST]), xlab= 'cgal_mros', cex.lab= 1.4, cex= 1.3, ylim= c(0,1), main= chromN, cex.main= 1.6)
	#abline(h= 0, lty= 3)
	plot(dat$Fst_mros_pach, type= 'l', ylab= expression(F[ST]), xlab= 'pach_mros', cex.lab= 1.4, cex= 1.3, ylim= c(0,1))
	#abline(h= 0, lty= 3)
	plot(dat$Fst_cgal_pach, type= 'l', ylab= expression(F[ST]), xlab= 'cgal_pach', cex.lab= 1.4, cex= 1.3, ylim= c(0,1))
	#abline(h= 0, lty= 3)
}

plotPi <- function(x, tchr){
	chromnum <- c(2, 7, 10, 18, 21)
	#specs <- c(rep('c. galanthus', 10), rep('m. rosina', 10), rep('pachinus', 10))
	chromnames <- paste('chrom', chromnum, sep= '')
	chromIX <- which(chromnum == tchr)
	chromN <- chromnames[chromIX]
	par(mfrow= c(3,1), oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
	#for(i in 1:length(x)){
		#specsX <- names(x)
	dat <- x[[chromIX]]
	plot(dat$pi_cgal, type= 'l', ylab= expression(pi), xlab= 'cgal', cex.lab= 1.4, cex= 1.3, ylim= c(0,1), main= chromN, cex.main= 1.6)
	#abline(h= 0, lty= 3)
	plot(dat$pi_mros, type= 'l', ylab= expression(pi), xlab= 'mros', cex.lab= 1.4, cex= 1.3, ylim= c(0,1))
	#abline(h= 0, lty= 3)
	plot(dat$pi_pach, type= 'l', ylab= expression(pi), xlab= 'pach', cex.lab= 1.4, cex= 1.3, ylim= c(0,1))
	#abline(h= 0, lty= 3)
}

getIndyP <- function(dat, distro= 'asymptotic', tstat= 'maximum', alt= 'less'){
	pops <- names(dat)
	s1 <- dat[[1]]$LDmat
	s2 <- dat[[2]]$LDmat
	n1 <- dat[[1]]$LDmatN
	n2 <- dat[[2]]$LDmatN
	#
	tst1 <- independence_test(log(s1) ~ log(n1), distribution= distro, teststat= tstat, alternative= alt)
	tst2 <- independence_test(log(s2) ~ log(n2), distribution= distro, teststat= tstat, alternative= alt)

	p1 <- pvalue(tst1)
	p2 <- pvalue(tst2)
	out <- list(p1, p2)
	names(out) <- pops
	return(out)
}


hist_logX <- function(dat, fname, ylim= "", breaks= 100){
	#yLim <- c(1,5e+05)
	pops <- names(dat)	
	pdf(fname)		#'cgmr2log2.pdf')
	par(mfrow= c(4,3))		#, mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))		#oma=c(4.5,4.5,4.5,4.5), 
	for(i in 1:2){
		#pops <- names(dat)
		hist(log(dat[[i]][[1]]), breaks= breaks, main= paste(pops[i], "selected sites"), xlab= expression(paste('log r'^'2')))
		hist(log(dat[[i]][[2]]), breaks= breaks, main= paste(pops[i], "neutral sites"), xlab= expression(paste('log r'^'2')))
		hist(log(dat[[i]][[3]]), breaks= breaks, main= paste(pops[i], "selected vs. neutral sites"), xlab= expression(paste('log r'^'2')))
		hist(log(dat[[i]][[4]]), breaks= breaks, main= paste(pops[i], "selected sites w/ 'other' s"), xlab= expression(paste('log r'^'2')))
		hist(log(dat[[i]][[5]]), breaks= breaks, main= paste(pops[i], "neutral sites w/ 'other' n"), xlab= expression(paste('log r'^'2')))
		hist(log(dat[[i]][[6]]), breaks= breaks, main= paste(pops[i], "selected w/ 'other' n sites"), xlab= expression(paste('log r'^'2')))
	}
	dev.off()
}



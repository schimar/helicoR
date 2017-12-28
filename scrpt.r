
setwd('/home/schimar/flaxmans/bf/helico/var/subchr/')
library(RColorBrewer)
library(fields)
#library(vioplot)

chromnum <- c(2, 7, 10, 18, 21)
chromnames <- paste('chrom', chromnum, sep= '')



################# 		chromStats 		####################

cs <- read.table("~/flaxmans/bf/helico/var/chromStats.txt", header= F, sep= ' ')
cs <- cs[order(cs$V1),]
names(cs) <- c('chrom', 'nVar', 'dpMean', 'dpSD', 'mqMean', 'mqSD')


scafs <- read.table('../hmel2.5.30f4nScafPerChrom.txt', header= F)

sTable1 <- cbind(cs, scafs[,2])
colnames(sTable1) <- c(colnames(cs), "nScafs")

# barplots


#
par(mfrow= c(1,3))

bp <- barplot(cs[,2], names.arg= seq(1,21), xlab= 'chromosome', ylab= 'number of variants', cex.lab= 1.4)

text(x= bp, y= 2e5, labels= scafs[,2], pos= 1, cex= 1.4)


barx <- barplot(cs$dpMean, names.arg= seq(1, 21, 1), xlab= 'chromosome', ylab= 'mean depth of coverage (with sd)', ylim= c(0, 500), cex.lab= 1.4)
error.bar(barx,cs$dpMean, 1.96*cs$dpSD/10)

barQ <- barplot(cs$mqMean, names.arg= seq(1, 21, 1), xlab= 'chromosome', ylab= 'mean mapping quality (with sd)', ylim= c(0, 65), cex.lab= 1.4)
error.bar(barQ,cs$mqMean, 1.96*cs$mqSD/10)


########################

# nVar at MAF thresholds

vars <- read.table('hmel2.5.30f4nVarThresh.txt', header= T)
vars <- vars[order(vars$chrom),]


ths <- c(0.0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25)
#barplot(as.matrix(vars[,2:8]), beside= T, names.arg= ths, xlab= 'chromosomes per threshold (1-21)', ylab= 'number of variants')

# better per chromosome

barplot(as.matrix(t(vars[,2:8])), beside= T, legend= ths, xlab= 'threshold * chromosome', names.arg= seq(1,21,1), ylab= 'number of variants', col= brewer.pal(7, "Set3"), ylim= c(0, 1000000))
grid(nx= 0, ny= 5, col= 'grey40')



# load ind/pop info

ids <- read.csv('../../kronforst2013ids.txt', sep= '\t', header= T)
ids <- ids[c(1:10, 13:32),]	# throw out h665 and i02_210
ids[] <- lapply(ids, function(x) if(is.factor(x)) factor(x) else x) # drop old factors 
# 



#########################################################

## pca 

chr2 <- read.table('hmel2.5.30chr2est.txt', header= T)
chr7 <- read.table('hmel2.5.30chr7est.txt', header= T)
chr10 <- read.table('hmel2.5.30chr10est.txt', header= T)
chr18 <- read.table('hmel2.5.30chr18est.txt', header= T)
chr21 <- read.table('hmel2.5.30chr21est.txt', header= T)

chrl <- list(chr2, chr7, chr10, chr18, chr21)
names(chrl) <- paste0('chr', chromnum)


inds <- colnames(chr2)

kids <- read.table("../../kronforst2013ids_noheader.txt", sep= '\t', header= F)
kids <- kids[kids$V1 %in% inds,]

kids <- kids[match(inds, kids[,1]),]

specs <- c(rep('c. galanthus', 10), rep('m. rosina', 10), rep('pachinus', 10))

# for each chrom
#cov(as.matrix(cmn))
#prcomp(cmncov, center= T, scale.= F)

chromnum <- c(2, 7, 10, 18, 21)

# calculate covariance matrix for each subchr
covl <- list()
for(i in 1:length(chrl)){
	covl[[i]] <- cov(as.matrix(chrl[[i]]))
}

# calc. pca
pcl <- list()
for(i in 1:length(covl)){
	pcl[[i]] <- prcomp(covl[[i]], center= T, scale.= F)
}

# calc. percentage of variance explained
povs <- list()
for(i in 1:length(pcl)){
	pov <- pcl[[i]]$sdev^2/sum(pcl[[i]]$sdev^2)#[1:2]
	povs[[i]] <- pov[1:2]
}

#PoV <- pca$sdev^2/sum(pca$sdev^2)



cols <- brewer.pal(3, "Set1")
colvec <- c(rep(cols[1], 10), rep(cols[2], 10), rep(cols[3], 10))


par(mfrow= c(1,5))
for(i in 1:length(pcl)){
	plot(pcl[[i]]$x[,1], pcl[[i]]$x[,2], main= paste('chrom ', chromnum[i]), xlab= paste('PC 1 (', round(povs[[i]][1]*100, 1), '% expl. var.)'), ylab= paste('PC 2 (', round(povs[[i]][2]*100, 1), '% expl. var.)'), type= 'n', cex.lab= 1.3)
	text(pcl[[i]]$x[,1], pcl[[i]]$x[,2], colnames(chrl[[i]]), cex= 1, col= colvec)
	if(i == 1){
		legend('topright', legend= unique(specs), fill= cols, cex= 1.3)
	}
}


###### dapc 
library(adegenet)

chr2clust <- find.clusters(t(chr2), scale= F, max.n.clust= 10)
chr2clustk2 <- chr2clust	# with 25 PCs and 2 clusters
# plus k = {3, 4} 

dapchr2k2 <- dapc(t(chr2), scale= F, grp= chr2clustk2$grp)

#dapchr2k2 <- dapchr2		# with 5 PCs and 1 DF reteined

## plots (make sure to change here, when using k > 2)
compoplot(dapchr2k2, posi= 'bottomright', txt.leg= paste("Cluster", 1:2), lab= "", ncol= 1, xlab= "individuals", col=funky(2))

library(gplots)
heatmap.2(dapchr2k2$posterior, dendrogram= 'none', trace= 'none', colsep= c(1,2), rowsep= seq(1, 30, 1))

############ allele freqs   (estpEM)

#rchrl <- lapply(chrl, round)
lsp <- system("ls pops/", intern= T)
#fls.p <- lsp[grep(".p", lsp)]#[c(1:14, 16)]


fstNames <- unlist(strsplit(fls.fst, '.windowed.weir.fst'))
#for(i in 1:length(fls.fst)){




# mean genotypes 

ids <- read.csv('../../kronforst2013ids.txt', sep= '\t', header= T)
ids <- ids[c(1:10, 13:32),]	# throw out h665 and i02_210
ids[] <- lapply(ids, function(x) if(is.factor(x)) factor(x) else x) # drop old factors 


# split into pops

#dat <- chr2
#chr2ls <- list()
#for(i in 1:length(unique(ids$species))){
#	specIX <- which(ids$species == unique(ids$species)[i])
#	chr2ls[[i]] <- dat[,specIX]
#}
#names(chr2ls) <- c('cgal', 'mros', 'pach')
#
## get mean allele freqs
#
#pbars <- lapply(chr2ls, rowMeans)
#names(pbars) <- names(chr2ls)


## population allele freqs., simple ML estimate
#G <- as.matrix(round(chr2))

pls <- pML(chrl, ids)

# for single g (chrX)
G <- as.matrix(chr2)

L <- dim(G)[1]
N <- length(unique(ids[,2]))
p <- matrix(NA, nrow= L, ncol= N)
for(i in 1:L){
    p[i,] <- tapply(X= G[i,], INDEX= ids[,2], mean)/2
}
colnames(p) <- unique(ids[,2])

# allele frequency diffs 

# pDiff <- getAFdiff(p)

# for chr list
pDifls <- list()
for(i in 1:length(pls)){
	pDifls[[i]] <- getAFdiff(pls[[i]])
}

dat <- pDifls[[5]]
par(mfrow= c(1,3))
boxplot(abs(dat$cgal_mros), main= colnames(dat)[1], ylab= 'abs(afDiff)')
boxplot(abs(dat$pach_mros), main= colnames(dat)[2], ylab= 'abs(afDiff)')
boxplot(abs(dat$cgal_pach), main= colnames(dat)[3], ylab= 'abs(afDiff)')


vioplot(abs(dat$cgal_mros), abs(dat$pach_mros), abs(dat$cgal_pach))



# get the xth quantile (< outlier)

pQls <- list()
for(i in 1:length(pDifls)){
	pQls[[i]] <- getPquant(pDifls[[i]], 0.98)
}



# clines ( = inverse of afDiffs
#length(which(1/abs(pDiff$cgal_mros) > quantile(1/abs(pDiff$cgal_mros), prob= 1- 0.05)))


#1/abs(pDiff$cgal_mros)[which(1/abs(pDiff$cgal_mros) > quantile(1/abs(pDiff$cgal_mros), prob= 1- 0.05))]



# fst 
pbar <- apply(p, 1, mean)
vp <- apply(p, 1, var)
fst <- mean(vp)/mean(pbar * (1 - pbar)) ## ratio of averages

# pairwise fst  (cgal | mros)
lp <- p[,c(1, 2)]
pbar_lp <- apply(lp, 1, mean)
vp_lp <- apply(lp, 1, var)
fst_lp <- mean(vp_lp)/mean(pbar_lp * (1- pbar_lp))




#########################################################

## nucleotide diversity  (vcftools)


# read files
lspi <- system("ls pgStats/pi/", intern= T)
sites.pi <- lspi[grep("sites.pi", lspi)]#[c(1:14, 16)]
piName <- unlist(strsplit(sites.pi, '_'))[seq(1, 30, 2)]
piNames <- paste('pi_', piName, sep= 
for(i in 1:length(sites.pi)){
    #oname = paste(estname[i], "est", sep= "")
	fname <- paste('pgStats/pi/', sites.pi[i], sep= '')
	print(fname)
    assign(piNames[i], read.table(fname, sep= '\t', header= T))
}


pi_cgal <- list(pi_cgalchr2, pi_cgalchr7, pi_cgalchr10, pi_cgalchr18, pi_cgalchr21)
names(pi_cgal) <- chromnames
pi_mros <- list(pi_mroschr2, pi_mroschr7, pi_mroschr10, pi_mroschr18, pi_mroschr21)
names(pi_mros) <- chromnames
pi_pach <- list(pi_pachchr2, pi_pachchr7, pi_pachchr10, pi_pachchr18, pi_pachchr21)
names(pi_pach) <- chromnames

popls <- list(pi_cgal, pi_mros, pi_pach)
names(popls) <- c('cgal', 'mros', 'pach')

# throw PIs in bins and output stats

dat <- pi_cgal$chrom2
binstats <- stats.bin(dat$POS,dat$PI,N=100)

matplot(binstats$centers, t(binstats$stats[ c("mean", "median","Q1", "Q3"),]), type="l",lty=c(1,2,2,2), col=c('red','blue','green','purple'), ylab="Pi diversity")

# mean = red median = blue Q1 = green Q3 = purple

mean_pis <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, c('cgal', 'mros', 'pach')))
sd_pis <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, c('cgal', 'mros', 'pach')))
for(i in 1:length(popls)){
	curSpec <- popls[[i]]
	for(j in 1:5){
		mean_pis[j, i]  <- mean(popls[[i]][[j]]$PI, na.rm= T)
		sd_pis[j, i] <- sd(popls[[i]][[j]]$PI, na.rm= T)
	}
}


barPi <- barplot(mean_pis, beside= T, col= brewer.pal(5, 'Set3'), ylim= c(0, 0.5), xlab= 'species*chroms', ylab= 'Nucleotide diversity')
legend('topright', legend= chromnames, fill= brewer.pal(5, 'Set3'))
error.bar(barPi, mean_pis, sd_pis) #1.96*cs$dpSD/10)




#########################################################

##    fst   (vcftools)  windows


#lsfst <- system("ls pgStats/fst/singleLoc/", intern= T)
#fls.fst <- lsfst[grep("weir.fst", lsfst)]#[c(1:14, 16)]
#
#
#fstNames <- unlist(strsplit(fls.fst, '.windowed.weir.fst'))
#for(i in 1:length(fls.fst)){
#    #oname = paste(estname[i], "est", sep= "")
#	fname <- paste('pgStats/fst/singleLoc/', fls.fst[i], sep= '')
#	print(fname)
#	obj <- read.table(fname, sep= '\t', header= T)
#	obj$WEIGHTED_FST[which(obj$WEIGHTED_FST < 0)] <- NA
#    assign(fstNames[i], obj)	#read.table(fname, sep= '\t', header= T))
#}
#
#
#cgal_mros <- list(cgal_mrosChr2, cgal_mrosChr7, cgal_mrosChr10, cgal_mrosChr18, cgal_mrosChr21) 
#names(cgal_mros) <- chromnames
#pach_mros <- list(pach_mrosChr2, pach_mrosChr7, pach_mrosChr10, pach_mrosChr18, pach_mrosChr21) 
#names(pach_mros) <- chromnames
#cgal_pach <- list(cgal_pachChr2, cgal_pachChr7, cgal_pachChr10, cgal_pachChr18, cgal_pachChr21) 
#names(cgal_pach) <- chromnames
#
#fst <- list(cgal_mros, pach_mros, cgal_pach)
#names(fst) <- c('cgal_mros', 'pach_mros', 'cgal_pach')

## plot fst chrom per spec


plotFstPerSpec(fst, 2)



###

# matrix of mean Fst (and sd)   NOT DONE YET!
mean_fst <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, names(fst)))
sd_fst <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, names(fst)))
for(i in 1:length(fst)){
	currC <- fst[[i]]
	for(j in 1:5){
		mean_fst[j, i]  <- mean(fst[[i]][[j]]$WEIR_AND_COCKERHAM_FST, na.rm= T)
		sd_fst[j, i] <- sd(fst[[i]][[j]]$WEIR_AND_COCKERHAM_FST, na.rm= T)
	}
}

barFst <- barplot(mean_fst, beside= T, col= brewer.pal(5, 'Set3'), ylim= c(0, 1), xlab= 'species*chroms', ylab= expression(F[ST]), cex= 1.3, cex.axis= 1.4)
legend('topright', legend= chromnames, fill= brewer.pal(5, 'Set3'))
error.bar(barFst, mean_fst, sd_fst) #1.96*cs$dpSD/10)



# now a full table for SI table 2 
fstM <- matrix(nrow=5, ncol= 6, dimnames= list(rownames(mean_fst), paste0(c("cgmrMean", "cgmrSD", "pamrMean", "pamrSD", "cgpaMean", "cgpaSD"), "Fst")))
for(i in 1:3){
	j <- i*2
	fstM[,j-1] <- mean_fst[,i]
	fstM[,j] <- sd_fst[,i]
}


####     per locus fst 


lsfst <- system("ls pgStats/fst/singleLoc/", intern= T)
fls.fst <- lsfst[grep("weir.fst", lsfst)]#[c(1:14, 16)]


fstNames <- unlist(strsplit(fls.fst, '.weir.fst'))
for(i in 1:length(fls.fst)){
    #oname = paste(estname[i], "est", sep= "")
	fname <- paste('pgStats/fst/singleLoc/', fls.fst[i], sep= '')
	print(fname)
	obj <- read.table(fname, sep= '\t', header= T)
	obj$WEIR_AND_COCKERHAM_FST[which(obj$WEIR_AND_COCKERHAM_FST < 0)] <- 0
    assign(fstNames[i], obj)	#read.table(fname, sep= '\t', header= T))
}


cgal_mros <- list(cgal_mrosChr2, cgal_mrosChr7, cgal_mrosChr10, cgal_mrosChr18, cgal_mrosChr21) 
names(cgal_mros) <- chromnames
pach_mros <- list(pach_mrosChr2, pach_mrosChr7, pach_mrosChr10, pach_mrosChr18, pach_mrosChr21) 
names(pach_mros) <- chromnames
cgal_pach <- list(cgal_pachChr2, cgal_pachChr7, cgal_pachChr10, cgal_pachChr18, cgal_pachChr21) 
names(cgal_pach) <- chromnames

fst <- list(cgal_mros, pach_mros, cgal_pach)
names(fst) <- c('cgal_mros', 'pach_mros', 'cgal_pach')

##

#uQuant <- function(dat, thresh= 2){
#	vec <- dat$WEIR_AND_COCKERHAM_FST
#	#dat <- tajDat$TajimaD
#	lowerq = quantile(vec, na.rm= T)[2]
#	upperq = quantile(vec, na.rm= T)[4]
#	iqr = upperq - lowerq #Or use IQR(data)
#	#Compute the bounds for an outlier, given thresh (1.5, 3,...)
#	upper = (iqr * thresh) + upperq
#	lower = lowerq - (iqr * thresh)
#	#return(dat[which(vec < lower | vec > upper),])		# or only return indices? (which())
#	return(which(vec < lower | vec > upper))	
#}

#n <- 0.05
#subset(data, V2 > quantile(V2, prob = 1 - n))



fstChr2 <- list(cgal_mrosChr2, pach_mrosChr2, cgal_pachChr2)
names(fstChr2) <- c('cgal_mros', 'pach_mros', 'cgal_pach')
fstChr7 <- list(cgal_mrosChr7, pach_mrosChr7, cgal_pachChr7)
names(fstChr7) <- c('cgal_mros', 'pach_mros', 'cgal_pach')
fstChr10 <- list(cgal_mrosChr10, pach_mrosChr10, cgal_pachChr10)
names(fstChr10) <- c('cgal_mros', 'pach_mros', 'cgal_pach')
fstChr18 <- list(cgal_mrosChr18, pach_mrosChr18, cgal_pachChr18)
names(fstChr18) <- c('cgal_mros', 'pach_mros', 'cgal_pach')
fstChr21 <- list(cgal_mrosChr21, pach_mrosChr21, cgal_pachChr21)
names(fstChr21) <- c('cgal_mros', 'pach_mros', 'cgal_pach')

fstls <- list(fstChr2, fstChr7, fstChr10, fstChr18, fstChr21) 


# fst outliers 

fstQls <- list()
for(i in 1:length(fstls)){
	fstQls[[i]] <- getFquant(fstls[[i]], 0.99)
}



# get the number of outliers for each specBYspec   (SI table 2)
fstOutM <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, paste0(names(fstQls[[1]]), "Fst")))
for(i in 1:3){
	fstOutM[,i] <- unlist(lapply(lapply(fstQls, '[[', 1), length))
}



# get the xth quantile (< outlier)


pQls <- list()
for(i in 1:length(pDifls)){
	pQls[[i]] <- getPquant(pDifls[[i]], 0.99)
}


pOutM <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, paste0(names(pQls[[1]]), "afDiffs")))
for(i in 1:3){
	pOutM[,i] <- unlist(lapply(lapply(pQls, '[[', 1), length))
}

p_fst_outN <- cbind(pOutM[,1], fstOutM[,1])
colnames(p_fst_outN) <- c('afDiffs', 'Fst')




#
getQld2 <- function(x, g, pops, n= 0.05){
	spec1 <- pops[1]
	spec2 <- pops[2]
	specs <- c(rep('cgal', 10), rep('mros', 10), rep('pach', 10))
	indXs <- c(1:10, 11:20, 21:30)
	specX1 <- which(specs == spec1)
	specX2 <- which(specs == spec2)
	selX <- which(x$WEIR_AND_COCKERHAM_FST > quantile(x$WEIR_AND_COCKERHAM_FST, prob= 1 - n, na.rm= T))
	neutX <- which(x$WEIR_AND_COCKERHAM_FST <= quantile(x$WEIR_AND_COCKERHAM_FST, prob= 1 - n, na.rm= T))
	ps <- g[selX,]
	#ps <- subset(g, WEIR_AND_COCKERHAM_FST > quantile(WEIR_AND_COCKERHAM_FST, prob = 1 - n))
	#pn <- subset(g, WEIR_AND_COCKERHAM_FST < quantile(WEIR_AND_COCKERHAM_FST, prob = 1 - n))
	invar1 <- which(apply(ps[,specX1], 1, sd) == 0)
	invar2 <- which(apply(ps[,specX2], 1, sd) == 0)
	LDmat1 <- cor(t(ps[-invar1,specX1]))^2
	LDmat2 <- cor(t(ps[-invar2,specX2]))^2
	#
	nps <- dim(LDmat1)[1]
	pn <- g[neutX,]
	#pn <- g[-uQuant(x, qtt),]
	iNvar1 <- which(apply(pn[,specX1], 1, sd) == 0)
	iNvar2 <- which(apply(pn[,specX2], 1, sd) == 0)
	ppn1 <- pn[-iNvar1, specX1]
	ppn2 <- pn[-iNvar2, specX2]
	snx <- sample(nrow(ppn1), nps, replace= F)
	pN1 <- ppn1[snx,]
	pN2 <- ppn2[snx,]
	LDmatN1 <- cor(t(pN1))^2
	LDmatN2 <- cor(t(pN2))^2
	LDlist <- list(list(LDmat1, LDmatN1), list(LDmat2, LDmatN2))
	names(LDlist) <- pops
	return(LDlist)
}

#ld21SN1 <- getQld(fst[[1]][[5]], chr21, c('cgal', 'mros'), qtt= 1.5)
#ld18SN2 <- getQld(fst[[2]][[4]], chr18, c('pach', 'mros'), qtt= 1.5)
#ldSN3 <- getQld2(fst[[3]][[1]], chr2, c('cgal', 'pach'), qtt= 1.5)

########

getQldPall <- function(pDiff, g, pops, qtt= 0.95){
	spec1 <- pops[1]
	spec2 <- pops[2]
	specs <- c(rep('cgal', 10), rep('mros', 10), rep('pach', 10))
	indXs <- c(1:10, 11:20, 21:30)
	specX1 <- which(specs == spec1)
	specX2 <- which(specs == spec2)
	selX <- getPquant(pDiff, qtt)#[[1]]
	ppairs <- names(selX)
	selXX <- selX[[which(ppairs == paste(pops[1], pops[2], sep= '_'))]]
	
	pps <- getPPS(gChr[selXX,], specX1, specX2)
	pps1 <- pps[[1]]
	pps2 <- pps[[2]]
	
	LDmat1 <- cor(t(pps1))^2
	LDmat2 <- cor(t(pps2))^2
	#
	nps <- dim(LDmat1)[1]
	pn <- g[-selXX,]
	iNvar1 <- which(apply(pn[,specX1], 1, sd) == 0)
	iNvar2 <- which(apply(pn[,specX2], 1, sd) == 0)
	ppn1 <- pn[-iNvar1, specX1]
	ppn2 <- pn[-iNvar2, specX2]
	snx <- sample(nrow(ppn1), nps, replace= F)
	pN1 <- ppn1[snx,]
	pN2 <- ppn2[snx,]
	LDmatN1 <- cor(t(pN1))^2
	LDmatN2 <- cor(t(pN2))^2
	#
	# s - n 
	mSN1 <- matrix(nrow= nps, ncol= nps)
	mSN2 <- matrix(nrow= nps, ncol= nps)
	for(i in 1:nps){
		selloc1 <- pps1[i,]
		selloc2 <- pps2[i,]
		mSN1[i,] <- cor(t(pN1), as.numeric(selloc1))^2
		mSN2[i,] <- cor(t(pN2), as.numeric(selloc2))^2
	}
	LDlist <- list(list(LDmat1, LDmatN1, mSN1), list(LDmat2, LDmatN2, mSN2))
	names(LDlist) <- pops
	return(LDlist)
}
# (NOTE: that the neutral dfs in the above function do not check for is.integer0() !) 



sn2 <- getQldPall(pDiff, chr2, c('cgal', 'mros'), qtt= 0.99)




##############################################


sn2 <- getQldPallChr(chrls= chrl, pQls= pQls, tChr= 2, pops= c('cgal', 'mros'))


#
#pdf("test.pdf") 
#plot(rnorm(10),rnorm(10)) 
#hist(sn2[[1]][[1]][upper.tri(sn2[[1]][[1]])])
#dev.off() 

#chr2cg_mr <- getQldP(pDiff2, chr2, c('cgal', 'mros'), qtt= 0.98)
#chr2pa_mr <- getQldP(pDiff2, chr2, c('pach', 'mros'), qtt= 0.98)
#chr2cg_pa <- getQldP(pDiff2, chr2, c('cgal', 'pach'), qtt= 0.98)

dat <- sn2		
pops <- names(dat)	
pdf('test1.pdf')
par(mfrow= c(4,3))		#, mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))		#oma=c(4.5,4.5,4.5,4.5), 
for(i in 1:2){
	#pops <- names(dat)
	hist(dat[[i]][[1]], main= paste(pops[i], "selected sites"), xlab= expression(paste('r'^'2')))
	hist(dat[[i]][[2]], main= paste(pops[i], "neutral sites"), xlab= expression(paste('r'^'2')))
	hist(dat[[i]][[3]], main= paste(pops[i], "selected vs. neutral sites"), xlab= expression(paste('r'^'2')))
	hist(dat[[i]][[4]], main= paste(pops[i], "selected sites w/ 'other' s"), xlab= expression(paste('r'^'2')))
	hist(dat[[i]][[5]], main= paste(pops[i], "neutral sites w/ 'other' n"), xlab= expression(paste('r'^'2')))
	hist(dat[[i]][[6]], main= paste(pops[i], "selected w/ 'other' n sites"), xlab= expression(paste('r'^'2')))
}
dev.off()



#x <- matrix( runif(12) , nrow = 3 )
#y <- runif(3)
#cor(x,y)

#apply(x, 2 , cor , y = as.numeric(y))
x <- lapply(chrl[-1], "[", sn2[[1]])

mapply(rep, 1:4, 4:1)
     
mapply(rep, times = 1:4, x = 4:1)

mapply(rep, times = 1:4, MoreArgs = list(x = 42))

mapply(function(x, y) seq_len(x) + y, 
	c(a =  1, b = 2, c = 3),  # names from first
    c(A = 10, B = 0, C = -10))



########
#fstLDlist <- list()
#ld  <- list()
#for(i in 1:length(fst)){
#	dat <- fst[[i]]
#	name <- names(fst)[i]
#	for(j in 1:length(dat)){
#		chrDat <- dat[[j]]
#		ld[[j]] <- getQld(chrDat, chrl[[j]], strsplit(name, '_'), qtt= 1.5)
#		#print(dim(chrDat))
#		#print(dim(chrl[[j]]))
#	}
#	fstLDlist[[i]] <- ld
#}



# check for outlier numbers  (also implement for afDiffs!!) 
for(i in 1:length(fst)){
	for(j in 1:length(fst[[i]])){
		#print(length(uQuant(fst[[i]][[j]], 2.0)))
		print(length(which(fst[[i]][[j]]$WEIR_AND_COCKERHAM_FST > quantile(fst[[i]][[j]]$WEIR_AND_COCKERHAM_FST, prob= 1 - 0.05, na.rm= T))))
	}
}

# now re-write getQldP to calc r^2 for s-s, n-n, s-n, 1) within chromosomes (done) and 2) between chromosomes!


par(mfrow= c(2,2))
hist(ldSN[[1]][upper.tri(ldSN[[1]])])
hist(ldSN[[2]][upper.tri(ldSN[[2]])])




### LD 

# table with first 10 SNPs
cor(chr2[1:10,])^2
## remove invariant and grab subset
chr2sub <- chr2[1001:2000,]
invar <- which(apply(chr2sub, 1, sd) == 0)
#LDmat <- (cor(t(chr2sub[-invar,]))^2)
LDmat <- (cor(t(chr2sub))^2)
hist(as.numeric(LDmat[upper.tri(LDmat)]))
image(LDmat)



nsnp<-dim(LDmat)[1]
Dmat<-LDmat
for(i in 1:nsnp){
	for(j in 1:nsnp){
		Dmat[i,j]<-abs(i-j)
	}
}

plot(Dmat[upper.tri(Dmat)],LDmat[upper.tri(LDmat)])
cor(Dmat[upper.tri(Dmat)],LDmat[upper.tri(LDmat)])




#par(mfrow= c(3,1), oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))



#########################################################

## Tajima's D    (vcftools)


lstajD <- system("ls pgStats/tajD/", intern= T)
fls.tajD <- lstajD[grep(".Tajima.D", lstajD)]				# change here (w1k= .Tajima.D; w10k= w10k.Tajima.D, w100bp.TajimaD)

# fls.tajD <- fls.tajD[seq(1,45,3)]		# for 1k (*Chr*.TajimaD)

#tajDNames <- unlist(strsplit(fls.tajD, '.windowed.weir.fst'))
for(i in 1:length(fls.tajD)){
    #oname = paste(estname[i], "est", sep= "")
	fname <- paste('pgStats/tajD/', fls.tajD[i], sep= '')
	print(fname)
    assign(fls.tajD[i], read.table(fname, sep= '\t', header= T))
}
#### w100bp

tajDcgal <- list(cgalChr2w100bp.Tajima.D, cgalChr7w100bp.Tajima.D, cgalChr10w100bp.Tajima.D, cgalChr18w100bp.Tajima.D, cgalChr21w100bp.Tajima.D)
names(tajDcgal) <- chromnames
tajDmros <- list(mrosChr2w100bp.Tajima.D, mrosChr7w100bp.Tajima.D, mrosChr10w100bp.Tajima.D, mrosChr18w100bp.Tajima.D, mrosChr21w100bp.Tajima.D)
names(tajDmros) <- chromnames
tajDpach <- list(pachChr2w100bp.Tajima.D, pachChr7w100bp.Tajima.D, pachChr10w100bp.Tajima.D, pachChr18w100bp.Tajima.D, pachChr21w100bp.Tajima.D)
names(tajDpach) <- chromnames
#
tajDs <- list(tajDcgal, tajDmros, tajDpach)
names(tajDs) <- c('cgal', 'mros', 'pach')



####  w1k 
tajDcgal <- list(cgalChr2.Tajima.D, cgalChr7.Tajima.D, cgalChr10.Tajima.D, cgalChr18.Tajima.D, cgalChr21.Tajima.D)
names(tajDcgal) <- chromnames
tajDmros <- list(mrosChr2.Tajima.D, mrosChr7.Tajima.D, mrosChr10.Tajima.D, mrosChr18.Tajima.D, mrosChr21.Tajima.D)
names(tajDmros) <- chromnames
tajDpach <- list(pachChr2.Tajima.D, pachChr7.Tajima.D, pachChr10.Tajima.D, pachChr18.Tajima.D, pachChr21.Tajima.D)
names(tajDpach) <- chromnames
#
tajDs <- list(tajDcgal, tajDmros, tajDpach)
names(tajDs) <- c('cgal', 'mros', 'pach')

####   w10k
tajDcgal <- list(cgalChr2w10k.Tajima.D, cgalChr7w10k.Tajima.D, cgalChr10w10k.Tajima.D, cgalChr18w10k.Tajima.D, cgalChr21w10k.Tajima.D)
names(tajDcgal) <- chromnames
tajDmros <- list(mrosChr2w10k.Tajima.D, mrosChr7w10k.Tajima.D, mrosChr10w10k.Tajima.D, mrosChr18w10k.Tajima.D, mrosChr21w10k.Tajima.D)
names(tajDmros) <- chromnames
tajDpach <- list(pachChr2w10k.Tajima.D, pachChr7w10k.Tajima.D, pachChr10w10k.Tajima.D, pachChr18w10k.Tajima.D, pachChr21w10k.Tajima.D)
names(tajDpach) <- chromnames
#
tajDs <- list(tajDcgal, tajDmros, tajDpach)
names(tajDs) <- c('cgal', 'mros', 'pach')


####
# matrix of mean Tajima's D (and sd)     (NOT SURE HOW USEFUL THIS IS!?)
mean_tajD <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, c('cgal', 'mros', 'pach')))
sd_tajD <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, c('cgal', 'mros', 'pach')))
for(i in 1:length(tajDs)){
	curSpec <- tajDs[[i]]
	for(j in 1:5){
		mean_tajD[j, i]  <- mean(tajDs[[i]][[j]]$TajimaD, na.rm= T)
		sd_tajD[j, i] <- sd(tajDs[[i]][[j]]$TajimaD, na.rm= T)
	}
}

barTajD <- barplot(mean_tajD, beside= T, col= brewer.pal(5, 'Set3'), ylim= c(0, 1.8), xlab= 'species*chroms', ylab= expression(F[ST]), cex= 1.3, cex.axis= 1.4)
legend('topleft', legend= chromnames, fill= brewer.pal(5, 'Set3'))
error.bar(barTajD, mean_tajD, sd_tajD) #1.96*cs$dpSD/10)

# or: 

par(mfrow= c(5,1), oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
plotVioChr(tajDs, 2)


###
par(mfrow= c(3,1), oma=c(4.5,4.5,4.5,4.5), mar=c(4.0,4.0,4.0,4.0) + 0.2, mgp= c(3,1,0))
plotVioSpec(tajDs, 'cgal')

# now a full table for SI table 2   (with fstM) 
tajDM <- matrix(nrow=5, ncol= 6, dimnames= list(rownames(mean_tajD), paste0(c("cgmrMean", "cgmrSD", "pamrMean", "pamrSD", "cgpaMean", "cgpaSD"), "TajD")))
for(i in 1:3){
	j <- i*2
	tajDM[,j-1] <- mean_tajD[,i]
	tajDM[,j] <- sd_tajD[,i]
}





#ggplot(, aes(x=dose, y=len)) + geom_violin(trim=FALSE)
#	
#	
#	par(mfrow=c(2,1))
#       mu<-2
#       si<-0.6
#       bimodal<-c(rnorm(1000,-mu,si),rnorm(1000,mu,si)) 
#       uniform<-runif(2000,-4,4)
#       normal<-rnorm(2000,0,3)
#       vioplot(bimodal,uniform,normal)
#       boxplot(bimodal,uniform,normal)
#       
#       # add to an existing plot
#       x <- rnorm(100)
#       y <- rnorm(100)
#       plot(x, y, xlim=c(-5,5), ylim=c(-5,5))
#       vioplot(x, col="tomato", horizontal=TRUE, at=-4, add=TRUE,lty=2, rectCol="gray")
#       vioplot(y, col="cyan", horizontal=FALSE, at=-4, add=TRUE,lty=2)


# tajD for each chromosome 

plotTajimaD(tajDs, spec= 'mros')


#specs <- c(rep('c. galanthus', 10), rep('m. rosina', 10), rep('pachinus', 10))

#######################################

# get the number of Tajima's D outliers (+/-) for each chromosome

dat2 <- tajDs$cgal$chrom2$TajimaD 
dat7 <- tajDs$cgal$chrom7$TajimaD
dat10 <- tajDs$cgal$chrom10$TajimaD
dat18 <- tajDs$cgal$chrom18$TajimaD
dat21 <- tajDs$cgal$chrom21$TajimaD


#cgalChr21w4bp <- read.table("pgStats/tajD/cgalChr21w4bp.Tajima.D", sep = '\t', header= T)

par(mfrow=c(3,3))
for(i in 1:3){
	for(j in 1:3){
		

# get all N_SNPS and TajimaD values for cgal
nSNPs <- unlist(lapply(tajDs$cgal, '[[', 3))
tajDvals <- unlist(lapply(tajDs$cgal, '[[', 4))


plot(nSNPs, tajDvals, xlab= 'Number of SNPs in 10kbp window', ylab= "Tajima's D")



#########################################################

## snp density 

lsSNPd <- system("ls pgStats/snpDensity/", intern= T)
fls.SNPd <- lsSNPd[grep(".snpden", lsSNPd)]#[c(1:14, 16)]
#SNPdNames <- unlist(strsplit(.fst, '.windowed.weir.fst'))
for(i in 1:length(fls.SNPd)){
    #oname = paste(estname[i], "est", sep= "")
	fname <- paste('pgStats/snpDensity/', fls.SNPd[i], sep= '')
	print(fname)
    assign(fls.SNPd[i], read.table(fname, sep= '\t', header= T))
}


vioplot(chr2.snpden$SNP_COUNT, chr7.snpden$SNP_COUNT, chr10.snpden$SNP_COUNT, chr18.snpden$SNP_COUNT, chr21.snpden$SNP_COUNT) 

# hist
par(mfrow= c(1,5))
hist(chr2.snpden$SNP_COUNT)
hist(chr7.snpden$SNP_COUNT)
hist(chr10.snpden$SNP_COUNT)
hist(chr18.snpden$SNP_COUNT)
hist(chr21.snpden$SNP_COUNT)

# plot 
par(mfrow= c(5,1), mar=c(1.5,2,1,1) + 0.1, oma= c(5,0,0,0), mgp= c(0,1,0))
plot(chr2.snpden$SNP_COUNT, type= 'l')
plot(chr7.snpden$SNP_COUNT, type= 'l')
plot(chr10.snpden$SNP_COUNT, type= 'l')
plot(chr18.snpden$SNP_COUNT, type= 'l')
plot(chr21.snpden$SNP_COUNT, type= 'l')



#########################################################

## div Stats    (from genomics.py; see github.com/simonhmartin/genomics_general)



lsDiv <- system("ls pgStats/div_genomicsPy/", intern= T)

# window == 1000
flsW1000m10 <- lsDiv[grep("w1000m10", lsDiv)]#[c(1:14, 16)]
flsDiv <- flsW1000m10
oName <- unlist(strsplit(flsDiv, '_'))[seq(1, 10, 2)]
oName <- paste0(unlist(strsplit(oName, '30'))[seq(2, 10, 2)], '_w1k')
for(i in 1:length(flsDiv)){
    #oname = paste(estname[i], "est", sep= "")
	fname <- paste('pgStats/div_genomicsPy/', flsDiv[i], sep= '')
	print(fname)
    assign(oName[i], read.csv(fname, header= T))
}

divW1k <- list(chr2div_w1k, chr7div_w1k, chr10div_w1k, chr18div_w1k, chr21div_w1k) 
names(divW1k) <- chromnames

####

# window == 10000
flsW10000m100 <- lsDiv[grep("w10000m100", lsDiv)]#[c(1:14, 16)]
flsDiv <- flsW10000m100
oName <- unlist(strsplit(flsDiv, '_'))[seq(1, 10, 2)]
oName <- paste0(unlist(strsplit(oName, '30'))[seq(2, 10, 2)], '_w10k')
for(i in 1:length(flsDiv)){
    #oname = paste(estname[i], "est", sep= "")
	fname <- paste('pgStats/div_genomicsPy/', flsDiv[i], sep= '')
	print(fname)
    assign(oName[i], read.csv(fname, header= T))
}

divW10k <- list(chr2div_w10k, chr7div_w10k, chr10div_w10k, chr18div_w10k, chr21div_w10k) 
names(divW10k) <- chromnames


########################## plots  (div_genomicsPy)

## Dxy


plotDxy(divW10k, 2)



## Fst 


plotFst(divW10k, 2)



## pi


plotPi(divW10k, 2)


## Tajima's D 


####

# get mean values for Fst, Dxy, and pi for each chrom   (still genomics.py)

# mean Fst (and sd)     divW1k
divObj <- divW1k


#getEstStat <- function(divObj, stat= 'Fst'){
#	cols <- names(divObj[[1]])

fst_cgal_mros <- lapply(divObj, '[[', 12)
meanSd1 <- list(unlist(lapply(fst_cgal_mros, mean, na.rm= T)), unlist(lapply(fst_cgal_mros, sd, na.rm= T)))

fst_pach_mros <- lapply(divObj, '[[', 14)
meanSd2 <- list(unlist(lapply(fst_pach_mros, mean, na.rm= T)), unlist(lapply(fst_pach_mros, sd, na.rm= T)))

fst_cgal_pach <- lapply(divObj, '[[', 13)
meanSd3 <- list(unlist(lapply(fst_cgal_pach, mean, na.rm= T)), unlist(lapply(fst_cgal_pach, sd, na.rm= T)))

fstL <- list(meanSd1, meanSd2, meanSd3) 

fstM <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, c('cgal_mros', 'pach_mros', 'cgal_pach')))
fstSD <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, c('cgal_mros', 'pach_mros', 'cgal_pach')))

for(i in 1:3){
	#for(j in 1:5){
		fstM[,i] <- fstL[[i]][[1]]
		fstSD[,i] <- fstL[[i]][[2]]
}

barFst <- barplot(fstM, beside= T, col= brewer.pal(5, 'Set3'), ylim= c(0, 1), xlab= 'species*chroms', ylab= expression(F[ST]), cex.names= 1.4, cex.axis= 1.4)
legend('topright', legend= chromnames, fill= brewer.pal(5, 'Set3'))
error.bar(barFst, fstM, fstSD) #1.96*cs$dpSD/10)



## Dxy

divObj <- divW1k

dxy_cgal_mros <- lapply(divObj, '[[', 9)
meanSd1 <- list(unlist(lapply(dxy_cgal_mros, mean, na.rm= T)), unlist(lapply(dxy_cgal_mros, sd, na.rm= T)))

dxy_pach_mros <- lapply(divObj, '[[', 11)
meanSd2 <- list(unlist(lapply(dxy_pach_mros, mean, na.rm= T)), unlist(lapply(dxy_pach_mros, sd, na.rm= T)))

dxy_cgal_pach <- lapply(divObj, '[[', 10)
meanSd3 <- list(unlist(lapply(dxy_cgal_pach, mean, na.rm= T)), unlist(lapply(dxy_cgal_pach, sd, na.rm= T)))

dxyL <- list(meanSd1, meanSd2, meanSd3) 

dxyM <- matrix(nrow= 5, ncol= 6, dimnames= list(chromnames, paste0(c("cgmrMean", "cgmrSD", "pamrMean", "pamrSD", "cgpaMean", "cgpaSD"), '_Dxy')))
#dxySD <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, c('cgal_mros', 'pach_mros', 'cgal_pach')))

for(i in 1:3){
	#for(j in 1:5){
	j <- i*2
		dxyM[,j-1] <- dxyL[[i]][[1]]
		dxyM[,j] <- dxyL[[i]][[2]]
}


###
#for(i in 1:3){
#	j <- i*2
#	tajDM[,j-1] <- mean_tajD[,i]
#	tajDM[,j] <- sd_tajD[,i]
#}
#
#
#barDxy <- barplot(dxyM, beside= T, col= brewer.pal(5, 'Set3'), ylim= c(0, 1), xlab= 'species*chroms', ylab= expression(D[XY]), cex.names= 1.4, cex.axis= 1.4)
#legend('topright', legend= chromnames, fill= brewer.pal(5, 'Set3'))
#error.bar(barDxy, dxyM, dxySD) #1.96*cs$dpSD/10)

## pi

divObj <- divW1k

pi_cgal <- lapply(divObj, '[[', 6)
meanSd1 <- list(unlist(lapply(pi_cgal, mean, na.rm= T)), unlist(lapply(pi_cgal, sd, na.rm= T)))

pi_mros <- lapply(divObj, '[[', 7)
meanSd2 <- list(unlist(lapply(pi_mros, mean, na.rm= T)), unlist(lapply(pi_mros, sd, na.rm= T)))

pi_pach <- lapply(divObj, '[[', 8)
meanSd3 <- list(unlist(lapply(pi_pach, mean, na.rm= T)), unlist(lapply(pi_pach, sd, na.rm= T)))

piL <- list(meanSd1, meanSd2, meanSd3) 

piM <- matrix(nrow= 5, ncol= 6, dimnames= list(chromnames, paste0(c("cgMean", "cgSD", "mrMean", "mrSD", "paMean", "paSD"), '_pi')))
#piSD <- matrix(nrow= 5, ncol= 3, dimnames= list(chromnames, c('cgal', 'mros', 'pach')))

for(i in 1:3){
	#for(j in 1:5){
	j <- i*2
	piM[,j-1] <- piL[[i]][[1]]
	piM[,j] <- piL[[i]][[2]]
}

barPi <- barplot(piM, beside= T, col= brewer.pal(5, 'Set3'), ylim= c(0, 1), xlab= 'species*chroms', ylab= expression(pi), cex.names= 1.4, cex.axis= 1.4)
legend('topright', legend= chromnames, fill= brewer.pal(5, 'Set3'))
error.bar(barPi, piM, piSD) #1.96*cs$dpSD/10)

#################################################################################################

#################################################################################################

# load genotypes
# gt <- read.csv('var/helico30f3maf10.txt', header= F, sep= ' ')

# 
#L <- dim(gt)[1]
#N <- length(unique(ids$species))
#
#colnames(gt) <- ids$sample
##p <- as.matrix(gt)
#
### population allele freqs., simple ML estimate
#p <- matrix(NA, nrow=L,ncol=N)		# NOTE: probably switch row & cols !!!
#for(i in 1:L){
#    p[i,]<-tapply(X=gt[i,],INDEX=ids[,2],mean)/2
#    }

#rownames(p)<-unique(ids[,2])


## Fst 
#pbar <- apply(p,1,mean)
#vp <- apply(p,1,var)
#fst <- mean(vp)/mean(pbar * (1 - pbar)) ## ratio of averages
#
### pairwise Fst for pop pairs
#lp <- p[c(3, 4),]
#pbar_lp <- apply(lp, 2, mean)
#vp_lp <- apply(lp, 2, var)
#fst_lp <- mean(vp_lp)/mean(pbar_lp * (1- pbar_lp))
#
### heterozygosity
#het<-2 * p * (1-p)
#mnH<-apply(het,1,mean)
#barplot(mnH,las=2)
#
#par(mfrow= c(1,2))
#hist(p, main= "reference allele")
#hist(1-p, main= "alternative allele")


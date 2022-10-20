### r/qtl to build a genetic map ###

### Data summary and filtering
library(qtl)
data(mapthis)
summary(mapthis)

# missing data
par(las=1, cex=0.8, mar=c(4.1, 4.1, 3, 1.1))
plotMissing(mapthis, main="Missing data")

# individuals / markers
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")




### Segregation Distortion
gt <- geno.table(mapthis)
bon.pval.cutoff <- 0.05/totmar(mapthis) ### Bonferroni corrected p-value cutoff
gt[gt$P.value < bon.pval.cutoff, ]
todrop <- rownames(gt[gt$P.value < bon.pval.cutoff, ])
mapthis <- drop.markers(mapthis, todrop) ### drop markers showing significant segregation distortion




### Plot of LOD scores versus estimated recombination fractions for all marker pairs
mapthis <- est.rf(mapthis)

rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
par(mar=c(4.1,4.1,0.6,0.6), las=1, cex=0.8)
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


### Linkage groups
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)


### Order markers
mapthis <- orderMarkers(mapthis, use.ripple=TRUE, window = 7, chr=5)
pull.map(mapthis, chr=5)


### switch alleles (specifically for heterozygous parents)
plotRF(mapthis, alternate.chrid=TRUE)

toswitch <- markernames(mapthis, chr=c(5, 7:11))
mapthis <- switchAlleles(mapthis, toswitch)
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
mapthis <- orderMarkers(mapthis, use.ripple=TRUE, window = 7, chr=5)


### Look for problem individuals
plot(countXO(mapthis), ylab="Number of crossovers")


### Estimate genotyping error rate, bugs contained
loglik <- err <- c(0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
for(i in seq(along = err)) {
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(mapthis, error.prob = err[i])
  loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)

### plot MLE
plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))


### Look for genotyping errors and zero out these suspicious genotypes (that is, make them missing (NA))
mapthis <- calc.errorlod(mapthis, error.prob=0.005)  ### determine error LOD
toperr <- top.errorlod(mapthis, cutoff=4)
print(toperr)

mapthis.clean <- mapthis
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  mapthis.clean$geno[[chr]]$data[mapthis$pheno$id==id, mar] <- NA
}


### Summary and plot map
summaryMap(mapthis)
plotRF(mapthis, alternate.chrid=TRUE)

par(las=1, mar=c(4.1,4.1,1.1,0.1), cex=0.3)
plotMap(mapthis, main="", show.marker.names=TRUE)

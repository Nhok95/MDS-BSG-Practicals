#
#
#

#install.packages("HardyWeinberg")
library(HardyWeinberg)

x <- c(23,48,29)
names(x) <- c("CC", "CG", "GG")

x

results.chi <- HWChisq(x)

results.chi.nocor <- HWChisq(x,cc=0)

results.exact <- HWExact(x)

y <- c(0,7,93)
names(y) <- c("CC", "CT", "TT")

y

results.chi <- HWChisq(y)

results.chi.nocor <- HWChisq(y,cc=0)

results.exact <- HWExact(y)

HWAlltests(x,include.permutation.test=TRUE)
HWAlltests(y,include.permutation.test=TRUE)

X <- rbind(x,y)
X
HWTernaryPlot(X)

x
genotypestring <- rep(c("CC","CG","GG"),x)
genotypestring

allelestring <- unlist(strsplit(genotypestring,""))
allelestring
table(allelestring)

x.perm <- sample(allelestring)
table(x.perm)

ind1 <- seq(1,length(x.perm),2)
ind2 <- seq(2,length(x.perm),2)

geno.perm <- paste(x.perm[ind1],x.perm[ind2],sep="")
table(geno.perm)
geno.perm[geno.perm=="GC"] <- "CG"

geno.counts <- as.vector(table(geno.perm))
names(geno.counts) <- names(table(geno.perm))
geno.counts

#
# a permutation procedure
#

nsimul <- 20000
chistats <- numeric(nsimul)
for(i in 1:nsimul) {
  x.perm <- sample(allelestring)
  ind1 <- seq(1,length(x.perm),2)
  ind2 <- seq(2,length(x.perm),2)
  geno.perm <- paste(x.perm[ind1],x.perm[ind2],sep="")
  geno.perm[geno.perm=="GC"] <- "CG"
  geno.counts <- as.vector(table(geno.perm))
  names(geno.counts) <- names(table(geno.perm))
  chistats[i] <- HWChisq(geno.counts,cc=0,verbose=FALSE)$chisq
}

chisq.observed.sample <- HWChisq(x,cc=0)$chisq 

permutation.pvalue <- sum(chistats >= chisq.observed.sample)/nsimul
permutation.pvalue

#
# simulation
#

set.seed(123)
X <- HWData(100,nm=100)
head(X)

out <- HWTernaryPlot(X,region=1)

out <- HWTernaryPlot(X,region=7,verbose=FALSE)

statistics <- HWChisqStats(X)
hist(statistics)
plot(density(statistics[!is.na(statistics)]))

X <- HWData(100,nm=10000)
statistics <- HWChisqStats(X)
hist(statistics,breaks=30)
plot(density(statistics[!is.na(statistics)]))



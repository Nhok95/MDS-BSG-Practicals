
# STR dataset

rm(list=ls())

library(HardyWeinberg)

# Load data
data(NistSTRs)

# the rownames of the object consist of identifiers for each individual
# successive columns represent the two alleles of an individual for each STR

X <- NistSTRs

# Question 2
n <- nrow(X) # number of individuals
p <- ncol(X)/2 # number of STRs
n
p

# There are 361 individuals and 29 STRs

# Question 3

# Function that determines the number of alleles for a STR.
n.alleles <- function(X, str.index) {
  allele.1 <- as.list(X[,str.index])
  allele.2 <- as.list(X[,(str.index+1)])
  return(length(table(unlist(c(allele.1, allele.2))))) # number of alleles
}

n.alleles.per.str.list <- list()
str.index <- 1
for (str.num in 1:p) {
  n.alleles.per.str.list  <- append(n.alleles.per.str.list, n.alleles(X, str.index))
  str.index <- str.index + 2
}
n.alleles.per.str <- unlist(n.alleles.per.str.list)

# Basic descriptive statistics of the number of alleles
mean(n.alleles.per.str)
sd(n.alleles.per.str)
median(n.alleles.per.str)
max(n.alleles.per.str)
min(n.alleles.per.str)

# Question 4
barplot(table(n.alleles.per.str), xlab="Number of alleles", ylab="Number of STRs")
# The most common number of alleles for an STR is 8

# Question 5: Compute the expected heterozygosity for each STR
exp.heter <- function(X, str.index) {
  allele.1 <- as.list(X[,str.index])
  allele.2 <- as.list(X[,(str.index+1)])
  t <- table(unlist(c(allele.1, allele.2)))
  sum.t <- sum(unname(t)) # we sum the counts
  exp.heter <- round(1 - sum(sapply(unname(t), function(x) (x / sum.t)^2 )), 3)
  return(exp.heter) # expected heterozygosity formula
}

exp.heter.per.str.list <- list()
str.index <- 1
for (str.num in 1:p) {
  exp.heter.per.str.list  <- append(exp.heter.per.str.list, exp.heter(X, str.index))
  str.index <- str.index + 2
}
exp.heter.per.str <- unlist(exp.heter.per.str.list)

hist(exp.heter.per.str, xlab="Expected heterozygosity", main="Histogram of the expected heterozygosity")

round(mean(exp.heter.per.str), 3) # average expected heterozygosity over all STRs

# Question 6: Compute the observed heterozygosity for each STR

# genotype frequencies need to be calculated!
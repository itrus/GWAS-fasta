# ivan.trus@gmail.com 07/03/2015

# Initialization
start.time <- Sys.time()
cat("\014")
library("seqinr")

# Reading samples (plus.fas) and controls (minus.fas).
# All sequences should have the same length
# (They should be aligned all together before sorting them to different files.)
cat("Reading sequences\n")
seq.minus <- read.alignment(file="minus.fas", format="fasta")
seq.plus <- read.alignment(file="plus.fas", format="fasta")

# Reading intervals
seq.start <- 1 # Start position
#
seq.length <- nchar(seq.minus$seq[[1]]) # End position
#seq.length <- 100

# Calculating the number of the most frequent aminoacid at each position of
# controls
cat("Calculating moda values in negatives\n")
moda.minus <- rep(0, seq.length)
count.moda.minus <- rep(0, seq.length)
allele <- rep("", seq.minus$nb)
for (k in seq.start:(seq.length + seq.start - 1)) {
  for (i in 1:seq.minus$nb) {
    allele[i] <- substr(seq.minus$seq[[i]], k, k)
  }
  # Possible amino-acids are given
  j <- count(allele, 1, alphabet=s2c("arndceqghilkmfpstwyv-*?bzjx"))
  moda.minus[k - seq.start + 1] <- which.max(j)
  count.moda.minus[k - seq.start + 1] <- j[[moda.minus[k - seq.start + 1]]]
}

# Calculation of frequncy of finding the same aminoacid in positive samples
cat("Calculating moda values in positives\n")
count.moda.plus <- rep(0, seq.length)
allele <- rep("", seq.plus$nb)
for (k in seq.start:(seq.length + seq.start - 1)) {
  for (i in 1:seq.plus$nb) {
    allele[i] <- substr(seq.plus$seq[[i]], k, k)
  }
  j <- count(allele, 1, alphabet=s2c("arndceqghilkmfpstwyv-*?bzjx"))
  count.moda.plus[k - seq.start + 1] <- j[moda.minus[k - seq.start + 1]]
}

# Calculation of P level (prob) for each position with Fisher's or Chi test
cat("Calculating P values\n")
prob <- rep(1, seq.length)
for (k in seq.start:(seq.length + seq.start - 1)) {
  j <- rbind(c(count.moda.plus[k - seq.start + 1],
              seq.plus$nb - count.moda.plus[k - seq.start + 1]),
             c(count.moda.minus[k - seq.start + 1],
              seq.minus$nb - count.moda.minus[k - seq.start + 1]))
  # Chi squared test (chisq.test) can be used instead of Fisher's test (fisher.test)
  # if all tested values (seq.minus$nb, seq.plus$nb, count.moda.plus,
  # count.moda.minus) >=5
  prob[k - seq.start + 1] <- chisq.test(as.table(j))$p.value
}

# Plotting of a Manhattan plot with three reference lines
# (P=0.001, P=0.01, P=0.05)
cat("Plotting chart\n")
plot(-log10(prob),
     col  = ifelse(prob <= 0.05, "red", "black"),
     xlab = paste("Site (", seq.start, ":", seq.length + seq.start - 1, ")",
                  sep = ""),
     lwd  = 2)
lines(c(1, seq.length + 1), c(-log10(0.05), -log10(0.05)), lwd=1)
lines(c(1, seq.length + 1), c(-log10(0.01), -log10(0.01)), lwd=2)
lines(c(1, seq.length + 1), c(-log10(0.001), -log10(0.001)), lwd=3)
# Plotting of hystogramm of all registered P-values
#hist(-log10(prob), ylim=c(0, 20), breaks=200, freq=T)

# Cleaning of workaround
print(Sys.time()-start.time)

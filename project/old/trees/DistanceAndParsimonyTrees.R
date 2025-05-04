#following along 25Feb2025-PP563: distance and parsimony based

#installing the packages required for both methods
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
install.packages("seqinr", dep = TRUE)

#load the packages
library(ape)
library(adegenet)
library(seqinr)
library(phangorn)

#convert fasta file to alignment object
MSA_MAFFT <- read.alignment(file="mafft-aligned-seqs.fasta", format="fasta")
MSA_TCoffee <- read.alignment(file="allseqs.aln", format="clustal")

#DISTANCE BASED TREE
##compute the distances between sequences
D_MAFFT <- dist.alignment(MSA_MAFFT)
D_TCoffee <- dist.alignment(MSA_TCoffee)

##make the trees then reorganize with the ladder
tre_MAFFT <- nj(D_MAFFT)
Ltre_MAFFT <- ladderize(tre_MAFFT)
tre_TCoffee <- nj(D_TCoffee)
Ltre_TCoffee <- ladderize(tre_TCoffee)

##plot!
plot(Ltre_MAFFT, cex=0.6)
plot(Ltre_TCoffee, cex=0.6)
  
#PARSIMONY BASED TREE
##make the alignment objects into phangorn objects
pMSA_MAFFT <- as.phyDat(MSA_MAFFT)  
pMSA_TCoffee <- as.phyDat(MSA_TCoffee)

##make the guiding tree to start searching tree space from
tre.iniM <- nj(dist.alignment(MSA_MAFFT))
parsimony(tre.iniM, pMSA_MAFFT)
tre.iniT <- nj(dist.alignment(MSA_TCoffee))
parsimony(tre.iniT, pMSA_TCoffee)

##search the tree space!
tre.parsM <- optim.parsimony(tre.iniM, pMSA_MAFFT)
tre.parsT <- optim.parsimony(tre.iniT, pMSA_TCoffee)

#ladderize then plot the trees
PtreM <- ladderize(tre.parsM)
plot(PtreM, cex=0.6)
PtreT <- ladderize(tre.parsT)
plot(PtreT, cex=0.6)

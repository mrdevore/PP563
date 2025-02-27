#following along 25Feb2025-PP563: distance and parsimony based

#installing the packages required for both methods
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

#load the packages
library(ape)
library(adegenet)
library(phangorn)

#add the dataset, this one is a toy data set for practice (influenza data from 1993-2008)
dna <- fasta2DNAbin(file="http://adegenet.r-forge.r-project.org/files/usflu.fasta")

# >dna : outputs what is stored in "dna"
# >as.character(dna)[1:5,1:10] shows first five sequences and their first 10 bases as a table

#can download auxillary data that may help with making sense of the results (annotations)
##this is a data frame that has information on the samples (where they're from, the year they were sampled, etc)
## >table(annot$year) : summarizes data in a table based on the column selected (year in this example)

--------------------------
  
#DISTANCE BASED TREE
##compute the distances between sequences using the TN93 model (models change the tree bc they adjust the rate of mutations)
D <- dist.dna(dna, model="TN93")

##make the tree
tre <- nj(D)

##reorganize the tree by laddering
tre <- ladderize(tre)

##plot the tree
plot(tre, cex=.6)
title("A simple NJ tree")

#can make the tree prettier by following the in depth tutorial linked in CSL's tutorial
##rooting the tree:
tre2 <- root(tre,out=1)
tre2 <- ladderize(tre2)
plot(tre2,cex=0.6) #cex=text size

---------------------------------

#PARSIMONY BASED TREE
##convert the (dna) to a "phangorn object"
dna2 <- as.phyDat(dna)

##make a tree to start searching the tree space from and compute the parsimony for that tree
tre.ini <- nj(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)

##search for the maximum parsimony tree
tre.pars <- optim.parsimony(tre.ini, dna2)

##plot the tree :)
plot(tre.pars, cex=0.6)

##ew! ladderize? then plot again
Ptre <- ladderize(tre.pars)
plot(Ptre, cex=0.6)

##plot the tree unrooted
plot(tre.pars, type="unr", cex=0.6)

#yay :-)

#can analyze trees using things from the in depth tutorial (bootstrap values, how well it fits, etc)

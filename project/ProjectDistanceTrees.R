#load required packages for distance trees
library(ape)
library(adegenet)
library(phangorn)
library(seqinr)

#T COFFEE FIRST
#load alignment as object
MSA_TCoffee <- read.alignment(file="allseqs.aln", format="clustal")
class(MSA_TCoffee)

#compute the distances
D_TCoffee <- dist.alignment(MSA_TCoffee)

#make the tree and ladderize
tree_TCoffee <- nj(D_TCoffee)
Ltree_TCoffee <- ladderize(tree_TCoffee)

#plot the tree
plot(Ltree_TCoffee, cex=0.6)

#is it rooted?
is.rooted(Ltree_TCoffee)

#no, let's root it with the outgroup
RootLtree_TCoffee <- root(Ltree_TCoffee, outgroup= "Faustovirus")

#that didn't work, need to add another argument
RootLtree_TCoffee <- root(Ltree_TCoffee, outgroup = "Faustovirus", resolve.root = TRUE)
is.rooted(RootLtree_TCoffee)
plot(RootLtree_TCoffee, cex = 0.8)

#MAFFT NOW

MSA_MAFFT <- read.alignment(file="mafft-aligned-seqs.fasta", format="fasta")
D_MAFFT <- dist.alignment(MSA_MAFFT)
tree_MAFFT <- nj(D_MAFFT)
Ltree_MAFFT <- ladderize(tree_MAFFT)
plot(Ltree_MAFFT, cex=0.6)
RootLtree_MAFFT <- root(Ltree_MAFFT, outgroup = "Faustovirus", resolve.root = TRUE)
plot(RootLtree_MAFFT, cex=0.8)

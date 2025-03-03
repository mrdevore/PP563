#Expansin protein alignment project

Description: I am interested in feeding site formation in root-knot nematode infections of plant roots. Feeding site establishment requires manipulation 
of root cells by the nematode to induce giant cell (GC) formation. This manipulation includes cell wall changes. Expansins are proteins that loosen 
noncovalent connections in plant cell walls. My data set contains sequences from genes encoding characterized and putative expansins and expansin-like 
proteins from plant-parasitic nematodes, bacteria, fungi, and protists.

Source paper: Danchin, E. G. J. et al. Multiple lateral gene transfers and duplications have promoted plant parasitism  ability in nematodes. Proceedings 
    of the National Academy of Sciences 107, 17651–17656 (2010).
    Link: https://www.pnas.org/doi/10.1073/pnas.1008486107
    
Sequence source links:
    M. hapla sequences: https://brcwebportal.cos.ncsu.edu/medicago-hapla/mhapla.php / www.hapla.org
    M. incognita sequences: https://parasite.wormbase.org/index.html
    Other sequences: https://www.ncbi.nlm.nih.gov
    Reference accession numbers from Danchin (2010) Supplementary Info, pg. 100-101

06 February 2025
    - Updated notebook-log.md to have background information and relevant links
    - Created a folder for sequence files: Desktop/PP563/project/sequences
    
11 February 2025   
    - Downloaded sequences by copying peptide FASTAs into .txt files, files named according to identifier given by Danchin (2010) supplementary Dataset 
        S06, Accession column
    - Attempted to download all sequences in Danchin S06 table, samples lacking peptide sequences or from obsolete websites were omitted
    - The website listed in Danchin (2010) for M. hapla sequences no longer supports genome browsing and sequence downloads. Since M. hapla is my organism
        of interest in my work, I ran a BLASTP against all nematode sequences in WormBase ParaSite using M. javanica Expansin B as my input sequence
                Input peptide seq: MFAYLFIPIIIIISSSFSHGVPLNQQFTGSFTFYNDKGFGACGQQINAETEMLVAISHTQWIGGNPNNDPICRNICLKVDYKGK
                Output: 51 results, none from M. hapla
        - Since none of the sequences were from M. hapla, I did not add any to my sequences folder. I will continue with the sequences downloaded 
            according to Danchin (2010)

27 February 2025
    - After reading the MSA paper, it seems that none of the reported tools are a good fit for my data. I am aligning protein sequences, and the tools
    used in the paper are meant to align genomic sequences. So, I will look for protein multiple sequence alignment softwares.
        - I will base my MSA software decisions on this paper: https://doi.org/10.1007/978-1-59745-398-1_25
            - Based on the decision tree in this paper (Fig 25.3) I will use T-Coffee and MAFFT (L-ins-i). These tools are recommended for protein
            sequences with low global homology because they are better at finding local homology. My protein sequences come from a wide range of organisms
            from across the tree of life and are highly variable, so a program targeting local homology should better answer my question of whether plant-parasitic 
            nematodes have proteins that function like known expansins.

    - Software notes:
        - T-Coffee
            Description: from developers: "T-Coffee is a multiple sequence alignment package. You can use T-Coffee to align sequences or to combine the output of 
                your favorite alignment methods (Clustal, Mafft, Probcons, Muscle, etc.) into one unique alignment (M-coffee). T-Coffee can align Protein, DNA and 
                RNA sequences. It is also able to combine sequence information with protein structural information (Expresso), profile information (PSI-Coffee) or 
                RNA secondary structures (R-Coffee).
            Strengths: Accuracy, uses many tools at once
            Weaknesses: Slow
            Assumptions: pairwise alignments can lead to a better tree
            User Choices: 

        - MAFFT (L-ins-i)
            Description: From developers: "MAFFT is a multiple sequence alignment program for unix-like operating systems.  It offers a range of multiple 
                alignment methods, L-INS-i (accurate; recommended for <200 sequences), FFT-NS-2 (fast; recommended for >2,000 sequences), etc. "
            Strengths: reported to be highly accurate, especially for 10-100 protein sequences; able to align many and long sequences; allows for large gaps
                in the alignmenet (good for my data!)
            Weaknesses: does not work well for two very long and very unrelated sequences; assumes alignable domains are in a conserved order rather than rearranged
            Assumptions: ratio of transversions to transitions is 2; all sequences are homologous
            User Choices: https://mafft.cbrc.jp/alignment/software/manual/manual.html
    
    - Ran T-Coffee and MAFFT(L-ins-i) on sequences, saved to alignments folder
        T-coffee script: t_coffee allseqq.fasta
        MAFFT script: mafft-linsi allseqs.fasta > mafft-aligned-seqs.fasta

28 February 2025 - Assembling distance and parsimony based trees for my data
    
    I am running into an issues while trying to add my alignments to objects. In class we learned how to add DNA alignments, however I am working with protein
    sequences. I am working through the adegenet package documentation to fix this.

    Error when trying to add MAFFT alignment:
        >MSA_MAFFT <- alignment2genind("mafft-aligned-seqs.fasta")
        Error in alignment2genind("mafft-aligned-seqs.fasta") : 
            x is not a alignment object

    I installed and loaded the "seqinr" package from adegenet, because it seems I will need that to load my objects.

    Able to load MAFFT alignment into a genind object for seqinr by following this code:
        > MSA_MAFFT <- fasta2DNAbin(file="mafft-aligned-seqs.fasta")
        > ExpMSA_MAFFT <- DNAbin2genind(MSA_MAFFT)

    Attempted to compute distances:
        > D <- dist.alignment(ExpMSA_MAFFT)
        Error in dist.alignment(ExpMSA_MAFFT) : 
             Object of class 'alignment' expected
    
    Everything was broken so I started over! Converted my MAFFT alignment FASTA file into an "alignment" object:
        > MSA_MAFFT <- read.alignment(file="maffta-aligned-seqs.fasta", format="fasta")

    Checked the object class:
        > class(MSA_MAFFT)
        [1] "alignment"

    Moving on! Going to add the T-Coffee alignment next. After looking at the as.alignment help page in R, found that the .aln file from T-Coffee output is in ClustalW format. I added it as an alignment object then double checked the class.
        > MSA_TCoffee <- read.alignment(file="allseqs.aln", format="clustal")
        > class(MSA_TCoffee)
        [1] "alignment"

    Next I need to compute the distances. In lecture we used dist.dna, however I am aligning protein sequences so I will look through the documentation and the help pages in R to find the right functions for my alignments.
        - I will try to use dist.alignment
        Code: it worked!
        > D_MAFFT <- dist.alignment(MSA_MAFFT)
        > D_TCoffee <- dist.alignment(MSA_TCoffee)

    Continuing 03 March 2025 - I am making the trees and ladderizing them
        > tre_MAFFT <- nj(D_MAFFT)
        > Ltre_MAFFT <- ladderize(tre_MAFFT)
        > tre_TCoffee <- nj(D_TCoffee)
        > Ltre_TCoffee <- ladderize(tre_TCoffee)
        > plot(Ltre_MAFFT, cex=0.6)
        > plot(Ltre_TCoffee, cex=0.6)
        Everything worked! I saved the trees to my project folder as jpeg files.
            
    Now for the parsimony tree: following the class example with a single change --> instead of dist.dna, I am using dist.alignment to calculate the distances between my aligned amino acid sequences
        Transforming to phangorn objects:
        > pMSA_MAFFT <- as.phyDat(MSA_MAFFT)  
        > pMSA_TCoffee <- as.phyDat(MSA_TCoffee)

        Making the initial trees for searching tree space and computing their parsimony scores.
        >tre.iniM <- nj(dist.alignment(MSA_MAFFT))
        >parsimony(tre.iniM, pMSA_MAFFT)
        [1] 140
        
        >tre.iniT <- nj(dist.alignment(MSA_TCoffee))
        >parsimony(tre.iniT, pMSA_TCoffee)
        [1] 137

        Searching the tree space for the optimal tree based on parsimony scores.
        > tre.parsM <- optim.parsimony(tre.iniM, pMSA_MAFFT)
        Final p-score 129 after  4 nni operations 
        > tre.parsT <- optim.parsimony(tre.iniT, pMSA_TCoffee)
        Final p-score 122 after  2 nni operations 

        Ladderize and plot the trees. Everything worked, so I saved them as jpeg files in my project folder.
        > PtreM <- ladderize(tre.parsM)
        > plot(PtreM, cex=0.6)
        > PtreT <- ladderize(tre.parsT)
        > plot(PtreT, cex=0.6)
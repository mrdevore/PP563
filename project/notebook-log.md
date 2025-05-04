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
        T-coffee script: t_coffee allseqs.fasta
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
        > MSA_MAFFT <- read.alignment(file="mafft-aligned-seqs.fasta", format="fasta")

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

27 March 2025 - Maximum Likelihood Tree with RAxML
 - RAxML
        Description: from developers: "RAxML (Randomized Axelerated Maximum Likelihood) is a program for sequential and parallel Maximum Likelihood based inference of large phylogenetic trees. It can also be used for post- analyses of sets of phylogenetic trees, analyses of alignments and, evolutionary placement of short reads...It has originally been derived from fastDNAml which in turn was derived from Joe Felsentein’s dnaml which is part of the PHYLIP package."
        Strengths: Fast with minimal computational needs, accurate
        Weaknesses: Need to check MSA first
        Assumptions: 
        User Choices: 

    Checked MSA readability:
        > raxml-ng --check --msa mafft-aligned-seqs.fasta --model JTT+G --prefix pT1
    Output:
      RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 27-Mar-2025 11:50:23 as follows:

        raxml-ng --check --msa mafft-aligned-seqs.fasta --model JTT+G --prefix pT1

        Analysis options:
        run mode: Alignment validation
        start tree(s): 
        random seed: 1743094223
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: mafft-aligned-seqs.fasta
        [00:00:00] Loaded alignment with 47 taxa and 1116 sites

        NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) 
        NOTE: was saved to: /Users/melettedevore/Desktop/PP563/project/alignments/pT1.raxml.reduced.phy

        ERROR: Following taxon name contains invalid characters: Minc3s01215g21782:Minc3s01215g21782 peptide: Minc3s01215g21782 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: Minc3s00722g16527:Minc3s00722g16527 peptide: Minc3s00722g16527 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: Minc3s05299g37967 peptide: Minc3s05299g37967 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: Minc3s00430g12183:Minc3s00430g12183 peptide: Minc3s00430g12183 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: EEP81937.1 predicted protein [Uncinocarpus reesii 1704]
        ERROR: Following taxon name contains invalid characters: EDR15021.1 expansin family protein [Laccaria bicolor S238N-H82]
        ERROR: Following taxon name contains invalid characters: EDN21755.1 hypothetical protein BC1G_14954 [Botryotinia fuckeliana B05.10]
        ERROR: Following taxon name contains invalid characters: EDN92841.1 hypothetical protein SS1G_08706 [Sclerotinia sclerotiorum 1980 UF-70]
        ERROR: Following taxon name contains invalid characters: EDK04957.1 hypothetical protein MGG_11224 [Magnaporthe oryzae 70-15]
        ERROR: Following taxon name contains invalid characters: Minc3s00143g05927:Minc3s00143g05927 peptide: Minc3s00143g05927 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: EAW17746.1 riboflavin aldehyde-forming enzyme [Aspergillus fischeri NRRL 181]
        ERROR: Following taxon name contains invalid characters: EAU93049.1 hypothetical protein CC1G_06769 [Coprinopsis cinerea okayama7#130]
        ERROR: Following taxon name contains invalid characters: XP_001244656.1 uncharacterized protein CIMG_04097 [Coccidioides immitis RS]
        ERROR: Following taxon name contains invalid characters: EAU90886.1 hypothetical protein CC1G_02273 [Coprinopsis cinerea okayama7#130]
        ERROR: Following taxon name contains invalid characters: EAU88039.2 B2-aldehyde-forming enzyme [Coprinopsis cinerea okayama7#130]
        ERROR: Following taxon name contains invalid characters: CAC84564.1 EXPB2 protein [Globodera rostochiensis]
        ERROR: Following taxon name contains invalid characters: BAG16533.1 expansin-like protein [Bursaphelenchus xylophilus]
        ERROR: Following taxon name contains invalid characters: EAU88038.1 hypothetical protein CC1G_10811 [Coprinopsis cinerea okayama7#130]
        ERROR: Following taxon name contains invalid characters: Minc3s00722g16529:Minc3s00722g16529 peptide: Minc3s00722g16529 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: Minc3s00264g08913:Minc3s00264g08913 peptide: Minc3s00264g08913 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: EAL65166.1 hypothetical protein DDB_G0284479 [Dictyostelium discoideum AX4]
        ERROR: Following taxon name contains invalid characters: Minc3s03843g34959:Minc3s03843g34959 peptide: Minc3s03843g34959 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: CAG79677.1 YALI0E17941p [Yarrowia lipolytica CLIB122]
        ERROR: Following taxon name contains invalid characters: CAP59537.1 putative avirulence protein precursor, partial [Meloidogyne javanica]
        ERROR: Following taxon name contains invalid characters: EAL91558.1 riboflavin aldehyde-forming enzyme [Aspergillus fumigatus Af293]
        ERROR: Following taxon name contains invalid characters: Minc3s01732g25926:Minc3s01732g25926 peptide: Minc3s01732g25926 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: CAP59535.1 putative pathogenicity factor precursor, partial [Meloidogyne javanica]
        ERROR: Following taxon name contains invalid characters: Minc3s01671g25495:Minc3s01671g25495 peptide: Minc3s01671g25495 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: EAL65201.2 speract/scavenger receptor domain-containing protein [Dictyostelium discoideum AX4]
        ERROR: Following taxon name contains invalid characters: CAC83611.1 EXPB1 protein [Globodera rostochiensis]
        ERROR: Following taxon name contains invalid characters: EAK86428.1 hypothetical protein UM05495.1 [Ustilago maydis 521]
        ERROR: Following taxon name contains invalid characters: XP_504084.2 uncharacterized protein YALI1_E21282g [Yarrowia lipolytica]
        ERROR: Following taxon name contains invalid characters: CAC42207.1 hypothetical protein [Amycolatopsis mediterranei U32]
        ERROR: Following taxon name contains invalid characters: EAL04023.1 hypothetical protein CaO19.4688 [Candida albicans SC5314]
        ERROR: Following taxon name contains invalid characters: BAH80450.1 putative riboflavin aldehyde-forming enzyme [Lentinula edodes]
        ERROR: Following taxon name contains invalid characters: Minc3s00489g13197 peptide: Minc3s00489g13197 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: BAE66341.1 unnamed protein product [Aspergillus oryzae RIB40]
        ERROR: Following taxon name contains invalid characters: ABV60414.1 expansin [Bursaphelenchus xylophilus]
        ERROR: Following taxon name contains invalid characters: BAG16534.1 expansin-like protein [Bursaphelenchus mucronatus]
        ERROR: Following taxon name contains invalid characters: CAP59538.1 putative avirulence protein precursor, partial [Meloidogyne javanica]
        ERROR: Following taxon name contains invalid characters: Minc3s02572g30704:Minc3s02572g30704 peptide: Minc3s02572g30704 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: ACU36800.1 Rare lipoprotein A [Actinosynnema mirum DSM 43827]
        ERROR: Following taxon name contains invalid characters: CAP59536.1 putative pathogenicity factor precursor, partial [Meloidogyne javanica]
        ERROR: Following taxon name contains invalid characters: ACN58322.1 expansin B2, partial [Globodera rostochiensis]
        ERROR: Following taxon name contains invalid characters: Minc3s02154g28582:Minc3s02154g28582 peptide: Minc3s02154g28582 pep:protein_coding
        ERROR: Following taxon name contains invalid characters: CAX41987.1 unnamed protein product [Candida dubliniensis CD36]
        ERROR: Following taxon name contains invalid characters: AAD32751.1 unknown [Streptomyces lavendulae]

        NOTE: Following symbols are not allowed in taxa names to ensure Newick compatibility:
        NOTE: " " (space), ";" (semicolon), ":" (colon), "," (comma), "()" (parentheses), "'" (quote). 
        NOTE: Please either correct the names manually, or use the reduced alignment file
        NOTE: generated by RAxML-NG (see above).

        ERROR: Alignment check failed (see details above)!

    Ran check again with reduced file made by RAxML:
        > raxml-ng --check --msa pT1.raxml.reduced.phy --model JTT+G --prefix pT1
    Output:
        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 27-Mar-2025 11:51:00 as follows:

        raxml-ng --check --msa pT1.raxml.reduced.phy --model JTT+G --prefix pT1

        Analysis options:
        run mode: Alignment validation
        start tree(s): 
        random seed: 1743094260
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: pT1.raxml.reduced.phy
        [00:00:00] Loaded alignment with 47 taxa and 1116 sites

        Alignment comprises 1 partitions and 1116 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 1116 / 1116
        Gaps: 77.77 %
        Invariant sites: 32.17 %


        Alignment can be successfully read by RAxML-NG.


        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/alignments/pT1.raxml.log

        Analysis started: 27-Mar-2025 11:51:00 / finished: 27-Mar-2025 11:51:00

        Elapsed time: 0.007 seconds
    
    Ran RAxML with parse as recommended for large datasets:
        > raxml-ng --parse --msa pT1.raxml.reduced.phy --model JTT+G --prefix pT2
    Output:       
        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 27-Mar-2025 11:56:17 as follows:

        raxml-ng --parse --msa pT1.raxml.reduced.phy --model JTT+G --prefix pT2

        Analysis options:
        run mode: Alignment parsing and compression
        start tree(s): 
        random seed: 1743094577
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: pT1.raxml.reduced.phy
        [00:00:00] Loaded alignment with 47 taxa and 1116 sites

        Alignment comprises 1 partitions and 892 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 1116 / 892
        Gaps: 77.77 %
        Invariant sites: 32.17 %


        NOTE: Binary MSA file created: pT2.raxml.rba

        * Estimated memory requirements                : 51 MB

        * Recommended number of threads / MPI processes: 9

        Please note that numbers given above are rough estimates only. 
        Actual memory consumption and parallel performance on your system may differ!

        Alignment can be successfully read by RAxML-NG.


        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/alignments/pT2.raxml.log

        Analysis started: 27-Mar-2025 11:56:17 / finished: 27-Mar-2025 11:56:17

        Elapsed time: 0.010 seconds 

    Ran the tree search: 
        > raxml-ng --msa pT1.raxml.reduced.phy --model JTT+G --prefix T3 --threads 2 --seed 2
    Output:
        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 27-Mar-2025 11:58:44 as follows:

        raxml-ng --msa pT1.raxml.reduced.phy --model JTT+G --prefix T3 --threads 2 --seed 2

        Analysis options:
        run mode: ML tree search
        start tree(s): random (10) + parsimony (10)
        random seed: 2
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        fast spr radius: AUTO
        spr subtree cutoff: 1.000000
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: pT1.raxml.reduced.phy
        [00:00:00] Loaded alignment with 47 taxa and 1116 sites

        Alignment comprises 1 partitions and 892 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 1116 / 892
        Gaps: 77.77 %
        Invariant sites: 32.17 %


        NOTE: Binary MSA file created: T3.raxml.rba

        [00:00:00] Generating 10 random starting tree(s) with 47 taxa
        [00:00:00] Generating 10 parsimony starting tree(s) with 47 taxa
        [00:00:00] Data distribution: max. partitions/sites/weight per thread: 1 / 446 / 35680

        Starting ML tree search with 20 distinct starting trees

        [00:00:00 -28495.709683] Initial branch length optimization
        [00:00:00 -24952.317804] Model parameter optimization (eps = 10.000000)
        [00:00:00 -24945.220271] AUTODETECT spr round 1 (radius: 5)
        [00:00:02 -20111.585196] AUTODETECT spr round 2 (radius: 10)
        [00:00:06 -19779.017101] AUTODETECT spr round 3 (radius: 15)
        [00:00:08 -19779.009917] SPR radius for FAST iterations: 10 (autodetect)
        [00:00:08 -19779.009917] Model parameter optimization (eps = 3.000000)
        [00:00:08 -19771.801338] FAST spr round 1 (radius: 10)
        [00:00:18 -19549.875633] FAST spr round 2 (radius: 10)
        [00:00:22 -19544.825062] FAST spr round 3 (radius: 10)
        [00:00:27 -19544.824910] Model parameter optimization (eps = 1.000000)
        [00:00:27 -19542.136733] SLOW spr round 1 (radius: 5)
        [00:00:40 -19540.439728] SLOW spr round 2 (radius: 5)
        [00:00:54 -19540.418599] SLOW spr round 3 (radius: 10)
        [00:01:04 -19540.417555] SLOW spr round 4 (radius: 15)
        [00:01:11 -19540.417387] SLOW spr round 5 (radius: 20)
        [00:01:12 -19540.417357] SLOW spr round 6 (radius: 25)
        [00:01:13 -19540.417351] Model parameter optimization (eps = 0.100000)

        [00:01:13] ML tree search #1, logLikelihood: -19540.356262

        [00:01:13 -29196.091804] Initial branch length optimization
        [00:01:13 -25223.210255] Model parameter optimization (eps = 10.000000)
        [00:01:14 -25219.035559] AUTODETECT spr round 1 (radius: 5)
        [00:01:19 -20120.869033] AUTODETECT spr round 2 (radius: 10)
        [00:01:23 -19785.002276] AUTODETECT spr round 3 (radius: 15)
        [00:01:25 -19784.987450] SPR radius for FAST iterations: 10 (autodetect)
        [00:01:25 -19784.987450] Model parameter optimization (eps = 3.000000)
        [00:01:26 -19748.236919] FAST spr round 1 (radius: 10)
        [00:01:33 -19552.823753] FAST spr round 2 (radius: 10)
        [00:01:37 -19544.642221] FAST spr round 3 (radius: 10)
        [00:01:42 -19543.386434] FAST spr round 4 (radius: 10)
        [00:01:47 -19543.050431] FAST spr round 5 (radius: 10)
        [00:01:52 -19542.687043] FAST spr round 6 (radius: 10)
        [00:01:58 -19542.679205] Model parameter optimization (eps = 1.000000)
        [00:01:58 -19541.553403] SLOW spr round 1 (radius: 5)
        [00:02:11 -19540.361437] SLOW spr round 2 (radius: 5)
        [00:02:23 -19540.360179] SLOW spr round 3 (radius: 10)
        [00:02:34 -19540.359964] SLOW spr round 4 (radius: 15)
        [00:02:40 -19540.359923] SLOW spr round 5 (radius: 20)
        [00:02:40 -19540.359915] SLOW spr round 6 (radius: 25)
        [00:02:41 -19540.359913] Model parameter optimization (eps = 0.100000)

        [00:02:41] ML tree search #2, logLikelihood: -19540.349382

        [00:02:41 -29564.149134] Initial branch length optimization
        [00:02:41 -24988.600702] Model parameter optimization (eps = 10.000000)
        [00:02:41 -24984.376177] AUTODETECT spr round 1 (radius: 5)
        [00:02:44 -20709.430421] AUTODETECT spr round 2 (radius: 10)
        [00:02:48 -20274.227283] AUTODETECT spr round 3 (radius: 15)
        [00:02:50 -20274.221146] SPR radius for FAST iterations: 10 (autodetect)
        [00:02:50 -20274.221146] Model parameter optimization (eps = 3.000000)
        [00:02:51 -20243.949068] FAST spr round 1 (radius: 10)
        [00:02:57 -19550.826452] FAST spr round 2 (radius: 10)
        [00:03:04 -19543.682013] FAST spr round 3 (radius: 10)
        [00:03:09 -19543.681842] Model parameter optimization (eps = 1.000000)
        [00:03:09 -19539.382300] SLOW spr round 1 (radius: 5)
        [00:03:22 -19539.321421] SLOW spr round 2 (radius: 10)
        [00:03:31 -19539.320238] SLOW spr round 3 (radius: 15)
        [00:03:38 -19539.320167] SLOW spr round 4 (radius: 20)
        [00:03:40 -19539.320156] SLOW spr round 5 (radius: 25)
        [00:03:41 -19539.320154] Model parameter optimization (eps = 0.100000)

        [00:03:41] ML tree search #3, logLikelihood: -19539.233355

        [00:03:41 -29765.514038] Initial branch length optimization
        [00:03:41 -25203.614948] Model parameter optimization (eps = 10.000000)
        [00:03:43 -24988.381774] AUTODETECT spr round 1 (radius: 5)
        [00:03:45 -20195.080538] AUTODETECT spr round 2 (radius: 10)
        [00:03:49 -19794.506474] AUTODETECT spr round 3 (radius: 15)
        [00:03:52 -19794.459146] SPR radius for FAST iterations: 10 (autodetect)
        [00:03:52 -19794.459146] Model parameter optimization (eps = 3.000000)
        [00:03:52 -19755.445573] FAST spr round 1 (radius: 10)
        [00:03:58 -19577.126986] FAST spr round 2 (radius: 10)
        [00:04:04 -19540.512408] FAST spr round 3 (radius: 10)
        [00:04:08 -19539.336165] FAST spr round 4 (radius: 10)
        [00:04:12 -19539.335775] Model parameter optimization (eps = 1.000000)
        [00:04:12 -19539.235164] SLOW spr round 1 (radius: 5)
        [00:04:23 -19539.231447] SLOW spr round 2 (radius: 10)
        [00:04:33 -19539.231372] SLOW spr round 3 (radius: 15)
        [00:04:40 -19539.231365] SLOW spr round 4 (radius: 20)
        [00:04:43 -19539.231364] SLOW spr round 5 (radius: 25)
        [00:04:43 -19539.231363] Model parameter optimization (eps = 0.100000)

        [00:04:43] ML tree search #4, logLikelihood: -19539.222148

        [00:04:43 -29227.095222] Initial branch length optimization
        [00:04:44 -25128.598257] Model parameter optimization (eps = 10.000000)
        [00:04:44 -25120.684033] AUTODETECT spr round 1 (radius: 5)
        [00:04:46 -20010.342154] AUTODETECT spr round 2 (radius: 10)
        [00:04:50 -19802.946972] AUTODETECT spr round 3 (radius: 15)
        [00:04:51 -19802.934771] SPR radius for FAST iterations: 10 (autodetect)
        [00:04:51 -19802.934771] Model parameter optimization (eps = 3.000000)
        [00:04:52 -19755.941140] FAST spr round 1 (radius: 10)
        [00:04:58 -19545.664672] FAST spr round 2 (radius: 10)
        [00:05:02 -19543.943839] FAST spr round 3 (radius: 10)
        [00:05:06 -19542.022280] FAST spr round 4 (radius: 10)
        [00:05:10 -19542.018412] Model parameter optimization (eps = 1.000000)
        [00:05:10 -19541.581543] SLOW spr round 1 (radius: 5)
        [00:05:24 -19540.381211] SLOW spr round 2 (radius: 5)
        [00:05:34 -19540.380651] SLOW spr round 3 (radius: 10)
        [00:05:44 -19540.380572] SLOW spr round 4 (radius: 15)
        [00:05:51 -19540.380556] SLOW spr round 5 (radius: 20)
        [00:05:51 -19540.380553] SLOW spr round 6 (radius: 25)
        [00:05:52 -19540.380553] Model parameter optimization (eps = 0.100000)

        [00:05:52] ML tree search #5, logLikelihood: -19540.352091

        [00:05:52 -29950.652560] Initial branch length optimization
        [00:05:52 -25112.695059] Model parameter optimization (eps = 10.000000)
        [00:05:53 -25106.974358] AUTODETECT spr round 1 (radius: 5)
        [00:05:56 -21180.342170] AUTODETECT spr round 2 (radius: 10)
        [00:06:00 -19932.899193] AUTODETECT spr round 3 (radius: 15)
        [00:06:02 -19932.888256] SPR radius for FAST iterations: 10 (autodetect)
        [00:06:02 -19932.888256] Model parameter optimization (eps = 3.000000)
        [00:06:03 -19899.009707] FAST spr round 1 (radius: 10)
        [00:06:09 -19542.516937] FAST spr round 2 (radius: 10)
        [00:06:14 -19541.887559] FAST spr round 3 (radius: 10)
        [00:06:18 -19541.285814] FAST spr round 4 (radius: 10)
        [00:06:23 -19541.194830] Model parameter optimization (eps = 1.000000)
        [00:06:24 -19539.286449] SLOW spr round 1 (radius: 5)
        [00:06:37 -19539.264477] SLOW spr round 2 (radius: 10)
        [00:06:46 -19539.263392] SLOW spr round 3 (radius: 15)
        [00:06:54 -19539.263231] SLOW spr round 4 (radius: 20)
        [00:06:56 -19539.263199] SLOW spr round 5 (radius: 25)
        [00:06:57 -19539.263193] Model parameter optimization (eps = 0.100000)

        [00:06:57] ML tree search #6, logLikelihood: -19539.226147

        [00:06:57 -30183.008073] Initial branch length optimization
        [00:06:57 -25646.736865] Model parameter optimization (eps = 10.000000)
        [00:06:57 -25637.891956] AUTODETECT spr round 1 (radius: 5)
        [00:07:00 -20734.743645] AUTODETECT spr round 2 (radius: 10)
        [00:07:04 -19964.095326] AUTODETECT spr round 3 (radius: 15)
        [00:07:05 -19964.068545] SPR radius for FAST iterations: 10 (autodetect)
        [00:07:05 -19964.068545] Model parameter optimization (eps = 3.000000)
        [00:07:06 -19930.223596] FAST spr round 1 (radius: 10)
        [00:07:13 -19548.433642] FAST spr round 2 (radius: 10)
        [00:07:17 -19545.086392] FAST spr round 3 (radius: 10)
        [00:07:21 -19544.074805] FAST spr round 4 (radius: 10)
        [00:07:26 -19543.477529] FAST spr round 5 (radius: 10)
        [00:07:30 -19543.406830] Model parameter optimization (eps = 1.000000)
        [00:07:30 -19539.365754] SLOW spr round 1 (radius: 5)
        [00:07:41 -19539.317861] SLOW spr round 2 (radius: 10)
        [00:07:48 -19539.315742] SLOW spr round 3 (radius: 15)
        [00:07:55 -19539.315396] SLOW spr round 4 (radius: 20)
        [00:07:57 -19539.315326] SLOW spr round 5 (radius: 25)
        [00:07:57 -19539.315312] Model parameter optimization (eps = 0.100000)

        [00:07:57] ML tree search #7, logLikelihood: -19539.232652

        [00:07:57 -29703.255571] Initial branch length optimization
        [00:07:58 -25649.419412] Model parameter optimization (eps = 10.000000)
        [00:07:58 -25640.818393] AUTODETECT spr round 1 (radius: 5)
        [00:08:00 -20048.772036] AUTODETECT spr round 2 (radius: 10)
        [00:08:04 -19870.497140] AUTODETECT spr round 3 (radius: 15)
        [00:08:06 -19870.462594] SPR radius for FAST iterations: 10 (autodetect)
        [00:08:06 -19870.462594] Model parameter optimization (eps = 3.000000)
        [00:08:07 -19812.885738] FAST spr round 1 (radius: 10)
        [00:08:12 -19545.292029] FAST spr round 2 (radius: 10)
        [00:08:17 -19543.461823] FAST spr round 3 (radius: 10)
        [00:08:21 -19541.699830] FAST spr round 4 (radius: 10)
        [00:08:25 -19541.125571] FAST spr round 5 (radius: 10)
        [00:08:30 -19540.643057] FAST spr round 6 (radius: 10)
        [00:08:34 -19540.340066] FAST spr round 7 (radius: 10)
        [00:08:39 -19540.321169] Model parameter optimization (eps = 1.000000)
        [00:08:39 -19539.365249] SLOW spr round 1 (radius: 5)
        [00:08:50 -19539.324924] SLOW spr round 2 (radius: 10)
        [00:08:59 -19539.323818] SLOW spr round 3 (radius: 15)
        [00:09:05 -19539.323790] SLOW spr round 4 (radius: 20)
        [00:09:07 -19539.323789] SLOW spr round 5 (radius: 25)
        [00:09:08 -19539.323789] Model parameter optimization (eps = 0.100000)

        [00:09:08] ML tree search #8, logLikelihood: -19539.233802

        [00:09:08 -29046.872574] Initial branch length optimization
        [00:09:08 -24410.901642] Model parameter optimization (eps = 10.000000)
        [00:09:08 -24406.522751] AUTODETECT spr round 1 (radius: 5)
        [00:09:11 -20502.261363] AUTODETECT spr round 2 (radius: 10)
        [00:09:15 -19675.340589] AUTODETECT spr round 3 (radius: 15)
        [00:09:17 -19675.330856] SPR radius for FAST iterations: 10 (autodetect)
        [00:09:17 -19675.330856] Model parameter optimization (eps = 3.000000)
        [00:09:18 -19653.804771] FAST spr round 1 (radius: 10)
        [00:09:24 -19542.719241] FAST spr round 2 (radius: 10)
        [00:09:28 -19541.219705] FAST spr round 3 (radius: 10)
        [00:09:34 -19540.532919] FAST spr round 4 (radius: 10)
        [00:09:42 -19540.224758] FAST spr round 5 (radius: 10)
        [00:09:49 -19540.206243] Model parameter optimization (eps = 1.000000)
        [00:09:50 -19539.347366] SLOW spr round 1 (radius: 5)
        [00:10:06 -19539.313222] SLOW spr round 2 (radius: 10)
        [00:10:22 -19539.312501] SLOW spr round 3 (radius: 15)
        [00:10:33 -19539.312489] SLOW spr round 4 (radius: 20)
        [00:10:37 -19539.312488] SLOW spr round 5 (radius: 25)
        [00:10:37 -19539.312488] Model parameter optimization (eps = 0.100000)

        [00:10:38] ML tree search #9, logLikelihood: -19539.232642

        [00:10:38 -29395.935314] Initial branch length optimization
        [00:10:38 -25121.346289] Model parameter optimization (eps = 10.000000)
        [00:10:38 -25115.539572] AUTODETECT spr round 1 (radius: 5)
        [00:10:43 -20573.069172] AUTODETECT spr round 2 (radius: 10)
        [00:10:49 -20035.822494] AUTODETECT spr round 3 (radius: 15)
        [00:10:52 -20035.789228] SPR radius for FAST iterations: 10 (autodetect)
        [00:10:52 -20035.789228] Model parameter optimization (eps = 3.000000)
        [00:10:53 -20002.717878] FAST spr round 1 (radius: 10)
        [00:11:00 -19548.269475] FAST spr round 2 (radius: 10)
        [00:11:08 -19542.398798] FAST spr round 3 (radius: 10)
        [00:11:15 -19541.928682] FAST spr round 4 (radius: 10)
        [00:11:20 -19541.927945] Model parameter optimization (eps = 1.000000)
        [00:11:20 -19539.311606] SLOW spr round 1 (radius: 5)
        [00:11:35 -19539.279905] SLOW spr round 2 (radius: 10)
        [00:11:45 -19539.278212] SLOW spr round 3 (radius: 15)
        [00:11:54 -19539.277940] SLOW spr round 4 (radius: 20)
        [00:11:57 -19539.277886] SLOW spr round 5 (radius: 25)
        [00:11:58 -19539.277874] Model parameter optimization (eps = 0.100000)

        [00:11:58] ML tree search #10, logLikelihood: -19539.228458

        [00:11:58 -21481.489844] Initial branch length optimization
        [00:11:59 -19638.439359] Model parameter optimization (eps = 10.000000)
        [00:11:59 -19608.088194] AUTODETECT spr round 1 (radius: 5)
        [00:12:01 -19541.190632] AUTODETECT spr round 2 (radius: 10)
        [00:12:06 -19541.183460] SPR radius for FAST iterations: 5 (autodetect)
        [00:12:06 -19541.183460] Model parameter optimization (eps = 3.000000)
        [00:12:06 -19540.224364] FAST spr round 1 (radius: 5)
        [00:12:10 -19539.433762] FAST spr round 2 (radius: 5)
        [00:12:14 -19539.305136] FAST spr round 3 (radius: 5)
        [00:12:20 -19539.304537] Model parameter optimization (eps = 1.000000)
        [00:12:20 -19539.231502] SLOW spr round 1 (radius: 5)
        [00:12:38 -19539.228570] SLOW spr round 2 (radius: 10)
        [00:12:50 -19539.228487] SLOW spr round 3 (radius: 15)
        [00:12:57 -19539.228480] SLOW spr round 4 (radius: 20)
        [00:12:59 -19539.228479] SLOW spr round 5 (radius: 25)
        [00:13:00 -19539.228479] Model parameter optimization (eps = 0.100000)

        [00:13:01] ML tree search #11, logLikelihood: -19539.221814

        [00:13:01 -21521.638540] Initial branch length optimization
        [00:13:01 -19653.717605] Model parameter optimization (eps = 10.000000)
        [00:13:02 -19625.397936] AUTODETECT spr round 1 (radius: 5)
        [00:13:06 -19549.241678] AUTODETECT spr round 2 (radius: 10)
        [00:13:10 -19549.238048] SPR radius for FAST iterations: 5 (autodetect)
        [00:13:10 -19549.238048] Model parameter optimization (eps = 3.000000)
        [00:13:10 -19548.077763] FAST spr round 1 (radius: 5)
        [00:13:15 -19539.315885] FAST spr round 2 (radius: 5)
        [00:13:20 -19539.288816] Model parameter optimization (eps = 1.000000)
        [00:13:20 -19539.229672] SLOW spr round 1 (radius: 5)
        [00:13:33 -19539.226598] SLOW spr round 2 (radius: 10)
        [00:13:42 -19539.226397] SLOW spr round 3 (radius: 15)
        [00:13:49 -19539.226360] SLOW spr round 4 (radius: 20)
        [00:13:52 -19539.226353] SLOW spr round 5 (radius: 25)
        [00:13:53 -19539.226351] Model parameter optimization (eps = 0.100000)

        [00:13:53] ML tree search #12, logLikelihood: -19539.221512

        [00:13:53 -21433.465982] Initial branch length optimization
        [00:13:53 -19618.364379] Model parameter optimization (eps = 10.000000)
        [00:13:53 -19589.868483] AUTODETECT spr round 1 (radius: 5)
        [00:13:56 -19545.158718] AUTODETECT spr round 2 (radius: 10)
        [00:14:02 -19545.149567] SPR radius for FAST iterations: 5 (autodetect)
        [00:14:02 -19545.149567] Model parameter optimization (eps = 3.000000)
        [00:14:02 -19543.781062] FAST spr round 1 (radius: 5)
        [00:14:07 -19541.962889] FAST spr round 2 (radius: 5)
        [00:14:11 -19540.620404] FAST spr round 3 (radius: 5)
        [00:14:16 -19539.432834] FAST spr round 4 (radius: 5)
        [00:14:21 -19539.432524] Model parameter optimization (eps = 1.000000)
        [00:14:21 -19539.246543] SLOW spr round 1 (radius: 5)
        [00:14:33 -19539.240392] SLOW spr round 2 (radius: 10)
        [00:14:40 -19539.240218] SLOW spr round 3 (radius: 15)
        [00:14:49 -19539.240197] SLOW spr round 4 (radius: 20)
        [00:14:52 -19539.240193] SLOW spr round 5 (radius: 25)
        [00:14:53 -19539.240192] Model parameter optimization (eps = 0.100000)

        [00:14:53] ML tree search #13, logLikelihood: -19539.223377

        [00:14:53 -21432.123363] Initial branch length optimization
        [00:14:53 -19625.997163] Model parameter optimization (eps = 10.000000)
        [00:14:53 -19600.548678] AUTODETECT spr round 1 (radius: 5)
        [00:14:56 -19548.392224] AUTODETECT spr round 2 (radius: 10)
        [00:15:01 -19548.388032] SPR radius for FAST iterations: 5 (autodetect)
        [00:15:01 -19548.388032] Model parameter optimization (eps = 3.000000)
        [00:15:01 -19546.135980] FAST spr round 1 (radius: 5)
        [00:15:07 -19542.521729] FAST spr round 2 (radius: 5)
        [00:15:14 -19541.737953] FAST spr round 3 (radius: 5)
        [00:15:19 -19541.737327] Model parameter optimization (eps = 1.000000)
        [00:15:19 -19541.544011] SLOW spr round 1 (radius: 5)
        [00:15:31 -19540.376530] SLOW spr round 2 (radius: 5)
        [00:15:41 -19540.359036] SLOW spr round 3 (radius: 10)
        [00:15:52 -19540.358960] SLOW spr round 4 (radius: 15)
        [00:15:59 -19540.358957] SLOW spr round 5 (radius: 20)
        [00:16:00 -19540.358956] SLOW spr round 6 (radius: 25)
        [00:16:01 -19540.358956] Model parameter optimization (eps = 0.100000)

        [00:16:01] ML tree search #14, logLikelihood: -19540.349365

        [00:16:01 -21488.038475] Initial branch length optimization
        [00:16:01 -19657.823140] Model parameter optimization (eps = 10.000000)
        [00:16:02 -19625.699762] AUTODETECT spr round 1 (radius: 5)
        [00:16:05 -19540.841230] AUTODETECT spr round 2 (radius: 10)
        [00:16:11 -19540.833290] SPR radius for FAST iterations: 5 (autodetect)
        [00:16:11 -19540.833290] Model parameter optimization (eps = 3.000000)
        [00:16:11 -19539.986135] FAST spr round 1 (radius: 5)
        [00:16:15 -19539.254986] FAST spr round 2 (radius: 5)
        [00:16:19 -19539.254663] Model parameter optimization (eps = 1.000000)
        [00:16:20 -19539.225116] SLOW spr round 1 (radius: 5)
        [00:16:33 -19539.224101] SLOW spr round 2 (radius: 10)
        [00:16:44 -19539.224052] SLOW spr round 3 (radius: 15)
        [00:16:50 -19539.224045] SLOW spr round 4 (radius: 20)
        [00:16:52 -19539.224044] SLOW spr round 5 (radius: 25)
        [00:16:53 -19539.224044] Model parameter optimization (eps = 0.100000)

        [00:16:53] ML tree search #15, logLikelihood: -19539.221227

        [00:16:53 -21601.009683] Initial branch length optimization
        [00:16:53 -19665.014508] Model parameter optimization (eps = 10.000000)
        [00:16:53 -19636.677971] AUTODETECT spr round 1 (radius: 5)
        [00:16:55 -19547.371100] AUTODETECT spr round 2 (radius: 10)
        [00:16:59 -19547.359066] SPR radius for FAST iterations: 5 (autodetect)
        [00:16:59 -19547.359066] Model parameter optimization (eps = 3.000000)
        [00:16:59 -19545.901369] FAST spr round 1 (radius: 5)
        [00:17:02 -19544.812265] FAST spr round 2 (radius: 5)
        [00:17:06 -19544.193581] FAST spr round 3 (radius: 5)
        [00:17:10 -19543.583625] FAST spr round 4 (radius: 5)
        [00:17:14 -19541.675261] FAST spr round 5 (radius: 5)
        [00:17:18 -19541.669878] Model parameter optimization (eps = 1.000000)
        [00:17:18 -19541.535150] SLOW spr round 1 (radius: 5)
        [00:17:28 -19540.372789] SLOW spr round 2 (radius: 5)
        [00:17:37 -19540.354521] SLOW spr round 3 (radius: 10)
        [00:17:45 -19540.354411] SLOW spr round 4 (radius: 15)
        [00:17:50 -19540.354397] SLOW spr round 5 (radius: 20)
        [00:17:51 -19540.354394] SLOW spr round 6 (radius: 25)
        [00:17:51 -19540.354394] Model parameter optimization (eps = 0.100000)

        [00:17:52] ML tree search #16, logLikelihood: -19540.348733

        [00:17:52 -21438.711566] Initial branch length optimization
        [00:17:52 -19647.063664] Model parameter optimization (eps = 10.000000)
        [00:17:52 -19620.622043] AUTODETECT spr round 1 (radius: 5)
        [00:17:54 -19542.981228] AUTODETECT spr round 2 (radius: 10)
        [00:17:58 -19542.974474] SPR radius for FAST iterations: 5 (autodetect)
        [00:17:58 -19542.974474] Model parameter optimization (eps = 3.000000)
        [00:17:58 -19541.449893] FAST spr round 1 (radius: 5)
        [00:18:02 -19539.379422] FAST spr round 2 (radius: 5)
        [00:18:06 -19539.349279] Model parameter optimization (eps = 1.000000)
        [00:18:06 -19539.237701] SLOW spr round 1 (radius: 5)
        [00:18:16 -19539.232813] SLOW spr round 2 (radius: 10)
        [00:18:23 -19539.232617] SLOW spr round 3 (radius: 15)
        [00:18:29 -19539.232589] SLOW spr round 4 (radius: 20)
        [00:18:31 -19539.232584] SLOW spr round 5 (radius: 25)
        [00:18:32 -19539.232583] Model parameter optimization (eps = 0.100000)

        [00:18:32] ML tree search #17, logLikelihood: -19539.222339

        [00:18:32 -21497.982372] Initial branch length optimization
        [00:18:32 -19636.702505] Model parameter optimization (eps = 10.000000)
        [00:18:32 -19602.484954] AUTODETECT spr round 1 (radius: 5)
        [00:18:35 -19549.127779] AUTODETECT spr round 2 (radius: 10)
        [00:18:38 -19549.119717] SPR radius for FAST iterations: 5 (autodetect)
        [00:18:38 -19549.119717] Model parameter optimization (eps = 3.000000)
        [00:18:38 -19548.022746] FAST spr round 1 (radius: 5)
        [00:18:42 -19541.912993] FAST spr round 2 (radius: 5)
        [00:18:45 -19540.663419] FAST spr round 3 (radius: 5)
        [00:18:49 -19540.410652] FAST spr round 4 (radius: 5)
        [00:18:53 -19540.404471] Model parameter optimization (eps = 1.000000)
        [00:18:53 -19540.355517] SLOW spr round 1 (radius: 5)
        [00:19:03 -19540.353409] SLOW spr round 2 (radius: 10)
        [00:19:12 -19540.353335] SLOW spr round 3 (radius: 15)
        [00:19:17 -19540.353327] SLOW spr round 4 (radius: 20)
        [00:19:17 -19540.353325] SLOW spr round 5 (radius: 25)
        [00:19:18 -19540.353325] Model parameter optimization (eps = 0.100000)

        [00:19:18] ML tree search #18, logLikelihood: -19540.348605

        [00:19:18 -21398.689759] Initial branch length optimization
        [00:19:18 -19618.183347] Model parameter optimization (eps = 10.000000)
        [00:19:18 -19588.544304] AUTODETECT spr round 1 (radius: 5)
        [00:19:21 -19546.456642] AUTODETECT spr round 2 (radius: 10)
        [00:19:24 -19546.448876] SPR radius for FAST iterations: 5 (autodetect)
        [00:19:24 -19546.448876] Model parameter optimization (eps = 3.000000)
        [00:19:25 -19544.890587] FAST spr round 1 (radius: 5)
        [00:19:28 -19541.822817] FAST spr round 2 (radius: 5)
        [00:19:32 -19541.015126] FAST spr round 3 (radius: 5)
        [00:19:36 -19541.014612] Model parameter optimization (eps = 1.000000)
        [00:19:36 -19540.806426] SLOW spr round 1 (radius: 5)
        [00:19:46 -19540.396377] SLOW spr round 2 (radius: 5)
        [00:19:56 -19540.392885] SLOW spr round 3 (radius: 10)
        [00:20:04 -19540.392257] SLOW spr round 4 (radius: 15)
        [00:20:09 -19540.392144] SLOW spr round 5 (radius: 20)
        [00:20:09 -19540.392123] SLOW spr round 6 (radius: 25)
        [00:20:10 -19540.392120] Model parameter optimization (eps = 0.100000)

        [00:20:10] ML tree search #19, logLikelihood: -19540.353600

        [00:20:10 -21480.318201] Initial branch length optimization
        [00:20:10 -19636.576063] Model parameter optimization (eps = 10.000000)
        [00:20:10 -19607.608647] AUTODETECT spr round 1 (radius: 5)
        [00:20:13 -19541.992968] AUTODETECT spr round 2 (radius: 10)
        [00:20:16 -19541.987898] SPR radius for FAST iterations: 5 (autodetect)
        [00:20:16 -19541.987898] Model parameter optimization (eps = 3.000000)
        [00:20:16 -19540.327022] FAST spr round 1 (radius: 5)
        [00:20:20 -19539.715583] FAST spr round 2 (radius: 5)
        [00:20:24 -19539.329928] FAST spr round 3 (radius: 5)
        [00:20:28 -19539.329018] Model parameter optimization (eps = 1.000000)
        [00:20:28 -19539.235162] SLOW spr round 1 (radius: 5)
        [00:20:38 -19539.230917] SLOW spr round 2 (radius: 10)
        [00:20:46 -19539.230806] SLOW spr round 3 (radius: 15)
        [00:20:52 -19539.230796] SLOW spr round 4 (radius: 20)
        [00:20:54 -19539.230795] SLOW spr round 5 (radius: 25)
        [00:20:55 -19539.230794] Model parameter optimization (eps = 0.100000)

        [00:20:55] ML tree search #20, logLikelihood: -19539.222087


        Optimized model parameters:

        Partition 0: noname
        Rate heterogeneity: GAMMA (4 cats, mean),  alpha: 1.911585 (ML),  weights&rates: (0.250000,0.282487) (0.250000,0.645270) (0.250000,1.067611) (0.250000,2.004632) 
        Base frequencies (model): 0.076748 0.051691 0.042645 0.051544 0.019803 0.040752 0.061830 0.073152 0.022944 0.053761 0.091904 0.058676 0.023826 0.040126 0.050901 0.068765 0.058565 0.014261 0.032102 0.066005 
        Substitution rates (model): 58.000000 54.000000 81.000000 56.000000 57.000000 105.000000 179.000000 27.000000 36.000000 30.000000 35.000000 54.000000 15.000000 194.000000 378.000000 475.000000 9.000000 11.000000 298.000000 45.000000 16.000000 113.000000 310.000000 29.000000 137.000000 328.000000 22.000000 38.000000 646.000000 44.000000 5.000000 74.000000 101.000000 64.000000 126.000000 20.000000 17.000000 528.000000 34.000000 86.000000 58.000000 81.000000 391.000000 47.000000 12.000000 263.000000 30.000000 10.000000 15.000000 503.000000 232.000000 8.000000 70.000000 16.000000 10.000000 49.000000 767.000000 130.000000 112.000000 11.000000 7.000000 26.000000 15.000000 4.000000 15.000000 59.000000 38.000000 4.000000 46.000000 31.000000 9.000000 5.000000 59.000000 69.000000 17.000000 23.000000 7.000000 31.000000 78.000000 14.000000 223.000000 42.000000 115.000000 209.000000 62.000000 323.000000 26.000000 597.000000 9.000000 72.000000 292.000000 43.000000 4.000000 164.000000 53.000000 51.000000 18.000000 24.000000 20.000000 119.000000 26.000000 12.000000 9.000000 181.000000 18.000000 5.000000 18.000000 30.000000 32.000000 10.000000 7.000000 45.000000 23.000000 6.000000 6.000000 27.000000 14.000000 5.000000 24.000000 201.000000 33.000000 55.000000 8.000000 47.000000 16.000000 56.000000 45.000000 33.000000 40.000000 115.000000 73.000000 46.000000 8.000000 573.000000 11.000000 229.000000 21.000000 479.000000 89.000000 10.000000 40.000000 245.000000 9.000000 32.000000 961.000000 14.000000 388.000000 248.000000 102.000000 59.000000 25.000000 52.000000 24.000000 180.000000 65.000000 4.000000 21.000000 47.000000 103.000000 10.000000 8.000000 14.000000 43.000000 16.000000 29.000000 226.000000 24.000000 18.000000 323.000000 17.000000 92.000000 12.000000 53.000000 536.000000 62.000000 285.000000 118.000000 6.000000 10.000000 23.000000 477.000000 35.000000 63.000000 38.000000 12.000000 21.000000 112.000000 71.000000 25.000000 16.000000 


        Final LogLikelihood: -19539.221227

        AIC score: 39262.442454 / AICc score: 39279.169726 / BIC score: 39724.053019
        Free parameters (model + branch lengths): 92

        Best ML tree saved to: /Users/melettedevore/Desktop/PP563/project/alignments/T3.raxml.bestTree
        All ML trees saved to: /Users/melettedevore/Desktop/PP563/project/alignments/T3.raxml.mlTrees
        Optimized model saved to: /Users/melettedevore/Desktop/PP563/project/alignments/T3.raxml.bestModel

        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/alignments/T3.raxml.log

        Analysis started: 27-Mar-2025 11:58:44 / finished: 27-Mar-2025 12:19:39

        Elapsed time: 1255.189 seconds

03 May 2025
    After looking through my sequences and data again, I don't think the sequences I have are answering the questions I have. In addition, I have too many taxa! So, I am starting over with fewer taxa and pulling better sequences from NCBI. I have 5 bacterial expansins, 5 fungal expansins, and 3 nematode expansins. My outgroup gene is an expansin from a virus.

    I made a file with the original names from NCBI and a file with the names I want to appear on the alignments and trees. The original names file is "AllSeqs_OrigNames" and the preferred names file is "AllSeqs", the files are stored in Desktop/PP563/Project.

   INFERENCE METHOD #1: DISTANCE-BASED
    T-Coffee Alignment:
       
        (base) Melettes-MBP-6:project melettedevore$ t_coffee AllSeqs.fasta

        PROGRAM: T-COFFEE Version_12.00.7fb08c2 (2018-12-11 09:27:12 - Revision 7fb08c2 - Build 211)
        -full_log      	S	[0] 
        -genepred_score	S	[0] 	nsd
        -run_name      	S	[0] 
        -mem_mode      	S	[0] 	mem
        -extend        	D	[1] 	1 
        -extend_mode   	S	[0] 	very_fast_triplet
        -max_n_pair    	D	[0] 	10 
        -seq_name_for_quadruplet	S	[0] 	all
        -compact       	S	[0] 	default
        -clean         	S	[0] 	no
        -do_self       	FL	[0] 	0
        -do_normalise  	D	[0] 	1000 
        -template_file 	S	[0] 
        -setenv        	S	[0] 	0
        -export        	S	[0] 	0
        -template_mode 	S	[0] 
        -flip          	D	[0] 	0 
        -remove_template_file	D	[0] 	0 
        -profile_template_file	S	[0] 
        -in            	S	[0] 
        -seq           	S	[1] 	AllSeqs.fasta
        -aln           	S	[0] 
        -method_limits 	S	[0] 
        -method        	S	[0] 
        -lib           	S	[0] 
        -profile       	S	[0] 
        -profile1      	S	[0] 
        -profile2      	S	[0] 
        -pdb           	S	[0] 
        -relax_lib     	D	[0] 	1 
        -filter_lib    	D	[0] 	0 
        -shrink_lib    	D	[0] 	0 
        -out_lib       	W_F	[0] 	no
        -out_lib_mode  	S	[0] 	primary
        -lib_only      	D	[0] 	0 
        -outseqweight  	W_F	[0] 	no
        -seq_source    	S	[0] 	ANY
        -cosmetic_penalty	D	[0] 	0 
        -gapopen       	D	[0] 	0 
        -gapext        	D	[0] 	0 
        -fgapopen      	D	[0] 	0 
        -fgapext       	D	[0] 	0 
        -nomatch       	D	[0] 	0 
        -newtree       	W_F	[0] 	default
        -tree          	W_F	[0] 	NO
        -usetree       	R_F	[0] 
        -tree_mode     	S	[0] 	nj
        -distance_matrix_mode	S	[0] 	ktup
        -distance_matrix_sim_mode	S	[0] 	idmat_sim1
        -quicktree     	FL	[0] 	0
        -outfile       	W_F	[0] 	default
        -maximise      	FL	[1] 	1
        -output        	S	[0] 	aln	html
        -len           	D	[0] 	0 
        -infile        	R_F	[0] 
        -matrix        	S	[0] 	default
        -tg_mode       	D	[0] 	1 
        -profile_mode  	S	[0] 	cw_profile_profile
        -profile_comparison	S	[0] 	profile
        -dp_mode       	S	[0] 	linked_pair_wise
        -ktuple        	D	[0] 	1 
        -ndiag         	D	[0] 	0 
        -diag_threshold	D	[0] 	0 
        -diag_mode     	D	[0] 	0 
        -sim_matrix    	S	[0] 	vasiliky
        -transform     	S	[0] 
        -extend_seq    	FL	[0] 	0
        -outorder      	S	[0] 	input
        -inorder       	S	[0] 	aligned
        -seqnos        	S	[0] 	off
        -case          	S	[0] 	keep
        -cpu           	D	[0] 	0 
        -ulimit        	D	[0] 	-1 
        -maxnseq       	D	[0] 	-1 
        -maxlen        	D	[0] 	-1 
        -sample_dp     	D	[0] 	0 
        -weight        	S	[0] 	default
        -seq_weight    	S	[0] 	no
        -align         	FL	[1] 	1
        -mocca         	FL	[0] 	0
        -domain        	FL	[0] 	0
        -start         	D	[0] 	0 
        -len           	D	[0] 	0 
        -scale         	D	[0] 	0 
        -mocca_interactive	FL	[0] 	0
        -method_evaluate_mode	S	[0] 	default
        -color_mode    	S	[0] 	new
        -aln_line_length	D	[0] 	0 
        -evaluate_mode 	S	[0] 	triplet
        -get_type      	FL	[0] 	0
        -clean_aln     	D	[0] 	0 
        -clean_threshold	D	[1] 	1 
        -clean_iteration	D	[1] 	1 
        -clean_evaluate_mode	S	[0] 	t_coffee_fast
        -extend_matrix 	FL	[0] 	0
        -prot_min_sim  	D	[40] 	40 
        -prot_max_sim  	D	[90] 	90 
        -prot_trim     	D	[20] 	20 
        -prot_min_cov  	D	[40] 	40 
        -pdb_type      	S	[0] 	d
        -pdb_min_sim   	D	[35] 	35 
        -pdb_max_sim   	D	[100] 	100 
        -pdb_min_cov   	D	[50] 	50 
        -pdb_blast_server	W_F	[0] 	EBI
        -blast         	W_F	[0] 
        -blast_server  	W_F	[0] 	EBI
        -pdb_db        	W_F	[0] 	pdb
        -protein_db    	W_F	[0] 	uniref50
        -method_log    	W_F	[0] 	no
        -struc_to_use  	S	[0] 
        -cache         	W_F	[0] 	use
        -print_cache   	FL	[0] 	0
        -align_pdb_param_file	W_F	[0] 	no
        -align_pdb_hasch_mode	W_F	[0] 	hasch_ca_trace_bubble
        -external_aligner	S	[0] 	NO
        -msa_mode      	S	[0] 	tree
        -et_mode       	S	[0] 	et
        -master        	S	[0] 	no
        -blast_nseq    	D	[0] 	0 
        -lalign_n_top  	D	[0] 	10 
        -iterate       	D	[0] 	0 
        -trim          	D	[0] 	0 
        -split         	D	[0] 	0 
        -trimfile      	S	[0] 	default
        -split         	D	[0] 	0 
        -split_nseq_thres	D	[0] 	0 
        -split_score_thres	D	[0] 	0 
        -check_pdb_status	D	[0] 	0 
        -clean_seq_name	D	[0] 	0 
        -seq_to_keep   	S	[0] 
        -dpa_master_aln	S	[0] 
        -dpa_maxnseq   	D	[0] 	0 
        -dpa_min_score1	D	[0] 
        -dpa_min_score2	D	[0] 
        -dpa_keep_tmpfile	FL	[0] 	0
        -dpa_debug     	D	[0] 	0 
        -multi_core    	S	[0] 	templates_jobs_relax_msa_evaluate
        -n_core        	D	[0] 	0 
        -max_n_proc    	D	[0] 	0 
        -lib_list      	S	[0] 
        -prune_lib_mode	S	[0] 	5
        -tip           	S	[0] 	none
        -rna_lib       	S	[0] 
        -no_warning    	D	[0] 	0 
        -run_local_script	D	[0] 	0 
        -plugins       	S	[0] 	default
        -proxy         	S	[0] 	unset
        -email         	S	[0] 
        -clean_overaln 	D	[0] 	0 
        -overaln_param 	S	[0] 
        -overaln_mode  	S	[0] 
        -overaln_model 	S	[0] 
        -overaln_threshold	D	[0] 	0 
        -overaln_target	D	[0] 	0 
        -overaln_P1    	D	[0] 	0 
        -overaln_P2    	D	[0] 	0 
        -overaln_P3    	D	[0] 	0 
        -overaln_P4    	D	[0] 	0 
        -exon_boundaries	S	[0] 
        -dump          	S	[0] 	no
        -display       	D	[0] 	100 

        INPUT FILES
            Input File (S) AllSeqs.fasta  Format fasta_seq
            Input File (M) proba_pair 

        Identify Master Sequences [no]:

        Master Sequences Identified
        INPUT SEQUENCES: 14 SEQUENCES  [PROTEIN]
        Input File AllSeqs.fasta                   Seq Abortiporus_biennis             Length  142 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Bacillus_subtilis               Length  232 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Bursaphelenchus_xylophilus_EXP1 Length  149 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Colletotrichum_destructivum     Length  219 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Faustovirus                     Length  253 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Gigaspora_rosea                 Length  210 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Globodera_rostochiensis_EXPB1   Length  271 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Meloidogyne_javanica_AP1        Length  289 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Rhizoctonia_solani              Length  225 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Stemphylium_lycopersici         Length  306 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Streptomyces_mirabilis          Length  315 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Xanthomonas_campestris          Length  590 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Xanthomonas_oryzae              Length  590 type PROTEIN Struct Unchecked
        Input File AllSeqs.fasta                   Seq Xylella_fastidiosa              Length  329 type PROTEIN Struct Unchecked

            Multi Core Mode: 4 processors:

            --- Process Method/Library/Aln SAllSeqs.fasta
            --- Process Method/Library/Aln Mproba_pair
            xxx Retrieved SAllSeqs.fasta
            xxx Retrieved Mproba_pair

            All Methods Retrieved

        MANUAL PENALTIES: gapopen=0 gapext=0

            Library Total Size: [151178]

        Library Relaxation: Multi_proc [4]
        
        !		[Relax Library][TOT=    3][100 %]

        Relaxation Summary: [151178]--->[61436]



        UN-WEIGHTED MODE: EVERY SEQUENCE WEIGHTS 1

        MAKE GUIDE TREE 
            [MODE=nj][DONE]

        PROGRESSIVE_ALIGNMENT [Tree Based]
        Group    1: Abortiporus_biennis
        Group    2: Bacillus_subtilis
        Group    3: Bursaphelenchus_xylophilus_EXP1
        Group    4: Colletotrichum_destructivum
        Group    5: Faustovirus
        Group    6: Gigaspora_rosea
        Group    7: Globodera_rostochiensis_EXPB1
        Group    8: Meloidogyne_javanica_AP1
        Group    9: Rhizoctonia_solani
        Group   10: Stemphylium_lycopersici
        Group   11: Streptomyces_mirabilis
        Group   12: Xanthomonas_campestris
        Group   13: Xanthomonas_oryzae
        Group   14: Xylella_fastidiosa


            Group   15: [Group    6 (   1 seq)] with [Group    1 (   1 seq)]-->[Len=  211][PID:16421]
            Group   16: [Group    5 (   1 seq)] with [Group   15 (   2 seq)]-->[Len=  309][PID:16421]
            Group   17: [Group    9 (   1 seq)] with [Group    8 (   1 seq)]-->[Len=  327][PID:16421]
            Group   18: [Group   11 (   1 seq)] with [Group   10 (   1 seq)]-->[Len=  425][PID:16421]
            Group   19: [Group   18 (   2 seq)] with [Group   17 (   2 seq)]-->[Len=  513][PID:16421]
            Group   20: [Group    4 (   1 seq)] with [Group   19 (   4 seq)]-->[Len=  523][PID:16421]
            Group   21: [Group   13 (   1 seq)] with [Group   12 (   1 seq)]-->[Len=  590][PID:16421]
            Group   22: [Group   21 (   2 seq)] with [Group    2 (   1 seq)]-->[Len=  591][PID:16421]
            Group   23: [Group   14 (   1 seq)] with [Group   22 (   3 seq)]-->[Len=  604][PID:16421]
            Group   24: [Group   23 (   4 seq)] with [Group   20 (   5 seq)]-->[Len=  765][PID:16421]
            Group   25: [Group    7 (   1 seq)] with [Group   24 (   9 seq)]-->[Len=  787][PID:16421]
            Group   26: [Group   25 (  10 seq)] with [Group   16 (   3 seq)]-->[Len=  850][PID:16421]
            Group   27: [Group    3 (   1 seq)] with [Group   26 (  13 seq)]-->[Len=  852][PID:16421]


        !		[Final Evaluation][TOT=  213][100 %]

        CLUSTAL FORMAT for T-COFFEE Version_12.00.7fb08c2 [http://www.tcoffee.org] [MODE:  ], CPU=0.00 sec, SCORE=672, Nseq=14, Len=852 

        Xanthomonas_campestris           MNPRAQWLVACVALAPVVACAMQVDAHGQLRDTSGHAVQLRGVTWPGFDS
        Bacillus_subtilis                MKKIMSAFVGIVLLT-IFCFSPQ---------------------------
        Xanthomonas_oryzae               MLSHIHFLIALATLVPVTAPAMQVSTQAPLVDATGQTLHIRGVTWPGFDR
        Xylella_fastidiosa               MPRFAGTD------------------------------------------
        Streptomyces_mirabilis           MKAASQA------------ARR----------------------------
        Abortiporus_biennis              MVPSHISLL-----------------------------------------
        Colletotrichum_destructivum      MHGLLSIGTALA-------VMAQ---------------------G-----
        Gigaspora_rosea                  MFSSITRIG-------------------------------NSSQSSSEKS
        Rhizoctonia_solani               MFSTAVAAFVLAA---AVNQAAA---------------------E-----
        Stemphylium_lycopersici          MKTSAVLASLLFG---SYAVAA----------------------------
        Globodera_rostochiensis_EXPB1    MSSSEAILCLLCL----LAVNFR----AQIVLAS-VTAKLEGKSWNGGGQ
        Bursaphelenchus_xylophilus_EXP1  MNSFYLLA-LLVP----LALSDN---------------------------
        Meloidogyne_javanica_AP1         MASLKLFITLMS----AFILSKT---------------------F-----
        Faustovirus                      MAAAANDNTT------------------------------SGLAVPSYDS
                                        *                                                 

        Xanthomonas_campestris           ASHVARGLQDITLQQALDQMQAAQINALRIPVCPASLQAQPLAARDVAGE
        Bacillus_subtilis                --------------------------------------------------
        Xanthomonas_oryzae               AGLAAAGMRNNTLAQLLDRMQASDINAVRVPVCAAVLQRAPVAAAEVAGD
        Xylella_fastidiosa               --------------------------------------------------
        Streptomyces_mirabilis           --------------------------------------------------
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      --------------------------------------------------
        Gigaspora_rosea                  L-------------------------------------------------
        Rhizoctonia_solani               --------------------------------------------------
        Stemphylium_lycopersici          --------------------------------------------------
        Globodera_rostochiensis_EXPB1    ---YVPNFKNNDGSKIA---------------CSVKFSLTPK--------
        Bursaphelenchus_xylophilus_EXP1  --------------------------------------I-----------
        Meloidogyne_javanica_AP1         --------------------------------------------------
        Faustovirus                      C-------------------------------------------------
                                                                                        

        Xanthomonas_campestris           PALRGATALQALDAVVRASTQRGMQVMLAFANGGCDDQAPQLGADTTTWT
        Bacillus_subtilis                --------------------------------------------------
        Xanthomonas_oryzae               STLRGLDSLQLLDAVVHAATQRGMQVVFAFADGGCDDRAPLLGAQQQAWT
        Xylella_fastidiosa               --------------------------------------------------
        Streptomyces_mirabilis           ----------------------------------------------RPRR
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      --------------------------------------------------
        Gigaspora_rosea                  --------------------------------------------------
        Rhizoctonia_solani               ----------------------------------------------PSRA
        Stemphylium_lycopersici          ----------------------------------------------PFDR
        Globodera_rostochiensis_EXPB1    --------------------------------------------------
        Bursaphelenchus_xylophilus_EXP1  --------------------------------------------------
        Meloidogyne_javanica_AP1         --------------------------------------------------
        Faustovirus                      --------------------------------------------------
                                                                                        

        Xanthomonas_campestris           RGLRMLAQRYRTNTGVIGV----DL-----GSSGYRNATWAGNPP--A-S
        Bacillus_subtilis                --------------------------------------------------
        Xanthomonas_oryzae               RGLVTLARRYGGNTNVLGI----DL-----GSSGYRNASWAGNAA--D-Q
        Xylella_fastidiosa               --------------------------------------------------
        Streptomyces_mirabilis           RGLVT-GTALAMVATVLLVSLIVVL-RS-GRDAGTGN------A------
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      ------------------------------T-------------------
        Gigaspora_rosea                  --PAL-PRPVQSKKSNNTS-K-----S-----TGWQKTKT----------
        Rhizoctonia_solani               RHLYA-ARAYAQDSHLLED-YD--TYHARYIAIGC------KDAKKNDTD
        Stemphylium_lycopersici          RALVY-KTEVVTETVVIYT----TV-Y----DDEYPEA-----T------
        Globodera_rostochiensis_EXPB1    ----------------KGTT----I-----------GSVWGANAV--S-G
        Bursaphelenchus_xylophilus_EXP1  --------------------------------------------------
        Meloidogyne_javanica_AP1         ------------------------------TSARCPKGTHPKNGE--KS-
        Faustovirus                      --VEF-MRIAQENKHNITN-A-----------------------------
                                                                                        

        Xanthomonas_campestris           DWNRVATR----------------------AAAAVLQHAPHWVIGVEGVG
        Bacillus_subtilis                --------------------------------------------------
        Xanthomonas_oryzae               NWNRVASR----------------------AAAMVLAQAPRWVVGVEGVG
        Xylella_fastidiosa               --------------------------------------------------
        Streptomyces_mirabilis           ----A---------------------------------------------
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      --------------------------------------------------
        Gigaspora_rosea                  --------------------------------------------------
        Rhizoctonia_solani               FWNNC---------------------------------------------
        Stemphylium_lycopersici          ----S---------------------------------------------
        Globodera_rostochiensis_EXPB1    ASNQYTLAPPADIAPGATHTNAGVNINGNGAPTLKLIEAKYFIDDVCGG-
        Bursaphelenchus_xylophilus_EXP1  --------------------------------------------------
        Meloidogyne_javanica_AP1         TINPI---------------------------------------------
        Faustovirus                      --------------------------------------------------
                                                                                        

        Xanthomonas_campestris           SNAVCSDTARKAPGSNLQPLACVPPRI-AAANL-----VLMPKLAGPDRD
        Bacillus_subtilis                --------------------------------------------------
        Xanthomonas_oryzae               SNAVCSDPGRKALGSNLQPFACVPLDI-DRRHL-----VLMPKLAGPDRD
        Xylella_fastidiosa               ------------------------------------------------MD
        Streptomyces_mirabilis           ---------------------ATPVAG-SQETA-----------------
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      --------------------------------------------------
        Gigaspora_rosea                  --------------TC-----LTPYYK-WSEWL-----EKKHPKVAPYKF
        Rhizoctonia_solani               ---------------------CHPL--------------AKDADS-----
        Stemphylium_lycopersici          ---------------------SAPVYE-QQKPSTTSTAVVYPTSS-----
        Globodera_rostochiensis_EXPB1    --------------------------------------------------
        Bursaphelenchus_xylophilus_EXP1  --------------------------------------------------
        Meloidogyne_javanica_AP1         ---------------------NNPPAVGTIHPP------VQPPSN-----
        Faustovirus                      ---------------------ITPNTG-TERWF-----IEFCKTAKPDQ-
                                                                                        

        Xanthomonas_campestris           A-----------------AEAFAARDFPDALPATWQRDFGTLAAT--HAV
        Bacillus_subtilis                --------------------------------------------------
        Xanthomonas_oryzae               T-----------------TDAFAAPGFAQALPAMWQRDFGQFAID--HTV
        Xylella_fastidiosa               T-----------------EAA-----------------LGAFSGK--QTV
        Streptomyces_mirabilis           ------------------------------VPPV----------------
        Abortiporus_biennis              -APIFALIFLLSMLLPTVH--AAPISIQPLQKRASL--------------
        Colletotrichum_destructivum      --------------------------------------------------
        Gigaspora_rosea                  AILISLGLFLLAILIIIIV--AAVGGFQNNDTGTTG--------------
        Rhizoctonia_solani               -----------------------------SVPAQCKA-------------
        Stemphylium_lycopersici          -----------------------------AVPAP----------------
        Globodera_rostochiensis_EXPB1    --------------------------------------------------
        Bursaphelenchus_xylophilus_EXP1  --------------------------------------------------
        Meloidogyne_javanica_AP1         -----------------------------KPPSPGYPS------------
        Faustovirus                      -----------------IV--HLLAGFTDDLPATTLA------AETVNKI
                                                                                        

        Xanthomonas_campestris           LPVSTGGGLGDGDPRDP----D-WQHALSTYLDTTG-------LRSAFLG
        Bacillus_subtilis                --------------------------------------------------
        Xanthomonas_oryzae               VPVSLGGGLGDGDPRDP----A-WQTALSGYLANAG-------MRSAFLG
        Xylella_fastidiosa               LPNSLDA-------TDA----EQLAHRIDA-LLAFG-------IRQGFYG
        Streptomyces_mirabilis           ---ASGHKPPTASPAAT-------------------------------AG
        Abortiporus_biennis              -----------------NDAH-----------------------------
        Colletotrichum_destructivum      --------------------------------------------------
        Gigaspora_rosea                  -----------------LAMS-----------------------------
        Rhizoctonia_solani               ---------SSKSYS----------------------------------G
        Stemphylium_lycopersici          ----------ASTSSAA-------------------------------PA
        Globodera_rostochiensis_EXPB1    -------------------------------AP-AG-------SCMGC--
        Bursaphelenchus_xylophilus_EXP1  --------------------------------------------------
        Meloidogyne_javanica_AP1         --PPSSNLPKQPSPP----------------------------------G
        Faustovirus                      LPYH------DSV---SKDAYK-YVHALCPYVHDRWNLNKDKCIGGSFGG
                                                                                        

        Xanthomonas_campestris           SWEIGDANNGGLLARD-------------------------------G-T
        Bacillus_subtilis                --------------------------------------------------
        Xanthomonas_oryzae               SWETGNANNGGLLAPD-------------------------------G-S
        Xylella_fastidiosa               SWMTSAQMPFGLLDND-------------------------------GRT
        Streptomyces_mirabilis           ATTPVP---SAATAPA-------------------------------R-T
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      --------------------------------------------------
        Gigaspora_rosea                  --------------------------------------------------
        Rhizoctonia_solani               SGDSSPSSSSAAYTPA-------------------------------P-T
        Stemphylium_lycopersici          YTPVVP---EKPSTPAYTPVAPESSTPAAETSTSIYTPEPVYTPKPSP-T
        Globodera_rostochiensis_EXPB1    ----LSNTKMDGP--I---------------------------------N
        Bursaphelenchus_xylophilus_EXP1  -------------------------------------------------T
        Meloidogyne_javanica_AP1         KYPSGPPNHNQPKQPS-------------------------------P-T
        Faustovirus                      LVDCYWRNMGDLTYKP----------------------------------
                                                                                        

        Xanthomonas_campestris           P-----------------RADKLQVLHRAWGVSTALPATASA--KRASAA
        Bacillus_subtilis                -------------------A------------------------------
        Xanthomonas_oryzae               P-----------------RADKLLILRHAWGMLPVMPAIATA--TGDSTK
        Xylella_fastidiosa               P-----------------RTSLIAQLHRWWGVSRVDVASENATTKNQTTT
        Streptomyces_mirabilis           P-----------------------------------PATTRQ--PASGTA
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      --------------------------------------------------
        Gigaspora_rosea                  --------------------------------------------------
        Rhizoctonia_solani               P-----------------------------------DSG-----------
        Stemphylium_lycopersici          PSPEPEPETSSVVPPPAPTT-------TEAAYTPVTPSLQQA--AAAPTA
        Globodera_rostochiensis_EXPB1    K-----------------NL------------------------------
        Bursaphelenchus_xylophilus_EXP1  P-----------------KAN-----------------------------
        Meloidogyne_javanica_AP1         P-----------------------------------PSSGPS--P-----
        Faustovirus                      ------------------QLN-----------------------------
                                                                                        

        Xanthomonas_campestris           DA---------SKATGWDS-------------T-----F---------SG
        Bacillus_subtilis                -------------SAAYDD-------------L-----H---------EG
        Xanthomonas_oryzae               NA---------SGTKGWNS-------------T-----F---------TG
        Xylella_fastidiosa               DTNGCVAGDNSVPLNGWDT-------------S-----F---------SG
        Streptomyces_mirabilis           SL---------AGRVRPET-------------T-----Y---------RG
        Abortiporus_biennis              --------------------------------------T---------GG
        Colletotrichum_destructivum      -----------------SA-------------A-----T---------SA
        Gigaspora_rosea                  --------------------------------------G---------SG
        Rhizoctonia_solani               ----------------SGE------------------TYT--------GG
        Stemphylium_lycopersici          --------------AAFDP-------------AEPGYTSTGVSTGSFTDV
        Globodera_rostochiensis_EXPB1    -----------------NK-------------P-----FK--------NS
        Bursaphelenchus_xylophilus_EXP1  -----------------KP-------------I-----K---------GG
        Meloidogyne_javanica_AP1         ----------------PSS-------------APS--VPA--------LP
        Faustovirus                      ------------------EFAKLVLSSLTSGFT-----F---------SI
                                                                                        

        Xanthomonas_campestris           TATYTGSGY-----------------------------------------
        Bacillus_subtilis                YATYTGSGY-----------------------------------------
        Xanthomonas_oryzae               IATPTGSGY-----------------------------------------
        Xylella_fastidiosa               VATYTYTGY-----------------------------------------
        Streptomyces_mirabilis           TATHYDAGT-----------------------------------------
        Abortiporus_biennis              RGTFFFP-------------------------------------------
        Colletotrichum_destructivum      IASWYGGNL-----------------------------------------
        Gigaspora_rosea                  DGTYYDPGV-----------------------------------------
        Rhizoctonia_solani               SGTYFYQNG-----------------------------------------
        Stemphylium_lycopersici          DMTVYDNAG-----------------------------------------
        Globodera_rostochiensis_EXPB1    VFTFYGA-G-----------------------------------------
        Bursaphelenchus_xylophilus_EXP1  EFTYYNA-V-----------------------------------------
        Meloidogyne_javanica_AP1         SGQYPPGSNVTSNSLTVPPLYTTNETCLPLTPDMPGVVDQWLNKPKTGRF
        Faustovirus                      TQSKMSNDS-----------------------------------------
                                                                                        

        Xanthomonas_campestris           ----------SGGALLLD---PIPRTAFITALNPVQLNFGG---------
        Bacillus_subtilis                ----------SGGAFLLD---PIPSDMEITAINPADLNYGG---------
        Xanthomonas_oryzae               ----------SGGAVLLD---PIPSDAFITALNPTQMXFGG---------
        Xylella_fastidiosa               ----------KGGALMLD---PIQSHVQITALNPTQLNLGG---------
        Streptomyces_mirabilis           ----------GDGACLYG---PS-DDRMTAAMNHTDYE------------
        Abortiporus_biennis              ----------GLGACGKT---NG-DNDMIVAVNQAVFQQF----------
        Colletotrichum_destructivum      ----------NGGNCALTGYTLP-SGVLGTALSYTNYN------------
        Gigaspora_rosea                  ----------GLGSCGWQ---NY-DSELVAAMNAPQYGVFA---------
        Rhizoctonia_solani               ----------VAGACGTV---HS-DSDYIVAVDYRRYGDLG---------
        Stemphylium_lycopersici          ----------TPGACGVK---LT-DDMMIAAIAQPAWDAKGGSTYDTA--
        Globodera_rostochiensis_EXPB1    ----------GRGACGLD---AG-VPKMSAAGSGNLFKPDGQWVDACRKD
        Bursaphelenchus_xylophilus_EXP1  ----------GDGACGRT---LT---GCVAAVAPTLFDPSAKWVASDLPD
        Meloidogyne_javanica_AP1         VWHALTNAHSYGGECGLG---IQ-YP-YHACVGSDLWSPVKKWVRSCLTD
        Faustovirus                      ----------GLKCCVVD---CY-DDAIRYTGSY-----Y----------
                                                                                        

        Xanthomonas_campestris           -------VKAALA-G-A------------------YLQVTG-PKG-TTTV
        Bacillus_subtilis                -------VIAALA-G-S------------------YLEVEG-PKG-KTTV
        Xanthomonas_oryzae               -------IKAALA-G-A------------------YLEVTG-PKG-KTTV
        Xylella_fastidiosa               -------IPAAMA-G-A------------------YLRVQG-PKG-STTV
        Streptomyces_mirabilis           --------SAKAC-G-A------------------YVRVRA-AGGASVTV
        Abortiporus_biennis              --------RDVICDN-K------------------SMEVRA-NGK-SATV
        Colletotrichum_destructivum      ----------GNC-G-T------------------CLNVKG-PKGNTIKV
        Gigaspora_rosea                  ----DP-GNSPVC-G-K------------------CIQVTG-PKG-TVKV
        Rhizoctonia_solani               -------KQSDLC-G-K------------------KVLVTNTSNGKSVTC
        Stemphylium_lycopersici          ---TGK-SSNPWC-G-T------------------EIDVHH--NGETVKV
        Globodera_rostochiensis_EXPB1    K--RTL-LDDPIC-KNI------------------CVKIDY--NGKTLTV
        Bursaphelenchus_xylophilus_EXP1  K--RYVLGDKVCQ-G-M------------------CVKVEY--KGKTASF
        Meloidogyne_javanica_AP1         IPNHYI-RNDQVC-VNK------------------CVKI--EYKNKTLIM
        Faustovirus                      --------NGYIC-T-NCRSAYHTSCFNKINNETLDAKMDF-NG------
                                                                            :           

        Xanthomonas_campestris           YVTDLYPE-G-ASG-GLDLSHNAFAAIGDM-VQGRIPI-SWKVVRAPVTG
        Bacillus_subtilis                YVTDLYPE-G-AQG-ALDLSPNAFRKICNM-KDGKINI-KWRVVKAPITG
        Xanthomonas_oryzae               YVTDLYPE-A-ASG-GLDLSYNAFAAIGNM-ADGDIPI-SWKVVRAPITG
        Xylella_fastidiosa               YVTDLYPT-G-SSG-GLDLSPNAFASIGNM-AQGRIPV-QWKVVSAPVSG
        Streptomyces_mirabilis           RITNECPLPC-APG-QLDLSAQAFAVLAAP-SLGRIPV-TWSLLSPSTSD
        Abortiporus_biennis              QIVDECPG-C-GPN-DIDLSPAAFKVFADE-SVGVLTV-DWNIF------
        Colletotrichum_destructivum      MVVDKCPEGC-GAG-QLDLFPNVFAALDNP-DKGLINV-QWEQVPCGITS
        Gigaspora_rosea                  KVVDKCPV-C-KSG-DIDMSSTAFKQIANL-DDGRVKI-TWKGC------
        Rhizoctonia_solani               TVADACPT-CDNEN-SLDMSEGAFKQLGDL-SLGLLKI-KYTVL------
        Stemphylium_lycopersici          TIMDLCPG-C-VGH-DIDLSHAAWKGLGLG-APDRFKA-SWNVV------
        Globodera_rostochiensis_EXPB1    PINNKCPE-C-TPS-HVDLSIDAFNYLEPR--GGLVG--KATGX------
        Bursaphelenchus_xylophilus_EXP1  PIEDKCPG-C-QEN-HVDLSENAFKILEXNLGVGKAKDATITYV------
        Meloidogyne_javanica_AP1         PIMDETDY-LK-PY-ELDISEPAFKYLEMNPNAGWCTA-TVTYI------
        Faustovirus                      QRINMCPV-C-HITWDLKNAEAHDKGY----------Q-K---W------
                                            :            :.                                

        Xanthomonas_campestris           NLQYRIKEGGSRWWAAIQVRNHAYPVVKLEVKQGS--TWKNLQKMDYNHF
        Bacillus_subtilis                NFTYRIKEGSSRWWAAIQVRNHKYPVMKMEYEKDG--KWINMEKMDYNHF
        Xanthomonas_oryzae               NVQYRIKEGSSRYWAAIQVRHHIYPVVKLEVKQGS--TWTSLPKTDYNHF
        Xylella_fastidiosa               NLIYRVKKGSSGWWAAIQVREHRYPVLKLEICQDG--TWLNLPKRNYNYF
        Streptomyces_mirabilis           TVSIRYKTGSSRWWCAIQVIGHRNPVARLEVRTGG--GWHQLPRTDYNYF
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      PLVVRNKEGTSKWWFSMQVMNHNYPVTKFEVSTDGGKSWQPTVRQDYNYW
        Gigaspora_rosea                  --------------------------------------------------
        Rhizoctonia_solani               --------------------------------------------------
        Stemphylium_lycopersici          --------------------------------------------------
        Globodera_rostochiensis_EXPB1    ----------------------RSPI------------------------
        Bursaphelenchus_xylophilus_EXP1  --------------------------------------------------
        Meloidogyne_javanica_AP1         --------------------------------------------------
        Faustovirus                      --------------------------------------------------
                                                                                        

        Xanthomonas_campestris           LGEQ---LGDQPLTLRITDIRGKVLTDTLPRLPEDGSKPAYFEPGHVQFP
        Bacillus_subtilis                VSTN---LGTGSLKVRMTDIRGKVVKDTIPKLPESGTSKAYTVPGHVQFP
        Xanthomonas_oryzae               VGTN---LGNKPLSIRITDIRGKVIADKIPALPEYGGSKAYFEPGHVQFP
        Xylella_fastidiosa               VGTR---LGNQPLSMRMTDIRGQTLIDTLPALPNKASSKAYSVNGNVQFS
        Streptomyces_mirabilis           LSDQ---GSGCGGAIRITDIYGEPLVLNGIALRPD-----AVQPTRVQFA
        Abortiporus_biennis              --------------------------------------------------
        Colletotrichum_destructivum      QRSSGDGFNVDTVTVRVTCSNGKQVTVNNIGTKEKAQ---FTASGN----
        Gigaspora_rosea                  --------------------------------------------------
        Rhizoctonia_solani               -------------------------N------------------------
        Stemphylium_lycopersici          --------------------------------------------------
        Globodera_rostochiensis_EXPB1    --------------------------------------------------
        Bursaphelenchus_xylophilus_EXP1  --------------------S----------------------CDKTQVT
        Meloidogyne_javanica_AP1         --------------------RCEKLG------------------------
        Faustovirus                      --------------------------------------------------
                                                                                        

        Xanthomonas_campestris           --
        Bacillus_subtilis                E-
        Xanthomonas_oryzae               --
        Xylella_fastidiosa               E-
        Streptomyces_mirabilis           QH
        Abortiporus_biennis              --
        Colletotrichum_destructivum      -C
        Gigaspora_rosea                  --
        Rhizoctonia_solani               --
        Stemphylium_lycopersici          --
        Globodera_rostochiensis_EXPB1    --
        Bursaphelenchus_xylophilus_EXP1  KC
        Meloidogyne_javanica_AP1         --
        Faustovirus                      --
                                        





        OUTPUT RESULTS
            #### File Type= GUIDE_TREE      Format= newick          Name= AllSeqs.dnd
            #### File Type= MSA             Format= aln             Name= AllSeqs.aln
            #### File Type= MSA             Format= html            Name= AllSeqs.html

        # Command Line: t_coffee AllSeqs.fasta  [PROGRAM:T-COFFEE]
        # T-COFFEE Memory Usage: Current= 36.313 Mb, Max= 38.970 Mb
        # Results Produced with T-COFFEE Version_12.00.7fb08c2 (2018-12-11 09:27:12 - Revision 7fb08c2 - Build 211)
        # T-COFFEE is available from http://www.tcoffee.org
        # Register on: https://groups.google.com/group/tcoffee/
    
        IN R STUDIO:
            > library(ape)
            > library(adegenet)
            > library(phangorn)
            > library(seqinr)
            > 
            > MSA_TCoffee <- read.alignment(file="allseqs.aln", format="clustal")
            > class(MSA_TCoffee)
            [1] "alignment"
            > #compute the distances
            > D_TCoffee <- dist.alignment(MSA_TCoffee)
            > #make the tree
            > tree_TCoffee <- nj(D_TCoffee)
            > Ltree_TCoffee <- ladderize(tree_TCoffee)
            > #plot and edit the tree
            > plot(Ltre_TCoffee, cex=0.6)
            Error in h(simpleError(msg, call)) : 
            error in evaluating the argument 'x' in selecting a method for function 'plot': object 'Ltre_TCoffee' not found
            > #plot and edit the tree
            > plot(Ltree_TCoffee, cex=0.6)
            > ?root
            > root(Ltree_TCoffee, Faustovirus)
            Error: object 'Faustovirus' not found
            > root(Ltree_TCoffee, outgroup = Faustovirus)
            Error: object 'Faustovirus' not found
            > root(Ltree_TCoffee, outgroup=Faustovirus)
            Error: object 'Faustovirus' not found
            > root(Ltree_TCoffee, outgroup= "Faustovirus")

            Phylogenetic tree with 14 tips and 12 internal nodes.

            Tip labels:
            Xanthomonas_campestris, Bacillus_subtilis, Xanthomonas_oryzae, Xylella_fastidiosa, Streptomyces_mirabilis, Abortiporus_biennis, ...

            Unrooted; includes branch length(s).
            > plot(Ltree_TCoffee)
            > tree_TCoffee <- nj(D_TCoffee)
            > Ltree_TCoffee <- ladderize(tree_TCoffee)
            > #is it rooted?
            > is.rooted(Ltree_TCoffee)
            [1] FALSE
            > #no, let's root it with the outgroup
            > root(Ltree_TCoffee, outgroup= "Faustovirus")

            Phylogenetic tree with 14 tips and 12 internal nodes.

            Tip labels:
            Xanthomonas_campestris, Bacillus_subtilis, Xanthomonas_oryzae, Xylella_fastidiosa, Streptomyces_mirabilis, Abortiporus_biennis, ...

            Unrooted; includes branch length(s).
            > tree_TCoffee <- nj(D_TCoffee)
            > Ltree_TCoffee <- ladderize(tree_TCoffee)
            > #no, let's root it with the outgroup
            > RootLtree_TCoffee <- root(Ltree_TCoffee, outgroup= "Faustovirus")
            > is.rooted(RootLtree_TCoffee)
            [1] FALSE
            > #that didn't work, need to add another argument
            > RootLtree_TCoffee <- root(Ltree_TCoffee, outgroup = "Faustovirus", resolve.root = TRUE)
            > is.rooted(RootLtree_TCoffee)
            [1] TRUE
            > plot(RootLtree_TCoffee)
            > plot(RootLtree_TCoffee, cex = 0.6)
            > plot(RootLtree_TCoffee, cex = 0.8)

            Tree saved: Desktop/PP563/Project/TCoffee_RootedDistanceTree.jpeg

    MAFFT Alignment:

        (base) Melettes-MBP-6:project melettedevore$ mafft-linsi AllSeqs.fasta > mafft-aligned-seqs.fast
        outputhat23=16
        treein = 0
        compacttree = 0
        stacksize: 8192 kb
        rescale = 1
        All-to-all alignment.
        tbfast-pair (aa) Version 7.525
        alg=L, model=BLOSUM62, 2.00, -0.10, +0.10, noshift, amax=0.0
        0 thread(s)

        outputhat23=16
        Loading 'hat3.seed' ... 
        done.
        Writing hat3 for iterative refinement
        rescale = 1
        Gap Penalty = -1.53, +0.00, +0.00
        tbutree = 1, compacttree = 0
        Constructing a UPGMA tree ... 
        10 / 14
        done.

        Progressive alignment ... 
        STEP    13 /13 
        done.
        tbfast (aa) Version 7.525
        alg=A, model=BLOSUM62, 1.53, -0.00, -0.00, noshift, amax=0.0
        1 thread(s)

        minimumweight = 0.000010
        autosubalignment = 0.000000
        nthread = 0
        randomseed = 0
        blosum 62 / kimura 200
        poffset = 0
        niter = 16
        sueff_global = 0.100000
        nadd = 16
        Loading 'hat3' ... done.
        rescale = 1

        10 / 14
        Segment   1/  1    1- 798
        STEP 008-013-1  rejected..   
        Converged.

        done
        dvtditr (aa) Version 7.525
        alg=A, model=BLOSUM62, 1.53, -0.00, -0.00, noshift, amax=0.0
        0 thread(s)


        Strategy:
        L-INS-i (Probably most accurate, very slow)
        Iterative refinement method (<16) with LOCAL pairwise alignment information

        If unsure which option to use, try 'mafft --auto input > output'.
        For more information, see 'mafft --help', 'mafft --man' and the mafft page.

        The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
        It tends to insert more gaps into gap-rich regions than previous versions.
        To disable this change, add the --leavegappyregion option.

        IN R STUDIO:
            > MSA_MAFFT <- read.alignment(file="mafft-aligned-seqs.fasta", format="fasta")
            > D_MAFFT <- dist.alignment(MSA_MAFFT)
            > tree_MAFFT <- nj(D_MAFFT)
            > Ltree_MAFFT <- ladderize(tree_MAFFT)
            > plot(Ltree_MAFFT, cex=0.6)
            > RootLtree_MAFFT <- root(Ltree_MAFFT, outgroup = "Faustovirus", resolve.root = TRUE)
            > plot(RootLtree_MAFFT, cex=0.8)

            Tree saved: Desktop/PP563/Project/MAFFT_RootedDistanceTree.jpeg
    
   INFERENCE METHOD #2: MAXIMUM LIKELIHOOD

   Converted T-Coffee MSA to fasta: 
        (base) Melettes-MBP-6:project melettedevore$ t_coffee -other_pg seq_reformat -in AllSeqs.aln -output phylip_aln > TCoffeeAllSeqsMSA.fasta
    
    Ran RaxML on T-Coffee Alignment:
        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --check --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix pT1

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 03-May-2025 17:52:45 as follows:

        raxml-ng --check --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix pT1

        Analysis options:
        run mode: Alignment validation
        start tree(s): 
        random seed: 1746312765
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: TCoffeeAllSeqsMSA.fasta
        [00:00:00] Loaded alignment with 14 taxa and 852 sites

        ERROR: Sequences 1 and 3 have identical name: Xanthomona

        ERROR: Duplicate sequence names found: 1

        ERROR: Alignment check failed (see details above)!

            **To fix this error I renamed one of the Xanthomonads to include part of the species name, this change was made in the TCoffee MSA Fasta file

        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --check --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix pT1

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 03-May-2025 17:53:44 as follows:

        raxml-ng --check --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix pT1

        Analysis options:
        run mode: Alignment validation
        start tree(s): 
        random seed: 1746312824
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: TCoffeeAllSeqsMSA.fasta
        [00:00:00] Loaded alignment with 14 taxa and 852 sites

        Alignment comprises 1 partitions and 852 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 852 / 852
        Gaps: 65.48 %
        Invariant sites: 37.21 %


        Alignment can be successfully read by RAxML-NG.


        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/pT1.raxml.log

        Analysis started: 03-May-2025 17:53:44 / finished: 03-May-2025 17:53:44

        Elapsed time: 0.007 seconds

        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --parse --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix pT2

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 03-May-2025 17:57:46 as follows:

        raxml-ng --parse --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix pT2

        Analysis options:
        run mode: Alignment parsing and compression
        start tree(s): 
        random seed: 1746313066
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: TCoffeeAllSeqsMSA.fasta
        [00:00:00] Loaded alignment with 14 taxa and 852 sites

        Alignment comprises 1 partitions and 660 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 852 / 660
        Gaps: 65.48 %
        Invariant sites: 37.21 %


        NOTE: Binary MSA file created: pT2.raxml.rba

        * Estimated memory requirements                : 11 MB

        * Recommended number of threads / MPI processes: 8

        Please note that numbers given above are rough estimates only. 
        Actual memory consumption and parallel performance on your system may differ!

        Alignment can be successfully read by RAxML-NG.


        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/pT2.raxml.log

        Analysis started: 03-May-2025 17:57:46 / finished: 03-May-2025 17:57:46

        Elapsed time: 0.005 seconds

        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix T3 --threads 2 --seed 2

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 03-May-2025 17:59:14 as follows:

        raxml-ng --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix T3 --threads 2 --seed 2

        Analysis options:
        run mode: ML tree search
        start tree(s): random (10) + parsimony (10)
        random seed: 2
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        fast spr radius: AUTO
        spr subtree cutoff: 1.000000
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: TCoffeeAllSeqsMSA.fasta
        [00:00:00] Loaded alignment with 14 taxa and 852 sites

        Alignment comprises 1 partitions and 660 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 852 / 660
        Gaps: 65.48 %
        Invariant sites: 37.21 %


        NOTE: Binary MSA file created: T3.raxml.rba

        [00:00:00] Generating 10 random starting tree(s) with 14 taxa
        [00:00:00] Generating 10 parsimony starting tree(s) with 14 taxa
        [00:00:00] Data distribution: max. partitions/sites/weight per thread: 1 / 330 / 26400

        Starting ML tree search with 20 distinct starting trees

        [00:00:00 -10948.385299] Initial branch length optimization
        [00:00:00 -9661.525996] Model parameter optimization (eps = 10.000000)
        [00:00:00 -9659.254433] AUTODETECT spr round 1 (radius: 5)
        [00:00:00 -9491.131985] AUTODETECT spr round 2 (radius: 10)
        [00:00:00 -9491.123641] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:00 -9491.123641] Model parameter optimization (eps = 3.000000)
        [00:00:00 -9484.734460] FAST spr round 1 (radius: 5)
        [00:00:01 -9473.028870] FAST spr round 2 (radius: 5)
        [00:00:02 -9473.028638] Model parameter optimization (eps = 1.000000)
        [00:00:02 -9472.201733] SLOW spr round 1 (radius: 5)
        [00:00:04 -9472.181961] SLOW spr round 2 (radius: 10)
        [00:00:07 -9472.181641] Model parameter optimization (eps = 0.100000)

        [00:00:07] ML tree search #1, logLikelihood: -9471.969852

        [00:00:07 -11042.853082] Initial branch length optimization
        [00:00:08 -9659.447401] Model parameter optimization (eps = 10.000000)
        [00:00:08 -9654.428426] AUTODETECT spr round 1 (radius: 5)
        [00:00:09 -9491.714853] AUTODETECT spr round 2 (radius: 10)
        [00:00:09 -9491.700933] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:09 -9491.700933] Model parameter optimization (eps = 3.000000)
        [00:00:09 -9487.430237] FAST spr round 1 (radius: 5)
        [00:00:11 -9473.183597] FAST spr round 2 (radius: 5)
        [00:00:11 -9473.182951] Model parameter optimization (eps = 1.000000)
        [00:00:11 -9472.257990] SLOW spr round 1 (radius: 5)
        [00:00:14 -9472.216438] SLOW spr round 2 (radius: 10)
        [00:00:15 -9472.216035] Model parameter optimization (eps = 0.100000)

        [00:00:15] ML tree search #2, logLikelihood: -9471.973265

        [00:00:15 -10919.469277] Initial branch length optimization
        [00:00:15 -9624.886832] Model parameter optimization (eps = 10.000000)
        [00:00:15 -9620.861663] AUTODETECT spr round 1 (radius: 5)
        [00:00:15 -9486.418542] AUTODETECT spr round 2 (radius: 10)
        [00:00:15 -9486.401311] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:15 -9486.401311] Model parameter optimization (eps = 3.000000)
        [00:00:15 -9481.237211] FAST spr round 1 (radius: 5)
        [00:00:16 -9472.950049] FAST spr round 2 (radius: 5)
        [00:00:17 -9472.948429] Model parameter optimization (eps = 1.000000)
        [00:00:17 -9472.190566] SLOW spr round 1 (radius: 5)
        [00:00:19 -9472.163430] SLOW spr round 2 (radius: 10)
        [00:00:20 -9472.163095] Model parameter optimization (eps = 0.100000)

        [00:00:20] ML tree search #3, logLikelihood: -9471.969471

        [00:00:20 -11110.547276] Initial branch length optimization
        [00:00:20 -9682.483338] Model parameter optimization (eps = 10.000000)
        [00:00:20 -9679.863656] AUTODETECT spr round 1 (radius: 5)
        [00:00:20 -9497.352543] AUTODETECT spr round 2 (radius: 10)
        [00:00:20 -9497.348633] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:20 -9497.348633] Model parameter optimization (eps = 3.000000)
        [00:00:20 -9488.766542] FAST spr round 1 (radius: 5)
        [00:00:21 -9472.696147] FAST spr round 2 (radius: 5)
        [00:00:22 -9472.695889] Model parameter optimization (eps = 1.000000)
        [00:00:22 -9472.131102] SLOW spr round 1 (radius: 5)
        [00:00:24 -9472.106961] SLOW spr round 2 (radius: 10)
        [00:00:24 -9472.106602] Model parameter optimization (eps = 0.100000)

        [00:00:25] ML tree search #4, logLikelihood: -9471.965267

        [00:00:25 -10878.558855] Initial branch length optimization
        [00:00:25 -9622.461964] Model parameter optimization (eps = 10.000000)
        [00:00:25 -9618.223290] AUTODETECT spr round 1 (radius: 5)
        [00:00:25 -9483.918288] AUTODETECT spr round 2 (radius: 10)
        [00:00:25 -9483.910442] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:25 -9483.910442] Model parameter optimization (eps = 3.000000)
        [00:00:25 -9478.588506] FAST spr round 1 (radius: 5)
        [00:00:26 -9472.827925] FAST spr round 2 (radius: 5)
        [00:00:27 -9472.402965] FAST spr round 3 (radius: 5)
        [00:00:27 -9472.402454] Model parameter optimization (eps = 1.000000)
        [00:00:28 -9472.058087] SLOW spr round 1 (radius: 5)
        [00:00:29 -9472.047333] SLOW spr round 2 (radius: 10)
        [00:00:30 -9472.047248] Model parameter optimization (eps = 0.100000)

        [00:00:30] ML tree search #5, logLikelihood: -9471.976197

        [00:00:30 -11187.526405] Initial branch length optimization
        [00:00:30 -9688.649568] Model parameter optimization (eps = 10.000000)
        [00:00:30 -9687.393458] AUTODETECT spr round 1 (radius: 5)
        [00:00:31 -9495.402509] AUTODETECT spr round 2 (radius: 10)
        [00:00:31 -9495.400811] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:31 -9495.400811] Model parameter optimization (eps = 3.000000)
        [00:00:31 -9486.659343] FAST spr round 1 (radius: 5)
        [00:00:32 -9475.024570] FAST spr round 2 (radius: 5)
        [00:00:33 -9472.831198] FAST spr round 3 (radius: 5)
        [00:00:33 -9472.826644] Model parameter optimization (eps = 1.000000)
        [00:00:34 -9472.154483] SLOW spr round 1 (radius: 5)
        [00:00:35 -9472.134760] SLOW spr round 2 (radius: 10)
        [00:00:36 -9472.134482] Model parameter optimization (eps = 0.100000)

        [00:00:36] ML tree search #6, logLikelihood: -9471.966371

        [00:00:36 -11129.406873] Initial branch length optimization
        [00:00:36 -9669.882281] Model parameter optimization (eps = 10.000000)
        [00:00:37 -9667.653349] AUTODETECT spr round 1 (radius: 5)
        [00:00:37 -9490.115637] AUTODETECT spr round 2 (radius: 10)
        [00:00:37 -9490.106768] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:37 -9490.106768] Model parameter optimization (eps = 3.000000)
        [00:00:37 -9481.717493] FAST spr round 1 (radius: 5)
        [00:00:38 -9472.671537] FAST spr round 2 (radius: 5)
        [00:00:39 -9472.671130] Model parameter optimization (eps = 1.000000)
        [00:00:39 -9472.120882] SLOW spr round 1 (radius: 5)
        [00:00:40 -9472.101358] SLOW spr round 2 (radius: 10)
        [00:00:42 -9472.101222] Model parameter optimization (eps = 0.100000)

        [00:00:42] ML tree search #7, logLikelihood: -9471.963973

        [00:00:42 -10921.746373] Initial branch length optimization
        [00:00:42 -9612.450670] Model parameter optimization (eps = 10.000000)
        [00:00:42 -9606.889602] AUTODETECT spr round 1 (radius: 5)
        [00:00:42 -9483.350616] AUTODETECT spr round 2 (radius: 10)
        [00:00:42 -9483.318365] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:42 -9483.318365] Model parameter optimization (eps = 3.000000)
        [00:00:42 -9478.467320] FAST spr round 1 (radius: 5)
        [00:00:43 -9472.707882] FAST spr round 2 (radius: 5)
        [00:00:44 -9472.707559] Model parameter optimization (eps = 1.000000)
        [00:00:44 -9472.125882] SLOW spr round 1 (radius: 5)
        [00:00:46 -9472.109552] SLOW spr round 2 (radius: 10)
        [00:00:47 -9472.109414] Model parameter optimization (eps = 0.100000)

        [00:00:47] ML tree search #8, logLikelihood: -9471.964939

        [00:00:47 -11108.092009] Initial branch length optimization
        [00:00:47 -9678.744591] Model parameter optimization (eps = 10.000000)
        [00:00:47 -9676.296901] AUTODETECT spr round 1 (radius: 5)
        [00:00:47 -9503.814703] AUTODETECT spr round 2 (radius: 10)
        [00:00:47 -9503.811891] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:47 -9503.811891] Model parameter optimization (eps = 3.000000)
        [00:00:47 -9496.814367] FAST spr round 1 (radius: 5)
        [00:00:48 -9472.520762] FAST spr round 2 (radius: 5)
        [00:00:49 -9472.519713] Model parameter optimization (eps = 1.000000)
        [00:00:49 -9472.085471] SLOW spr round 1 (radius: 5)
        [00:00:51 -9472.072360] SLOW spr round 2 (radius: 10)
        [00:00:52 -9472.072229] Model parameter optimization (eps = 0.100000)

        [00:00:53] ML tree search #9, logLikelihood: -9471.981657

        [00:00:53 -10885.503662] Initial branch length optimization
        [00:00:53 -9660.124087] Model parameter optimization (eps = 10.000000)
        [00:00:53 -9658.505984] AUTODETECT spr round 1 (radius: 5)
        [00:00:53 -9484.271692] AUTODETECT spr round 2 (radius: 10)
        [00:00:53 -9484.266252] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:53 -9484.266252] Model parameter optimization (eps = 3.000000)
        [00:00:53 -9476.727103] FAST spr round 1 (radius: 5)
        [00:00:54 -9472.996262] FAST spr round 2 (radius: 5)
        [00:00:55 -9472.995024] Model parameter optimization (eps = 1.000000)
        [00:00:55 -9472.197861] SLOW spr round 1 (radius: 5)
        [00:00:57 -9472.175463] SLOW spr round 2 (radius: 10)
        [00:00:57 -9472.175304] Model parameter optimization (eps = 0.100000)

        [00:00:58] ML tree search #10, logLikelihood: -9471.969632

        [00:00:58 -10330.296723] Initial branch length optimization
        [00:00:58 -9501.178552] Model parameter optimization (eps = 10.000000)
        [00:00:58 -9493.270678] AUTODETECT spr round 1 (radius: 5)
        [00:00:58 -9476.812950] AUTODETECT spr round 2 (radius: 10)
        [00:00:58 -9476.808001] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:58 -9476.808001] Model parameter optimization (eps = 3.000000)
        [00:00:58 -9472.412711] FAST spr round 1 (radius: 5)
        [00:00:59 -9472.347790] Model parameter optimization (eps = 1.000000)
        [00:00:59 -9472.044436] SLOW spr round 1 (radius: 5)
        [00:01:01 -9472.032668] SLOW spr round 2 (radius: 10)
        [00:01:02 -9472.032525] Model parameter optimization (eps = 0.100000)

        [00:01:02] ML tree search #11, logLikelihood: -9471.972742

        [00:01:02 -10407.303204] Initial branch length optimization
        [00:01:02 -9503.863139] Model parameter optimization (eps = 10.000000)
        [00:01:02 -9497.201010] AUTODETECT spr round 1 (radius: 5)
        [00:01:02 -9477.593692] AUTODETECT spr round 2 (radius: 10)
        [00:01:02 -9477.584903] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:02 -9477.584903] Model parameter optimization (eps = 3.000000)
        [00:01:02 -9473.115096] FAST spr round 1 (radius: 5)
        [00:01:03 -9472.335053] FAST spr round 2 (radius: 5)
        [00:01:04 -9472.332852] Model parameter optimization (eps = 1.000000)
        [00:01:04 -9472.041756] SLOW spr round 1 (radius: 5)
        [00:01:06 -9472.031021] SLOW spr round 2 (radius: 10)
        [00:01:07 -9472.030930] Model parameter optimization (eps = 0.100000)

        [00:01:07] ML tree search #12, logLikelihood: -9471.972287

        [00:01:07 -10383.606514] Initial branch length optimization
        [00:01:07 -9504.138507] Model parameter optimization (eps = 10.000000)
        [00:01:07 -9495.393342] AUTODETECT spr round 1 (radius: 5)
        [00:01:07 -9476.459987] AUTODETECT spr round 2 (radius: 10)
        [00:01:07 -9476.451995] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:07 -9476.451995] Model parameter optimization (eps = 3.000000)
        [00:01:07 -9472.721739] FAST spr round 1 (radius: 5)
        [00:01:08 -9472.238586] FAST spr round 2 (radius: 5)
        [00:01:09 -9472.237861] Model parameter optimization (eps = 1.000000)
        [00:01:09 -9472.020460] SLOW spr round 1 (radius: 5)
        [00:01:11 -9472.011766] SLOW spr round 2 (radius: 10)
        [00:01:12 -9472.011677] Model parameter optimization (eps = 0.100000)

        [00:01:12] ML tree search #13, logLikelihood: -9471.968091

        [00:01:12 -10346.023462] Initial branch length optimization
        [00:01:12 -9490.521686] Model parameter optimization (eps = 10.000000)
        [00:01:12 -9482.023697] AUTODETECT spr round 1 (radius: 5)
        [00:01:12 -9476.572317] AUTODETECT spr round 2 (radius: 10)
        [00:01:12 -9476.569638] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:12 -9476.569638] Model parameter optimization (eps = 3.000000)
        [00:01:12 -9472.369046] FAST spr round 1 (radius: 5)
        [00:01:13 -9472.305639] Model parameter optimization (eps = 1.000000)
        [00:01:13 -9472.035636] SLOW spr round 1 (radius: 5)
        [00:01:15 -9472.025385] SLOW spr round 2 (radius: 10)
        [00:01:16 -9472.025272] Model parameter optimization (eps = 0.100000)

        [00:01:16] ML tree search #14, logLikelihood: -9471.970960

        [00:01:16 -10409.505415] Initial branch length optimization
        [00:01:16 -9511.963784] Model parameter optimization (eps = 10.000000)
        [00:01:16 -9504.537838] AUTODETECT spr round 1 (radius: 5)
        [00:01:16 -9478.775365] AUTODETECT spr round 2 (radius: 10)
        [00:01:16 -9478.774706] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:16 -9478.774706] Model parameter optimization (eps = 3.000000)
        [00:01:16 -9474.951270] FAST spr round 1 (radius: 5)
        [00:01:17 -9472.395443] FAST spr round 2 (radius: 5)
        [00:01:18 -9472.395392] Model parameter optimization (eps = 1.000000)
        [00:01:18 -9472.057028] SLOW spr round 1 (radius: 5)
        [00:01:20 -9472.043234] SLOW spr round 2 (radius: 10)
        [00:01:21 -9472.043102] Model parameter optimization (eps = 0.100000)

        [00:01:21] ML tree search #15, logLikelihood: -9471.975522

        [00:01:21 -10340.411000] Initial branch length optimization
        [00:01:21 -9496.336912] Model parameter optimization (eps = 10.000000)
        [00:01:21 -9488.657343] AUTODETECT spr round 1 (radius: 5)
        [00:01:21 -9476.668059] AUTODETECT spr round 2 (radius: 10)
        [00:01:21 -9476.665751] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:21 -9476.665751] Model parameter optimization (eps = 3.000000)
        [00:01:21 -9472.363170] FAST spr round 1 (radius: 5)
        [00:01:22 -9472.303065] Model parameter optimization (eps = 1.000000)
        [00:01:22 -9472.035103] SLOW spr round 1 (radius: 5)
        [00:01:24 -9472.024528] SLOW spr round 2 (radius: 10)
        [00:01:25 -9472.024445] Model parameter optimization (eps = 0.100000)

        [00:01:25] ML tree search #16, logLikelihood: -9471.970749

        [00:01:25 -10421.720069] Initial branch length optimization
        [00:01:25 -9525.026206] Model parameter optimization (eps = 10.000000)
        [00:01:25 -9518.306134] AUTODETECT spr round 1 (radius: 5)
        [00:01:25 -9479.166088] AUTODETECT spr round 2 (radius: 10)
        [00:01:25 -9479.163990] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:25 -9479.163990] Model parameter optimization (eps = 3.000000)
        [00:01:25 -9473.655421] FAST spr round 1 (radius: 5)
        [00:01:26 -9472.846821] FAST spr round 2 (radius: 5)
        [00:01:27 -9472.420033] FAST spr round 3 (radius: 5)
        [00:01:28 -9472.419377] Model parameter optimization (eps = 1.000000)
        [00:01:28 -9472.061850] SLOW spr round 1 (radius: 5)
        [00:01:30 -9472.048942] SLOW spr round 2 (radius: 10)
        [00:01:31 -9472.048812] Model parameter optimization (eps = 0.100000)

        [00:01:31] ML tree search #17, logLikelihood: -9471.976038

        [00:01:31 -10351.922093] Initial branch length optimization
        [00:01:31 -9493.558264] Model parameter optimization (eps = 10.000000)
        [00:01:31 -9485.877049] AUTODETECT spr round 1 (radius: 5)
        [00:01:31 -9477.518204] AUTODETECT spr round 2 (radius: 10)
        [00:01:31 -9477.513618] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:31 -9477.513618] Model parameter optimization (eps = 3.000000)
        [00:01:31 -9472.470516] FAST spr round 1 (radius: 5)
        [00:01:32 -9472.383537] Model parameter optimization (eps = 1.000000)
        [00:01:32 -9472.053546] SLOW spr round 1 (radius: 5)
        [00:01:34 -9472.039779] SLOW spr round 2 (radius: 10)
        [00:01:35 -9472.039547] Model parameter optimization (eps = 0.100000)

        [00:01:35] ML tree search #18, logLikelihood: -9471.974305

        [00:01:35 -10384.504194] Initial branch length optimization
        [00:01:35 -9502.056645] Model parameter optimization (eps = 10.000000)
        [00:01:35 -9492.806155] AUTODETECT spr round 1 (radius: 5)
        [00:01:35 -9477.125251] AUTODETECT spr round 2 (radius: 10)
        [00:01:36 -9477.118019] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:36 -9477.118019] Model parameter optimization (eps = 3.000000)
        [00:01:36 -9474.684657] FAST spr round 1 (radius: 5)
        [00:01:36 -9472.939761] FAST spr round 2 (radius: 5)
        [00:01:37 -9472.938826] Model parameter optimization (eps = 1.000000)
        [00:01:37 -9472.190916] SLOW spr round 1 (radius: 5)
        [00:01:39 -9472.160661] SLOW spr round 2 (radius: 10)
        [00:01:40 -9472.160200] Model parameter optimization (eps = 0.100000)

        [00:01:40] ML tree search #19, logLikelihood: -9471.968975

        [00:01:40 -10338.397028] Initial branch length optimization
        [00:01:40 -9493.744161] Model parameter optimization (eps = 10.000000)
        [00:01:40 -9484.655106] AUTODETECT spr round 1 (radius: 5)
        [00:01:40 -9476.347898] AUTODETECT spr round 2 (radius: 10)
        [00:01:41 -9476.338818] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:41 -9476.338818] Model parameter optimization (eps = 3.000000)
        [00:01:41 -9472.346948] FAST spr round 1 (radius: 5)
        [00:01:41 -9472.276953] Model parameter optimization (eps = 1.000000)
        [00:01:41 -9472.029998] SLOW spr round 1 (radius: 5)
        [00:01:43 -9472.019497] SLOW spr round 2 (radius: 10)
        [00:01:44 -9472.019406] Model parameter optimization (eps = 0.100000)

        [00:01:44] ML tree search #20, logLikelihood: -9471.969976


        Optimized model parameters:

        Partition 0: noname
        Rate heterogeneity: GAMMA (4 cats, mean),  alpha: 1.800587 (ML),  weights&rates: (0.250000,0.268250) (0.250000,0.631982) (0.250000,1.064093) (0.250000,2.035674) 
        Base frequencies (model): 0.076748 0.051691 0.042645 0.051544 0.019803 0.040752 0.061830 0.073152 0.022944 0.053761 0.091904 0.058676 0.023826 0.040126 0.050901 0.068765 0.058565 0.014261 0.032102 0.066005 
        Substitution rates (model): 58.000000 54.000000 81.000000 56.000000 57.000000 105.000000 179.000000 27.000000 36.000000 30.000000 35.000000 54.000000 15.000000 194.000000 378.000000 475.000000 9.000000 11.000000 298.000000 45.000000 16.000000 113.000000 310.000000 29.000000 137.000000 328.000000 22.000000 38.000000 646.000000 44.000000 5.000000 74.000000 101.000000 64.000000 126.000000 20.000000 17.000000 528.000000 34.000000 86.000000 58.000000 81.000000 391.000000 47.000000 12.000000 263.000000 30.000000 10.000000 15.000000 503.000000 232.000000 8.000000 70.000000 16.000000 10.000000 49.000000 767.000000 130.000000 112.000000 11.000000 7.000000 26.000000 15.000000 4.000000 15.000000 59.000000 38.000000 4.000000 46.000000 31.000000 9.000000 5.000000 59.000000 69.000000 17.000000 23.000000 7.000000 31.000000 78.000000 14.000000 223.000000 42.000000 115.000000 209.000000 62.000000 323.000000 26.000000 597.000000 9.000000 72.000000 292.000000 43.000000 4.000000 164.000000 53.000000 51.000000 18.000000 24.000000 20.000000 119.000000 26.000000 12.000000 9.000000 181.000000 18.000000 5.000000 18.000000 30.000000 32.000000 10.000000 7.000000 45.000000 23.000000 6.000000 6.000000 27.000000 14.000000 5.000000 24.000000 201.000000 33.000000 55.000000 8.000000 47.000000 16.000000 56.000000 45.000000 33.000000 40.000000 115.000000 73.000000 46.000000 8.000000 573.000000 11.000000 229.000000 21.000000 479.000000 89.000000 10.000000 40.000000 245.000000 9.000000 32.000000 961.000000 14.000000 388.000000 248.000000 102.000000 59.000000 25.000000 52.000000 24.000000 180.000000 65.000000 4.000000 21.000000 47.000000 103.000000 10.000000 8.000000 14.000000 43.000000 16.000000 29.000000 226.000000 24.000000 18.000000 323.000000 17.000000 92.000000 12.000000 53.000000 536.000000 62.000000 285.000000 118.000000 6.000000 10.000000 23.000000 477.000000 35.000000 63.000000 38.000000 12.000000 21.000000 112.000000 71.000000 25.000000 16.000000 


        Final LogLikelihood: -9471.963973

        AIC score: 18995.927945 / AICc score: 18997.629764 / BIC score: 19119.365195
        Free parameters (model + branch lengths): 26

        Best ML tree saved to: /Users/melettedevore/Desktop/PP563/project/T3.raxml.bestTree
        All ML trees saved to: /Users/melettedevore/Desktop/PP563/project/T3.raxml.mlTrees
        Optimized model saved to: /Users/melettedevore/Desktop/PP563/project/T3.raxml.bestModel

        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/T3.raxml.log

        Analysis started: 03-May-2025 17:59:14 / finished: 03-May-2025 18:00:59

        Elapsed time: 104.510 seconds

    * Re-ran with bootstrap analysis

        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --all --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix T4 --threads 2 --seed 2

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 04-May-2025 14:51:53 as follows:

        raxml-ng --all --msa TCoffeeAllSeqsMSA.fasta --model JTT+G --prefix T4 --threads 2 --seed 2

        Analysis options:
        run mode: ML tree search + bootstrapping (Felsenstein Bootstrap)
        start tree(s): random (10) + parsimony (10)
        bootstrap replicates: max: 1000 + bootstopping (autoMRE, cutoff: 0.030000)
        random seed: 2
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: TCoffeeAllSeqsMSA.fasta
        [00:00:00] Loaded alignment with 14 taxa and 852 sites

        Alignment comprises 1 partitions and 660 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 852 / 660
        Gaps: 65.48 %
        Invariant sites: 37.21 %


        NOTE: Binary MSA file created: T4.raxml.rba

        [00:00:00] Generating 10 random starting tree(s) with 14 taxa
        [00:00:00] Generating 10 parsimony starting tree(s) with 14 taxa
        [00:00:00] Data distribution: max. partitions/sites/weight per thread: 1 / 330 / 26400

        Starting ML tree search with 20 distinct starting trees

        [00:00:05] ML tree search #1, logLikelihood: -9471.969852
        [00:00:10] ML tree search #2, logLikelihood: -9471.973265
        [00:00:15] ML tree search #3, logLikelihood: -9471.969471
        [00:00:20] ML tree search #4, logLikelihood: -9471.965267
        [00:00:26] ML tree search #5, logLikelihood: -9471.976197
        [00:00:32] ML tree search #6, logLikelihood: -9471.966371
        [00:00:37] ML tree search #7, logLikelihood: -9471.963973
        [00:00:45] ML tree search #8, logLikelihood: -9471.964939
        [00:00:53] ML tree search #9, logLikelihood: -9471.981657
        [00:01:01] ML tree search #10, logLikelihood: -9471.969632
        [00:01:07] ML tree search #11, logLikelihood: -9471.972742
        [00:01:13] ML tree search #12, logLikelihood: -9471.972287
        [00:01:20] ML tree search #13, logLikelihood: -9471.968091
        [00:01:25] ML tree search #14, logLikelihood: -9471.970960
        [00:01:31] ML tree search #15, logLikelihood: -9471.975522
        [00:01:34] ML tree search #16, logLikelihood: -9471.970749
        [00:01:40] ML tree search #17, logLikelihood: -9471.976038
        [00:01:44] ML tree search #18, logLikelihood: -9471.974305
        [00:01:49] ML tree search #19, logLikelihood: -9471.968975
        [00:01:53] ML tree search #20, logLikelihood: -9471.969976

        [00:01:53] ML tree search completed, best tree logLH: -9471.963973

        [00:01:53] Starting bootstrapping analysis with 1000 replicates.

        [00:01:59] Bootstrap tree #1, logLikelihood: -9637.902843
        [00:02:03] Bootstrap tree #2, logLikelihood: -9451.870641
        [00:02:07] Bootstrap tree #3, logLikelihood: -9330.580191
        [00:02:12] Bootstrap tree #4, logLikelihood: -9495.831411
        [00:02:19] Bootstrap tree #5, logLikelihood: -9752.312335
        [00:02:27] Bootstrap tree #6, logLikelihood: -9504.084536
        [00:02:30] Bootstrap tree #7, logLikelihood: -9337.897401
        [00:02:35] Bootstrap tree #8, logLikelihood: -9723.908464
        [00:02:41] Bootstrap tree #9, logLikelihood: -9234.888685
        [00:02:45] Bootstrap tree #10, logLikelihood: -9200.180925
        [00:02:49] Bootstrap tree #11, logLikelihood: -9698.326112
        [00:02:54] Bootstrap tree #12, logLikelihood: -8811.267494
        [00:02:57] Bootstrap tree #13, logLikelihood: -9599.996809
        [00:03:03] Bootstrap tree #14, logLikelihood: -9497.856585
        [00:03:08] Bootstrap tree #15, logLikelihood: -9722.551311
        [00:03:12] Bootstrap tree #16, logLikelihood: -9384.852088
        [00:03:17] Bootstrap tree #17, logLikelihood: -9642.355053
        [00:03:20] Bootstrap tree #18, logLikelihood: -9324.424970
        [00:03:23] Bootstrap tree #19, logLikelihood: -9206.533796
        [00:03:27] Bootstrap tree #20, logLikelihood: -9467.764387
        [00:03:31] Bootstrap tree #21, logLikelihood: -9442.452816
        [00:03:35] Bootstrap tree #22, logLikelihood: -9775.211335
        [00:03:40] Bootstrap tree #23, logLikelihood: -9407.447045
        [00:03:45] Bootstrap tree #24, logLikelihood: -9740.186026
        [00:03:50] Bootstrap tree #25, logLikelihood: -9495.662014
        [00:03:53] Bootstrap tree #26, logLikelihood: -9119.082520
        [00:03:58] Bootstrap tree #27, logLikelihood: -9797.439332
        [00:04:02] Bootstrap tree #28, logLikelihood: -9249.998258
        [00:04:05] Bootstrap tree #29, logLikelihood: -9611.342083
        [00:04:09] Bootstrap tree #30, logLikelihood: -9770.996339
        [00:04:13] Bootstrap tree #31, logLikelihood: -9067.899038
        [00:04:17] Bootstrap tree #32, logLikelihood: -9343.216489
        [00:04:21] Bootstrap tree #33, logLikelihood: -9637.884063
        [00:04:25] Bootstrap tree #34, logLikelihood: -9397.850196
        [00:04:31] Bootstrap tree #35, logLikelihood: -9155.237596
        [00:04:34] Bootstrap tree #36, logLikelihood: -9143.814336
        [00:04:38] Bootstrap tree #37, logLikelihood: -9193.364473
        [00:04:43] Bootstrap tree #38, logLikelihood: -9319.444345
        [00:04:46] Bootstrap tree #39, logLikelihood: -9287.842731
        [00:04:51] Bootstrap tree #40, logLikelihood: -9960.978365
        [00:04:55] Bootstrap tree #41, logLikelihood: -9278.905110
        [00:04:59] Bootstrap tree #42, logLikelihood: -9688.338584
        [00:05:03] Bootstrap tree #43, logLikelihood: -9503.007093
        [00:05:08] Bootstrap tree #44, logLikelihood: -9492.673393
        [00:05:11] Bootstrap tree #45, logLikelihood: -9617.076434
        [00:05:14] Bootstrap tree #46, logLikelihood: -9809.475465
        [00:05:18] Bootstrap tree #47, logLikelihood: -9475.665589
        [00:05:22] Bootstrap tree #48, logLikelihood: -9758.515552
        [00:05:25] Bootstrap tree #49, logLikelihood: -9659.317084
        [00:05:30] Bootstrap tree #50, logLikelihood: -9653.409907
        [00:05:33] Bootstrap tree #51, logLikelihood: -9358.940544
        [00:05:38] Bootstrap tree #52, logLikelihood: -9392.641045
        [00:05:44] Bootstrap tree #53, logLikelihood: -9391.762978
        [00:05:47] Bootstrap tree #54, logLikelihood: -8929.064694
        [00:05:50] Bootstrap tree #55, logLikelihood: -9328.429129
        [00:05:55] Bootstrap tree #56, logLikelihood: -9890.812196
        [00:06:00] Bootstrap tree #57, logLikelihood: -9499.974071
        [00:06:03] Bootstrap tree #58, logLikelihood: -9554.753117
        [00:06:07] Bootstrap tree #59, logLikelihood: -9745.524762
        [00:06:11] Bootstrap tree #60, logLikelihood: -9219.441649
        [00:06:16] Bootstrap tree #61, logLikelihood: -9751.118974
        [00:06:20] Bootstrap tree #62, logLikelihood: -9609.215424
        [00:06:24] Bootstrap tree #63, logLikelihood: -9196.046762
        [00:06:29] Bootstrap tree #64, logLikelihood: -9951.289253
        [00:06:35] Bootstrap tree #65, logLikelihood: -9539.935777
        [00:06:38] Bootstrap tree #66, logLikelihood: -9659.598886
        [00:06:43] Bootstrap tree #67, logLikelihood: -8886.732770
        [00:06:48] Bootstrap tree #68, logLikelihood: -9450.193810
        [00:06:52] Bootstrap tree #69, logLikelihood: -9722.072981
        [00:06:57] Bootstrap tree #70, logLikelihood: -9405.904860
        [00:07:01] Bootstrap tree #71, logLikelihood: -9311.956622
        [00:07:05] Bootstrap tree #72, logLikelihood: -8564.437463
        [00:07:08] Bootstrap tree #73, logLikelihood: -9376.768991
        [00:07:11] Bootstrap tree #74, logLikelihood: -9449.730982
        [00:07:16] Bootstrap tree #75, logLikelihood: -9449.012280
        [00:07:20] Bootstrap tree #76, logLikelihood: -9246.587511
        [00:07:24] Bootstrap tree #77, logLikelihood: -9139.714797
        [00:07:28] Bootstrap tree #78, logLikelihood: -9282.074859
        [00:07:32] Bootstrap tree #79, logLikelihood: -9783.092324
        [00:07:39] Bootstrap tree #80, logLikelihood: -9199.392574
        [00:07:44] Bootstrap tree #81, logLikelihood: -9740.395759
        [00:07:48] Bootstrap tree #82, logLikelihood: -9270.086754
        [00:07:52] Bootstrap tree #83, logLikelihood: -9321.895435
        [00:07:56] Bootstrap tree #84, logLikelihood: -9182.881692
        [00:08:00] Bootstrap tree #85, logLikelihood: -9517.418890
        [00:08:05] Bootstrap tree #86, logLikelihood: -9744.152884
        [00:08:09] Bootstrap tree #87, logLikelihood: -9536.229785
        [00:08:13] Bootstrap tree #88, logLikelihood: -9352.772009
        [00:08:17] Bootstrap tree #89, logLikelihood: -9613.673455
        [00:08:20] Bootstrap tree #90, logLikelihood: -9447.705545
        [00:08:24] Bootstrap tree #91, logLikelihood: -9893.755398
        [00:08:28] Bootstrap tree #92, logLikelihood: -9508.562219
        [00:08:31] Bootstrap tree #93, logLikelihood: -9452.701914
        [00:08:35] Bootstrap tree #94, logLikelihood: -9755.748703
        [00:08:39] Bootstrap tree #95, logLikelihood: -9535.274496
        [00:08:42] Bootstrap tree #96, logLikelihood: -9883.377041
        [00:08:46] Bootstrap tree #97, logLikelihood: -9129.195197
        [00:08:50] Bootstrap tree #98, logLikelihood: -9401.177649
        [00:08:55] Bootstrap tree #99, logLikelihood: -9501.684181
        [00:08:58] Bootstrap tree #100, logLikelihood: -9350.402448
        [00:09:06] Bootstrap tree #101, logLikelihood: -9451.383264
        [00:09:10] Bootstrap tree #102, logLikelihood: -9567.027553
        [00:09:14] Bootstrap tree #103, logLikelihood: -9176.197410
        [00:09:18] Bootstrap tree #104, logLikelihood: -9324.750830
        [00:09:23] Bootstrap tree #105, logLikelihood: -9192.476567
        [00:09:26] Bootstrap tree #106, logLikelihood: -9573.363681
        [00:09:30] Bootstrap tree #107, logLikelihood: -9657.166313
        [00:09:34] Bootstrap tree #108, logLikelihood: -9734.441070
        [00:09:38] Bootstrap tree #109, logLikelihood: -9128.118083
        [00:09:42] Bootstrap tree #110, logLikelihood: -9629.569254
        [00:09:46] Bootstrap tree #111, logLikelihood: -9759.098318
        [00:09:50] Bootstrap tree #112, logLikelihood: -9695.772997
        [00:09:54] Bootstrap tree #113, logLikelihood: -9178.681990
        [00:10:00] Bootstrap tree #114, logLikelihood: -9732.778824
        [00:10:03] Bootstrap tree #115, logLikelihood: -9271.114905
        [00:10:07] Bootstrap tree #116, logLikelihood: -9362.564309
        [00:10:14] Bootstrap tree #117, logLikelihood: -9678.880276
        [00:10:17] Bootstrap tree #118, logLikelihood: -9539.367044
        [00:10:20] Bootstrap tree #119, logLikelihood: -9433.355170
        [00:10:26] Bootstrap tree #120, logLikelihood: -9071.068597
        [00:10:29] Bootstrap tree #121, logLikelihood: -9553.697123
        [00:10:33] Bootstrap tree #122, logLikelihood: -9495.526914
        [00:10:37] Bootstrap tree #123, logLikelihood: -9370.979003
        [00:10:41] Bootstrap tree #124, logLikelihood: -9553.941695
        [00:10:46] Bootstrap tree #125, logLikelihood: -9519.697479
        [00:10:51] Bootstrap tree #126, logLikelihood: -9014.661608
        [00:10:54] Bootstrap tree #127, logLikelihood: -9699.416940
        [00:10:58] Bootstrap tree #128, logLikelihood: -9334.345032
        [00:11:02] Bootstrap tree #129, logLikelihood: -9532.498510
        [00:11:08] Bootstrap tree #130, logLikelihood: -9374.376060
        [00:11:13] Bootstrap tree #131, logLikelihood: -9627.479973
        [00:11:18] Bootstrap tree #132, logLikelihood: -9296.234848
        [00:11:22] Bootstrap tree #133, logLikelihood: -9607.996318
        [00:11:25] Bootstrap tree #134, logLikelihood: -9517.860298
        [00:11:28] Bootstrap tree #135, logLikelihood: -8652.440868
        [00:11:33] Bootstrap tree #136, logLikelihood: -9074.308717
        [00:11:37] Bootstrap tree #137, logLikelihood: -9569.527586
        [00:11:40] Bootstrap tree #138, logLikelihood: -9279.722078
        [00:11:44] Bootstrap tree #139, logLikelihood: -9752.764157
        [00:11:48] Bootstrap tree #140, logLikelihood: -9301.515104
        [00:11:52] Bootstrap tree #141, logLikelihood: -9134.880314
        [00:11:55] Bootstrap tree #142, logLikelihood: -9312.331451
        [00:12:02] Bootstrap tree #143, logLikelihood: -9279.991879
        [00:12:06] Bootstrap tree #144, logLikelihood: -9704.200350
        [00:12:10] Bootstrap tree #145, logLikelihood: -10082.147237
        [00:12:13] Bootstrap tree #146, logLikelihood: -9961.811295
        [00:12:17] Bootstrap tree #147, logLikelihood: -9265.516564
        [00:12:21] Bootstrap tree #148, logLikelihood: -9524.493264
        [00:12:26] Bootstrap tree #149, logLikelihood: -9631.574543
        [00:12:29] Bootstrap tree #150, logLikelihood: -9542.273480
        [00:12:33] Bootstrap tree #151, logLikelihood: -9389.285699
        [00:12:37] Bootstrap tree #152, logLikelihood: -9225.859524
        [00:12:42] Bootstrap tree #153, logLikelihood: -9485.897296
        [00:12:48] Bootstrap tree #154, logLikelihood: -9428.547152
        [00:12:51] Bootstrap tree #155, logLikelihood: -10103.364550
        [00:12:55] Bootstrap tree #156, logLikelihood: -9260.984287
        [00:13:02] Bootstrap tree #157, logLikelihood: -9935.928849
        [00:13:23] Bootstrap tree #158, logLikelihood: -9549.014387
        [00:13:33] Bootstrap tree #159, logLikelihood: -9441.748837
        [00:13:39] Bootstrap tree #160, logLikelihood: -9333.377361
        [00:13:44] Bootstrap tree #161, logLikelihood: -8959.477328
        [00:13:51] Bootstrap tree #162, logLikelihood: -9304.945079
        [00:13:57] Bootstrap tree #163, logLikelihood: -9239.118001
        [00:14:02] Bootstrap tree #164, logLikelihood: -9703.566736
        [00:14:08] Bootstrap tree #165, logLikelihood: -9138.667032
        [00:14:12] Bootstrap tree #166, logLikelihood: -9741.768965
        [00:14:16] Bootstrap tree #167, logLikelihood: -10211.589478
        [00:14:23] Bootstrap tree #168, logLikelihood: -9647.125303
        [00:14:29] Bootstrap tree #169, logLikelihood: -9177.973585
        [00:14:34] Bootstrap tree #170, logLikelihood: -9535.320073
        [00:14:40] Bootstrap tree #171, logLikelihood: -8951.832750
        [00:14:45] Bootstrap tree #172, logLikelihood: -9481.177251
        [00:14:50] Bootstrap tree #173, logLikelihood: -9004.424594
        [00:14:56] Bootstrap tree #174, logLikelihood: -9196.290802
        [00:15:01] Bootstrap tree #175, logLikelihood: -9317.092210
        [00:15:06] Bootstrap tree #176, logLikelihood: -9386.967085
        [00:15:09] Bootstrap tree #177, logLikelihood: -9854.635667
        [00:15:14] Bootstrap tree #178, logLikelihood: -9298.903406
        [00:15:19] Bootstrap tree #179, logLikelihood: -9456.460554
        [00:15:23] Bootstrap tree #180, logLikelihood: -9141.610739
        [00:15:27] Bootstrap tree #181, logLikelihood: -9626.015927
        [00:15:32] Bootstrap tree #182, logLikelihood: -9364.114204
        [00:15:36] Bootstrap tree #183, logLikelihood: -9611.120504
        [00:15:41] Bootstrap tree #184, logLikelihood: -9243.183985
        [00:15:47] Bootstrap tree #185, logLikelihood: -9740.625383
        [00:15:52] Bootstrap tree #186, logLikelihood: -9041.989820
        [00:15:56] Bootstrap tree #187, logLikelihood: -8953.198544
        [00:16:00] Bootstrap tree #188, logLikelihood: -9472.334466
        [00:16:04] Bootstrap tree #189, logLikelihood: -9176.173533
        [00:16:11] Bootstrap tree #190, logLikelihood: -9499.321086
        [00:16:16] Bootstrap tree #191, logLikelihood: -9512.848690
        [00:16:19] Bootstrap tree #192, logLikelihood: -9032.214685
        [00:16:24] Bootstrap tree #193, logLikelihood: -9065.123124
        [00:16:27] Bootstrap tree #194, logLikelihood: -9382.634542
        [00:16:31] Bootstrap tree #195, logLikelihood: -9752.310995
        [00:16:35] Bootstrap tree #196, logLikelihood: -9949.842930
        [00:16:39] Bootstrap tree #197, logLikelihood: -9227.400457
        [00:16:44] Bootstrap tree #198, logLikelihood: -9274.634156
        [00:16:48] Bootstrap tree #199, logLikelihood: -9432.199561
        [00:16:52] Bootstrap tree #200, logLikelihood: -8957.253612
        [00:16:56] Bootstrap tree #201, logLikelihood: -9748.946146
        [00:17:01] Bootstrap tree #202, logLikelihood: -9047.983815
        [00:17:05] Bootstrap tree #203, logLikelihood: -9298.651129
        [00:17:09] Bootstrap tree #204, logLikelihood: -9389.286047
        [00:17:13] Bootstrap tree #205, logLikelihood: -9638.478991
        [00:17:17] Bootstrap tree #206, logLikelihood: -9401.653464
        [00:17:21] Bootstrap tree #207, logLikelihood: -9446.388908
        [00:17:25] Bootstrap tree #208, logLikelihood: -9285.272745
        [00:17:29] Bootstrap tree #209, logLikelihood: -9083.672453
        [00:17:36] Bootstrap tree #210, logLikelihood: -9574.533706
        [00:17:40] Bootstrap tree #211, logLikelihood: -9378.655330
        [00:17:44] Bootstrap tree #212, logLikelihood: -9491.756620
        [00:17:47] Bootstrap tree #213, logLikelihood: -9098.993511
        [00:17:50] Bootstrap tree #214, logLikelihood: -9161.261421
        [00:17:54] Bootstrap tree #215, logLikelihood: -9161.430315
        [00:17:58] Bootstrap tree #216, logLikelihood: -9104.111386
        [00:18:04] Bootstrap tree #217, logLikelihood: -9544.767987
        [00:18:07] Bootstrap tree #218, logLikelihood: -9474.673714
        [00:18:11] Bootstrap tree #219, logLikelihood: -9283.021423
        [00:18:15] Bootstrap tree #220, logLikelihood: -9569.987014
        [00:18:19] Bootstrap tree #221, logLikelihood: -9540.041118
        [00:18:24] Bootstrap tree #222, logLikelihood: -9744.836478
        [00:18:28] Bootstrap tree #223, logLikelihood: -8993.257689
        [00:18:32] Bootstrap tree #224, logLikelihood: -9477.499519
        [00:18:37] Bootstrap tree #225, logLikelihood: -9636.901843
        [00:18:42] Bootstrap tree #226, logLikelihood: -9472.109035
        [00:18:46] Bootstrap tree #227, logLikelihood: -9937.900726
        [00:18:49] Bootstrap tree #228, logLikelihood: -9564.176941
        [00:18:53] Bootstrap tree #229, logLikelihood: -9630.219553
        [00:18:58] Bootstrap tree #230, logLikelihood: -9032.963355
        [00:19:02] Bootstrap tree #231, logLikelihood: -9424.603049
        [00:19:05] Bootstrap tree #232, logLikelihood: -9783.870409
        [00:19:09] Bootstrap tree #233, logLikelihood: -9645.913624
        [00:19:14] Bootstrap tree #234, logLikelihood: -9339.298403
        [00:19:18] Bootstrap tree #235, logLikelihood: -9790.277802
        [00:19:21] Bootstrap tree #236, logLikelihood: -9495.261534
        [00:19:27] Bootstrap tree #237, logLikelihood: -9904.167587
        [00:19:30] Bootstrap tree #238, logLikelihood: -9342.078199
        [00:19:33] Bootstrap tree #239, logLikelihood: -9373.230586
        [00:19:37] Bootstrap tree #240, logLikelihood: -9373.438462
        [00:19:41] Bootstrap tree #241, logLikelihood: -9121.946951
        [00:19:45] Bootstrap tree #242, logLikelihood: -9363.171863
        [00:19:51] Bootstrap tree #243, logLikelihood: -10090.838411
        [00:19:56] Bootstrap tree #244, logLikelihood: -9290.030136
        [00:20:01] Bootstrap tree #245, logLikelihood: -9679.893874
        [00:20:05] Bootstrap tree #246, logLikelihood: -9478.796574
        [00:20:08] Bootstrap tree #247, logLikelihood: -8858.207082
        [00:20:13] Bootstrap tree #248, logLikelihood: -9517.975241
        [00:20:17] Bootstrap tree #249, logLikelihood: -9413.355184
        [00:20:20] Bootstrap tree #250, logLikelihood: -9552.160285
        [00:20:24] Bootstrap tree #251, logLikelihood: -9954.181444
        [00:20:31] Bootstrap tree #252, logLikelihood: -9563.330134
        [00:20:35] Bootstrap tree #253, logLikelihood: -9668.304909
        [00:20:38] Bootstrap tree #254, logLikelihood: -9510.829426
        [00:20:42] Bootstrap tree #255, logLikelihood: -8999.777785
        [00:20:47] Bootstrap tree #256, logLikelihood: -9250.827469
        [00:20:51] Bootstrap tree #257, logLikelihood: -9536.928314
        [00:20:57] Bootstrap tree #258, logLikelihood: -9140.636850
        [00:21:01] Bootstrap tree #259, logLikelihood: -9544.803208
        [00:21:05] Bootstrap tree #260, logLikelihood: -9538.047134
        [00:21:11] Bootstrap tree #261, logLikelihood: -9847.249647
        [00:21:14] Bootstrap tree #262, logLikelihood: -9778.144736
        [00:21:19] Bootstrap tree #263, logLikelihood: -9331.889608
        [00:21:22] Bootstrap tree #264, logLikelihood: -9422.000367
        [00:21:26] Bootstrap tree #265, logLikelihood: -9031.866401
        [00:21:30] Bootstrap tree #266, logLikelihood: -9471.956009
        [00:21:35] Bootstrap tree #267, logLikelihood: -9126.879507
        [00:21:41] Bootstrap tree #268, logLikelihood: -9350.937929
        [00:21:44] Bootstrap tree #269, logLikelihood: -9562.118140
        [00:21:49] Bootstrap tree #270, logLikelihood: -9350.806473
        [00:21:52] Bootstrap tree #271, logLikelihood: -9379.901856
        [00:21:56] Bootstrap tree #272, logLikelihood: -9527.561200
        [00:21:59] Bootstrap tree #273, logLikelihood: -9555.737271
        [00:22:03] Bootstrap tree #274, logLikelihood: -9124.257842
        [00:22:07] Bootstrap tree #275, logLikelihood: -9019.617290
        [00:22:12] Bootstrap tree #276, logLikelihood: -9118.747187
        [00:22:18] Bootstrap tree #277, logLikelihood: -9763.349158
        [00:22:21] Bootstrap tree #278, logLikelihood: -9349.644707
        [00:22:25] Bootstrap tree #279, logLikelihood: -9544.768091
        [00:22:29] Bootstrap tree #280, logLikelihood: -9896.381429
        [00:22:33] Bootstrap tree #281, logLikelihood: -9373.580223
        [00:22:36] Bootstrap tree #282, logLikelihood: -9505.184649
        [00:22:40] Bootstrap tree #283, logLikelihood: -9159.176017
        [00:22:43] Bootstrap tree #284, logLikelihood: -9445.734549
        [00:22:46] Bootstrap tree #285, logLikelihood: -9700.945953
        [00:22:51] Bootstrap tree #286, logLikelihood: -9134.810706
        [00:22:54] Bootstrap tree #287, logLikelihood: -9582.527169
        [00:22:59] Bootstrap tree #288, logLikelihood: -10105.932349
        [00:23:04] Bootstrap tree #289, logLikelihood: -9718.364937
        [00:23:09] Bootstrap tree #290, logLikelihood: -9888.511677
        [00:23:12] Bootstrap tree #291, logLikelihood: -9477.093219
        [00:23:16] Bootstrap tree #292, logLikelihood: -9950.634726
        [00:23:20] Bootstrap tree #293, logLikelihood: -9858.090879
        [00:23:24] Bootstrap tree #294, logLikelihood: -9351.684215
        [00:23:30] Bootstrap tree #295, logLikelihood: -9796.620833
        [00:23:34] Bootstrap tree #296, logLikelihood: -10056.520284
        [00:23:38] Bootstrap tree #297, logLikelihood: -9827.558850
        [00:23:42] Bootstrap tree #298, logLikelihood: -9046.859952
        [00:23:46] Bootstrap tree #299, logLikelihood: -9601.756751
        [00:23:50] Bootstrap tree #300, logLikelihood: -9784.588775
        [00:23:53] Bootstrap tree #301, logLikelihood: -8995.515383
        [00:23:59] Bootstrap tree #302, logLikelihood: -9239.059885
        [00:24:16] Bootstrap tree #303, logLikelihood: -9175.143656
        [00:24:19] Bootstrap tree #304, logLikelihood: -9518.933002
        [00:24:26] Bootstrap tree #305, logLikelihood: -9482.461200
        [00:24:34] Bootstrap tree #306, logLikelihood: -9506.140258
        [00:24:44] Bootstrap tree #307, logLikelihood: -9696.000106
        [00:24:48] Bootstrap tree #308, logLikelihood: -9031.985796
        [00:24:52] Bootstrap tree #309, logLikelihood: -9453.114175
        [00:28:58] Bootstrap tree #310, logLikelihood: -10030.020625
        [00:29:09] Bootstrap tree #311, logLikelihood: -9589.869962
        [00:29:15] Bootstrap tree #312, logLikelihood: -9320.880067
        [00:29:22] Bootstrap tree #313, logLikelihood: -9630.198471
        [00:29:28] Bootstrap tree #314, logLikelihood: -9130.087648
        [00:29:34] Bootstrap tree #315, logLikelihood: -9480.376069
        [00:29:42] Bootstrap tree #316, logLikelihood: -10231.036390
        [00:29:46] Bootstrap tree #317, logLikelihood: -9322.918787
        [00:29:50] Bootstrap tree #318, logLikelihood: -9388.116006
        [00:29:56] Bootstrap tree #319, logLikelihood: -9208.083289
        [00:29:59] Bootstrap tree #320, logLikelihood: -9145.626624
        [00:30:04] Bootstrap tree #321, logLikelihood: -9865.418326
        [00:30:08] Bootstrap tree #322, logLikelihood: -8852.225127
        [00:30:12] Bootstrap tree #323, logLikelihood: -9066.412385
        [00:30:16] Bootstrap tree #324, logLikelihood: -9631.274142
        [00:30:19] Bootstrap tree #325, logLikelihood: -9268.152404
        [00:30:24] Bootstrap tree #326, logLikelihood: -9482.886307
        [00:30:28] Bootstrap tree #327, logLikelihood: -9497.976614
        [00:30:32] Bootstrap tree #328, logLikelihood: -9357.569144
        [00:30:36] Bootstrap tree #329, logLikelihood: -9398.562864
        [00:30:41] Bootstrap tree #330, logLikelihood: -9783.530166
        [00:30:46] Bootstrap tree #331, logLikelihood: -9513.186167
        [00:30:49] Bootstrap tree #332, logLikelihood: -9129.467096
        [00:30:53] Bootstrap tree #333, logLikelihood: -9641.198520
        [00:30:56] Bootstrap tree #334, logLikelihood: -9717.851505
        [00:31:01] Bootstrap tree #335, logLikelihood: -9415.123254
        [00:31:04] Bootstrap tree #336, logLikelihood: -10288.719656
        [00:31:11] Bootstrap tree #337, logLikelihood: -9089.737734
        [00:31:15] Bootstrap tree #338, logLikelihood: -9590.446132
        [00:31:19] Bootstrap tree #339, logLikelihood: -9312.292059
        [00:31:24] Bootstrap tree #340, logLikelihood: -9228.348746
        [00:31:27] Bootstrap tree #341, logLikelihood: -9557.712410
        [00:31:32] Bootstrap tree #342, logLikelihood: -9306.401392
        [00:31:36] Bootstrap tree #343, logLikelihood: -9396.068268
        [00:31:40] Bootstrap tree #344, logLikelihood: -9737.283239
        [00:31:44] Bootstrap tree #345, logLikelihood: -9696.769891
        [00:31:49] Bootstrap tree #346, logLikelihood: -9630.065616
        [00:31:54] Bootstrap tree #347, logLikelihood: -9440.185726
        [00:32:00] Bootstrap tree #348, logLikelihood: -9367.853394
        [00:32:05] Bootstrap tree #349, logLikelihood: -9455.981838
        [00:32:10] Bootstrap tree #350, logLikelihood: -9240.980702
        [00:32:15] Bootstrap tree #351, logLikelihood: -9847.183458
        [00:32:19] Bootstrap tree #352, logLikelihood: -9359.572259
        [00:32:23] Bootstrap tree #353, logLikelihood: -9761.924254
        [00:32:27] Bootstrap tree #354, logLikelihood: -9353.965610
        [00:32:31] Bootstrap tree #355, logLikelihood: -9351.634096
        [00:32:34] Bootstrap tree #356, logLikelihood: -9582.454962
        [00:32:38] Bootstrap tree #357, logLikelihood: -9248.605604
        [00:32:43] Bootstrap tree #358, logLikelihood: -9512.823623
        [00:32:46] Bootstrap tree #359, logLikelihood: -9588.332995
        [00:32:49] Bootstrap tree #360, logLikelihood: -9291.352959
        [00:32:53] Bootstrap tree #361, logLikelihood: -9671.466349
        [00:32:57] Bootstrap tree #362, logLikelihood: -9696.994067
        [00:33:00] Bootstrap tree #363, logLikelihood: -9506.344851
        [00:33:04] Bootstrap tree #364, logLikelihood: -8904.809202
        [00:33:08] Bootstrap tree #365, logLikelihood: -9498.421532
        [00:33:14] Bootstrap tree #366, logLikelihood: -9759.148387
        [00:33:19] Bootstrap tree #367, logLikelihood: -9584.683181
        [00:33:23] Bootstrap tree #368, logLikelihood: -9474.311790
        [00:33:28] Bootstrap tree #369, logLikelihood: -9114.294468
        [00:33:34] Bootstrap tree #370, logLikelihood: -9729.351223
        [00:33:38] Bootstrap tree #371, logLikelihood: -9731.758861
        [00:33:41] Bootstrap tree #372, logLikelihood: -9601.354720
        [00:33:45] Bootstrap tree #373, logLikelihood: -9643.231179
        [00:33:49] Bootstrap tree #374, logLikelihood: -9540.801734
        [00:33:53] Bootstrap tree #375, logLikelihood: -9299.121876
        [00:33:58] Bootstrap tree #376, logLikelihood: -9818.376701
        [00:34:01] Bootstrap tree #377, logLikelihood: -8937.231669
        [00:34:06] Bootstrap tree #378, logLikelihood: -9725.750751
        [00:34:10] Bootstrap tree #379, logLikelihood: -9454.995708
        [00:34:14] Bootstrap tree #380, logLikelihood: -9358.422101
        [00:34:17] Bootstrap tree #381, logLikelihood: -9559.481382
        [00:34:21] Bootstrap tree #382, logLikelihood: -9634.232482
        [00:34:25] Bootstrap tree #383, logLikelihood: -9606.169772
        [00:34:30] Bootstrap tree #384, logLikelihood: -9034.231476
        [00:34:33] Bootstrap tree #385, logLikelihood: -10220.575482
        [00:34:38] Bootstrap tree #386, logLikelihood: -9386.365183
        [00:34:43] Bootstrap tree #387, logLikelihood: -9560.139005
        [00:34:47] Bootstrap tree #388, logLikelihood: -9066.134365
        [00:34:53] Bootstrap tree #389, logLikelihood: -9714.327676
        [00:34:56] Bootstrap tree #390, logLikelihood: -9145.034676
        [00:35:00] Bootstrap tree #391, logLikelihood: -9574.564852
        [00:35:04] Bootstrap tree #392, logLikelihood: -8642.623509
        [00:35:08] Bootstrap tree #393, logLikelihood: -9379.708292
        [00:35:12] Bootstrap tree #394, logLikelihood: -9544.425475
        [00:35:16] Bootstrap tree #395, logLikelihood: -9468.155532
        [00:35:21] Bootstrap tree #396, logLikelihood: -9110.861364
        [00:35:26] Bootstrap tree #397, logLikelihood: -9292.733672
        [00:35:30] Bootstrap tree #398, logLikelihood: -9547.577819
        [00:35:34] Bootstrap tree #399, logLikelihood: -9742.596333
        [00:35:38] Bootstrap tree #400, logLikelihood: -9248.305640
        [00:35:38] Bootstrapping converged after 400 replicates.

        Optimized model parameters:

        Partition 0: noname
        Rate heterogeneity: GAMMA (4 cats, mean),  alpha: 1.800587 (ML),  weights&rates: (0.250000,0.268250) (0.250000,0.631982) (0.250000,1.064093) (0.250000,2.035674) 
        Base frequencies (model): 0.076748 0.051691 0.042645 0.051544 0.019803 0.040752 0.061830 0.073152 0.022944 0.053761 0.091904 0.058676 0.023826 0.040126 0.050901 0.068765 0.058565 0.014261 0.032102 0.066005 
        Substitution rates (model): 58.000000 54.000000 81.000000 56.000000 57.000000 105.000000 179.000000 27.000000 36.000000 30.000000 35.000000 54.000000 15.000000 194.000000 378.000000 475.000000 9.000000 11.000000 298.000000 45.000000 16.000000 113.000000 310.000000 29.000000 137.000000 328.000000 22.000000 38.000000 646.000000 44.000000 5.000000 74.000000 101.000000 64.000000 126.000000 20.000000 17.000000 528.000000 34.000000 86.000000 58.000000 81.000000 391.000000 47.000000 12.000000 263.000000 30.000000 10.000000 15.000000 503.000000 232.000000 8.000000 70.000000 16.000000 10.000000 49.000000 767.000000 130.000000 112.000000 11.000000 7.000000 26.000000 15.000000 4.000000 15.000000 59.000000 38.000000 4.000000 46.000000 31.000000 9.000000 5.000000 59.000000 69.000000 17.000000 23.000000 7.000000 31.000000 78.000000 14.000000 223.000000 42.000000 115.000000 209.000000 62.000000 323.000000 26.000000 597.000000 9.000000 72.000000 292.000000 43.000000 4.000000 164.000000 53.000000 51.000000 18.000000 24.000000 20.000000 119.000000 26.000000 12.000000 9.000000 181.000000 18.000000 5.000000 18.000000 30.000000 32.000000 10.000000 7.000000 45.000000 23.000000 6.000000 6.000000 27.000000 14.000000 5.000000 24.000000 201.000000 33.000000 55.000000 8.000000 47.000000 16.000000 56.000000 45.000000 33.000000 40.000000 115.000000 73.000000 46.000000 8.000000 573.000000 11.000000 229.000000 21.000000 479.000000 89.000000 10.000000 40.000000 245.000000 9.000000 32.000000 961.000000 14.000000 388.000000 248.000000 102.000000 59.000000 25.000000 52.000000 24.000000 180.000000 65.000000 4.000000 21.000000 47.000000 103.000000 10.000000 8.000000 14.000000 43.000000 16.000000 29.000000 226.000000 24.000000 18.000000 323.000000 17.000000 92.000000 12.000000 53.000000 536.000000 62.000000 285.000000 118.000000 6.000000 10.000000 23.000000 477.000000 35.000000 63.000000 38.000000 12.000000 21.000000 112.000000 71.000000 25.000000 16.000000 


        Final LogLikelihood: -9471.963973

        AIC score: 18995.927945 / AICc score: 18997.629764 / BIC score: 19119.365195
        Free parameters (model + branch lengths): 26

        Best ML tree saved to: /Users/melettedevore/Desktop/PP563/project/T4.raxml.bestTree
        All ML trees saved to: /Users/melettedevore/Desktop/PP563/project/T4.raxml.mlTrees
        Best ML tree with Felsenstein bootstrap (FBP) support values saved to: /Users/melettedevore/Desktop/PP563/project/T4.raxml.support
        Optimized model saved to: /Users/melettedevore/Desktop/PP563/project/T4.raxml.bestModel
        Bootstrap trees saved to: /Users/melettedevore/Desktop/PP563/project/T4.raxml.bootstraps

        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/T4.raxml.log

        Analysis started: 04-May-2025 14:51:53 / finished: 04-May-2025 15:27:32

        Elapsed time: 2138.955 seconds



    Ran RaxML on MAFFT Alignment:

        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --check --msa mafft-aligned-seqs.fasta --model JTT+G --prefix pM1

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 03-May-2025 18:04:38 as follows:

        raxml-ng --check --msa mafft-aligned-seqs.fasta --model JTT+G --prefix pM1

        Analysis options:
        run mode: Alignment validation
        start tree(s): 
        random seed: 1746313478
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: mafft-aligned-seqs.fasta
        [00:00:00] Loaded alignment with 14 taxa and 760 sites

        Alignment comprises 1 partitions and 760 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 760 / 760
        Gaps: 61.31 %
        Invariant sites: 27.11 %


        Alignment can be successfully read by RAxML-NG.


        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/pM1.raxml.log

        Analysis started: 03-May-2025 18:04:38 / finished: 03-May-2025 18:04:38

        Elapsed time: 0.004 seconds

        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --parse --msa mafft-aligned-seqs.fasta --model JTT+G --prefix pM2

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 03-May-2025 18:05:10 as follows:

        raxml-ng --parse --msa mafft-aligned-seqs.fasta --model JTT+G --prefix pM2

        Analysis options:
        run mode: Alignment parsing and compression
        start tree(s): 
        random seed: 1746313510
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: mafft-aligned-seqs.fasta
        [00:00:00] Loaded alignment with 14 taxa and 760 sites

        Alignment comprises 1 partitions and 639 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 760 / 639
        Gaps: 61.31 %
        Invariant sites: 27.11 %


        NOTE: Binary MSA file created: pM2.raxml.rba

        * Estimated memory requirements                : 11 MB

        * Recommended number of threads / MPI processes: 8

        Please note that numbers given above are rough estimates only. 
        Actual memory consumption and parallel performance on your system may differ!

        Alignment can be successfully read by RAxML-NG.


        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/pM2.raxml.log

        Analysis started: 03-May-2025 18:05:10 / finished: 03-May-2025 18:05:10

        Elapsed time: 0.006 seconds

        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --msa mafft-aligned-seqs.fasta --model JTT+G --prefix M3 --threads 2 --seed 2

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 03-May-2025 18:06:28 as follows:

        raxml-ng --msa mafft-aligned-seqs.fasta --model JTT+G --prefix M3 --threads 2 --seed 2

        Analysis options:
        run mode: ML tree search
        start tree(s): random (10) + parsimony (10)
        random seed: 2
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        fast spr radius: AUTO
        spr subtree cutoff: 1.000000
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: mafft-aligned-seqs.fasta
        [00:00:00] Loaded alignment with 14 taxa and 760 sites

        Alignment comprises 1 partitions and 639 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 760 / 639
        Gaps: 61.31 %
        Invariant sites: 27.11 %


        NOTE: Binary MSA file created: M3.raxml.rba

        [00:00:00] Generating 10 random starting tree(s) with 14 taxa
        [00:00:00] Generating 10 parsimony starting tree(s) with 14 taxa
        [00:00:00] Data distribution: max. partitions/sites/weight per thread: 1 / 320 / 25600

        Starting ML tree search with 20 distinct starting trees

        [00:00:00 -10733.144249] Initial branch length optimization
        [00:00:00 -9519.558844] Model parameter optimization (eps = 10.000000)
        [00:00:00 -9514.394527] AUTODETECT spr round 1 (radius: 5)
        [00:00:00 -9274.662928] AUTODETECT spr round 2 (radius: 10)
        [00:00:00 -9274.647285] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:00 -9274.647285] Model parameter optimization (eps = 3.000000)
        [00:00:00 -9262.700525] FAST spr round 1 (radius: 5)
        [00:00:01 -9257.439985] FAST spr round 2 (radius: 5)
        [00:00:02 -9255.257702] FAST spr round 3 (radius: 5)
        [00:00:02 -9254.817912] FAST spr round 4 (radius: 5)
        [00:00:03 -9254.817520] Model parameter optimization (eps = 1.000000)
        [00:00:03 -9254.510521] SLOW spr round 1 (radius: 5)
        [00:00:05 -9254.502811] SLOW spr round 2 (radius: 10)
        [00:00:06 -9254.502762] Model parameter optimization (eps = 0.100000)

        [00:00:06] ML tree search #1, logLikelihood: -9254.469284

        [00:00:06 -10785.428012] Initial branch length optimization
        [00:00:06 -9533.730980] Model parameter optimization (eps = 10.000000)
        [00:00:06 -9529.998862] AUTODETECT spr round 1 (radius: 5)
        [00:00:06 -9272.810831] AUTODETECT spr round 2 (radius: 10)
        [00:00:06 -9272.789857] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:06 -9272.789857] Model parameter optimization (eps = 3.000000)
        [00:00:06 -9260.116372] FAST spr round 1 (radius: 5)
        [00:00:07 -9254.872327] FAST spr round 2 (radius: 5)
        [00:00:08 -9254.872249] Model parameter optimization (eps = 1.000000)
        [00:00:08 -9254.517722] SLOW spr round 1 (radius: 5)
        [00:00:09 -9254.509735] SLOW spr round 2 (radius: 10)
        [00:00:10 -9254.509697] Model parameter optimization (eps = 0.100000)

        [00:00:10] ML tree search #2, logLikelihood: -9254.469948

        [00:00:10 -10606.163298] Initial branch length optimization
        [00:00:10 -9473.720752] Model parameter optimization (eps = 10.000000)
        [00:00:10 -9467.933948] AUTODETECT spr round 1 (radius: 5)
        [00:00:10 -9289.919122] AUTODETECT spr round 2 (radius: 10)
        [00:00:11 -9289.909700] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:11 -9289.909700] Model parameter optimization (eps = 3.000000)
        [00:00:11 -9279.807036] FAST spr round 1 (radius: 5)
        [00:00:11 -9255.003701] FAST spr round 2 (radius: 5)
        [00:00:12 -9255.003610] Model parameter optimization (eps = 1.000000)
        [00:00:12 -9254.533692] SLOW spr round 1 (radius: 5)
        [00:00:14 -9254.525331] SLOW spr round 2 (radius: 10)
        [00:00:15 -9254.525217] Model parameter optimization (eps = 0.100000)

        [00:00:15] ML tree search #3, logLikelihood: -9254.471990

        [00:00:15 -10865.716468] Initial branch length optimization
        [00:00:15 -9534.246069] Model parameter optimization (eps = 10.000000)
        [00:00:15 -9529.334089] AUTODETECT spr round 1 (radius: 5)
        [00:00:15 -9294.035465] AUTODETECT spr round 2 (radius: 10)
        [00:00:15 -9294.027284] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:15 -9294.027284] Model parameter optimization (eps = 3.000000)
        [00:00:15 -9282.535937] FAST spr round 1 (radius: 5)
        [00:00:16 -9256.094570] FAST spr round 2 (radius: 5)
        [00:00:17 -9255.224274] FAST spr round 3 (radius: 5)
        [00:00:18 -9255.224067] Model parameter optimization (eps = 1.000000)
        [00:00:18 -9254.570771] SLOW spr round 1 (radius: 5)
        [00:00:19 -9254.551645] SLOW spr round 2 (radius: 10)
        [00:00:20 -9254.551494] Model parameter optimization (eps = 0.100000)

        [00:00:20] ML tree search #4, logLikelihood: -9254.475545

        [00:00:20 -10655.666327] Initial branch length optimization
        [00:00:20 -9506.307425] Model parameter optimization (eps = 10.000000)
        [00:00:20 -9501.993037] AUTODETECT spr round 1 (radius: 5)
        [00:00:20 -9281.093778] AUTODETECT spr round 2 (radius: 10)
        [00:00:21 -9281.084085] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:21 -9281.084085] Model parameter optimization (eps = 3.000000)
        [00:00:21 -9269.934187] FAST spr round 1 (radius: 5)
        [00:00:21 -9255.406834] FAST spr round 2 (radius: 5)
        [00:00:22 -9254.966848] FAST spr round 3 (radius: 5)
        [00:00:23 -9254.966545] Model parameter optimization (eps = 1.000000)
        [00:00:23 -9254.531096] SLOW spr round 1 (radius: 5)
        [00:00:25 -9254.520944] SLOW spr round 2 (radius: 10)
        [00:00:26 -9254.520832] Model parameter optimization (eps = 0.100000)

        [00:00:26] ML tree search #5, logLikelihood: -9254.471678

        [00:00:26 -10923.891062] Initial branch length optimization
        [00:00:26 -9564.512160] Model parameter optimization (eps = 10.000000)
        [00:00:26 -9563.211397] AUTODETECT spr round 1 (radius: 5)
        [00:00:26 -9289.194646] AUTODETECT spr round 2 (radius: 10)
        [00:00:26 -9289.183324] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:26 -9289.183324] Model parameter optimization (eps = 3.000000)
        [00:00:26 -9273.275937] FAST spr round 1 (radius: 5)
        [00:00:27 -9256.101218] FAST spr round 2 (radius: 5)
        [00:00:28 -9255.582103] FAST spr round 3 (radius: 5)
        [00:00:29 -9255.142268] FAST spr round 4 (radius: 5)
        [00:00:29 -9255.141949] Model parameter optimization (eps = 1.000000)
        [00:00:29 -9254.556697] SLOW spr round 1 (radius: 5)
        [00:00:31 -9254.541621] SLOW spr round 2 (radius: 10)
        [00:00:32 -9254.541416] Model parameter optimization (eps = 0.100000)

        [00:00:32] ML tree search #6, logLikelihood: -9254.473591

        [00:00:32 -10834.997706] Initial branch length optimization
        [00:00:32 -9533.244042] Model parameter optimization (eps = 10.000000)
        [00:00:32 -9529.487895] AUTODETECT spr round 1 (radius: 5)
        [00:00:32 -9287.892390] AUTODETECT spr round 2 (radius: 10)
        [00:00:32 -9287.883667] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:32 -9287.883667] Model parameter optimization (eps = 3.000000)
        [00:00:33 -9274.341718] FAST spr round 1 (radius: 5)
        [00:00:33 -9256.759206] FAST spr round 2 (radius: 5)
        [00:00:34 -9255.161527] FAST spr round 3 (radius: 5)
        [00:00:35 -9255.161488] Model parameter optimization (eps = 1.000000)
        [00:00:35 -9254.554777] SLOW spr round 1 (radius: 5)
        [00:00:37 -9254.543899] SLOW spr round 2 (radius: 10)
        [00:00:38 -9254.543795] Model parameter optimization (eps = 0.100000)

        [00:00:38] ML tree search #7, logLikelihood: -9254.473402

        [00:00:38 -10650.632356] Initial branch length optimization
        [00:00:38 -9463.381074] Model parameter optimization (eps = 10.000000)
        [00:00:38 -9454.321962] AUTODETECT spr round 1 (radius: 5)
        [00:00:38 -9270.185097] AUTODETECT spr round 2 (radius: 10)
        [00:00:38 -9270.180815] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:38 -9270.180815] Model parameter optimization (eps = 3.000000)
        [00:00:38 -9262.391627] FAST spr round 1 (radius: 5)
        [00:00:39 -9255.565598] FAST spr round 2 (radius: 5)
        [00:00:40 -9254.774914] FAST spr round 3 (radius: 5)
        [00:00:41 -9254.773947] Model parameter optimization (eps = 1.000000)
        [00:00:41 -9254.503991] SLOW spr round 1 (radius: 5)
        [00:00:42 -9254.498002] SLOW spr round 2 (radius: 10)
        [00:00:43 -9254.497948] Model parameter optimization (eps = 0.100000)

        [00:00:43] ML tree search #8, logLikelihood: -9254.468644

        [00:00:43 -10880.613511] Initial branch length optimization
        [00:00:43 -9537.424647] Model parameter optimization (eps = 10.000000)
        [00:00:43 -9534.109044] AUTODETECT spr round 1 (radius: 5)
        [00:00:44 -9276.652329] AUTODETECT spr round 2 (radius: 10)
        [00:00:44 -9276.645594] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:44 -9276.645594] Model parameter optimization (eps = 3.000000)
        [00:00:44 -9264.315606] FAST spr round 1 (radius: 5)
        [00:00:45 -9255.460845] FAST spr round 2 (radius: 5)
        [00:00:45 -9255.009169] FAST spr round 3 (radius: 5)
        [00:00:46 -9255.008823] Model parameter optimization (eps = 1.000000)
        [00:00:46 -9254.535938] SLOW spr round 1 (radius: 5)
        [00:00:48 -9254.525722] SLOW spr round 2 (radius: 10)
        [00:00:49 -9254.525622] Model parameter optimization (eps = 0.100000)

        [00:00:49] ML tree search #9, logLikelihood: -9254.472174

        [00:00:49 -10615.858606] Initial branch length optimization
        [00:00:49 -9527.449568] Model parameter optimization (eps = 10.000000)
        [00:00:49 -9524.162974] AUTODETECT spr round 1 (radius: 5)
        [00:00:49 -9277.798689] AUTODETECT spr round 2 (radius: 10)
        [00:00:49 -9277.790637] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:49 -9277.790637] Model parameter optimization (eps = 3.000000)
        [00:00:49 -9264.117021] FAST spr round 1 (radius: 5)
        [00:00:50 -9255.081305] FAST spr round 2 (radius: 5)
        [00:00:51 -9255.081018] Model parameter optimization (eps = 1.000000)
        [00:00:51 -9254.544293] SLOW spr round 1 (radius: 5)
        [00:00:53 -9254.534029] SLOW spr round 2 (radius: 10)
        [00:00:53 -9254.533971] Model parameter optimization (eps = 0.100000)

        [00:00:54] ML tree search #10, logLikelihood: -9254.472887

        [00:00:54 -10105.928987] Initial branch length optimization
        [00:00:54 -9306.938435] Model parameter optimization (eps = 10.000000)
        [00:00:54 -9285.398024] AUTODETECT spr round 1 (radius: 5)
        [00:00:54 -9257.242791] AUTODETECT spr round 2 (radius: 10)
        [00:00:54 -9257.240124] SPR radius for FAST iterations: 5 (autodetect)
        [00:00:54 -9257.240124] Model parameter optimization (eps = 3.000000)
        [00:00:54 -9255.952153] FAST spr round 1 (radius: 5)
        [00:00:55 -9255.066774] FAST spr round 2 (radius: 5)
        [00:00:56 -9254.626609] FAST spr round 3 (radius: 5)
        [00:00:56 -9254.626515] Model parameter optimization (eps = 1.000000)
        [00:00:57 -9254.485332] SLOW spr round 1 (radius: 5)
        [00:00:58 -9254.481890] SLOW spr round 2 (radius: 10)
        [00:00:59 -9254.481873] Model parameter optimization (eps = 0.100000)

        [00:00:59] ML tree search #11, logLikelihood: -9254.466824

        [00:00:59 -10125.459109] Initial branch length optimization
        [00:00:59 -9303.609795] Model parameter optimization (eps = 10.000000)
        [00:00:59 -9282.717579] AUTODETECT spr round 1 (radius: 5)
        [00:01:00 -9257.490844] AUTODETECT spr round 2 (radius: 10)
        [00:01:00 -9257.486627] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:00 -9257.486627] Model parameter optimization (eps = 3.000000)
        [00:01:00 -9256.188441] FAST spr round 1 (radius: 5)
        [00:01:01 -9254.605526] FAST spr round 2 (radius: 5)
        [00:01:01 -9254.605401] Model parameter optimization (eps = 1.000000)
        [00:01:01 -9254.482641] SLOW spr round 1 (radius: 5)
        [00:01:03 -9254.479657] SLOW spr round 2 (radius: 10)
        [00:01:04 -9254.479631] Model parameter optimization (eps = 0.100000)

        [00:01:04] ML tree search #12, logLikelihood: -9254.466603

        [00:01:04 -10056.467248] Initial branch length optimization
        [00:01:04 -9284.197567] Model parameter optimization (eps = 10.000000)
        [00:01:04 -9261.543093] AUTODETECT spr round 1 (radius: 5)
        [00:01:04 -9255.282498] AUTODETECT spr round 2 (radius: 10)
        [00:01:04 -9255.281169] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:04 -9255.281169] Model parameter optimization (eps = 3.000000)
        [00:01:05 -9254.574195] FAST spr round 1 (radius: 5)
        [00:01:05 -9254.559729] Model parameter optimization (eps = 1.000000)
        [00:01:05 -9254.476519] SLOW spr round 1 (radius: 5)
        [00:01:07 -9254.474903] SLOW spr round 2 (radius: 10)
        [00:01:08 -9254.474892] Model parameter optimization (eps = 0.100000)

        [00:01:08] ML tree search #13, logLikelihood: -9254.466015

        [00:01:08 -10073.790112] Initial branch length optimization
        [00:01:08 -9287.907801] Model parameter optimization (eps = 10.000000)
        [00:01:08 -9267.173232] AUTODETECT spr round 1 (radius: 5)
        [00:01:08 -9256.721350] AUTODETECT spr round 2 (radius: 10)
        [00:01:08 -9256.720314] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:08 -9256.720314] Model parameter optimization (eps = 3.000000)
        [00:01:08 -9255.167780] FAST spr round 1 (radius: 5)
        [00:01:09 -9254.687794] FAST spr round 2 (radius: 5)
        [00:01:10 -9254.687273] Model parameter optimization (eps = 1.000000)
        [00:01:10 -9254.492898] SLOW spr round 1 (radius: 5)
        [00:01:12 -9254.488480] SLOW spr round 2 (radius: 10)
        [00:01:13 -9254.488417] Model parameter optimization (eps = 0.100000)

        [00:01:13] ML tree search #14, logLikelihood: -9254.467841

        [00:01:13 -10108.337160] Initial branch length optimization
        [00:01:13 -9292.376712] Model parameter optimization (eps = 10.000000)
        [00:01:13 -9271.268342] AUTODETECT spr round 1 (radius: 5)
        [00:01:13 -9256.032577] AUTODETECT spr round 2 (radius: 10)
        [00:01:13 -9256.030056] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:13 -9256.030056] Model parameter optimization (eps = 3.000000)
        [00:01:13 -9254.687762] FAST spr round 1 (radius: 5)
        [00:01:14 -9254.656523] Model parameter optimization (eps = 1.000000)
        [00:01:14 -9254.487826] SLOW spr round 1 (radius: 5)
        [00:01:15 -9254.484071] SLOW spr round 2 (radius: 10)
        [00:01:16 -9254.484034] Model parameter optimization (eps = 0.100000)

        [00:01:16] ML tree search #15, logLikelihood: -9254.467177

        [00:01:16 -10095.793018] Initial branch length optimization
        [00:01:16 -9306.401300] Model parameter optimization (eps = 10.000000)
        [00:01:17 -9284.788951] AUTODETECT spr round 1 (radius: 5)
        [00:01:17 -9256.395576] AUTODETECT spr round 2 (radius: 10)
        [00:01:17 -9256.392243] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:17 -9256.392243] Model parameter optimization (eps = 3.000000)
        [00:01:17 -9255.132947] FAST spr round 1 (radius: 5)
        [00:01:18 -9254.585587] FAST spr round 2 (radius: 5)
        [00:01:18 -9254.584807] Model parameter optimization (eps = 1.000000)
        [00:01:18 -9254.480539] SLOW spr round 1 (radius: 5)
        [00:01:20 -9254.477661] SLOW spr round 2 (radius: 10)
        [00:01:21 -9254.477639] Model parameter optimization (eps = 0.100000)

        [00:01:21] ML tree search #16, logLikelihood: -9254.466375

        [00:01:21 -10122.549101] Initial branch length optimization
        [00:01:21 -9328.450984] Model parameter optimization (eps = 10.000000)
        [00:01:21 -9309.194355] AUTODETECT spr round 1 (radius: 5)
        [00:01:21 -9257.817999] AUTODETECT spr round 2 (radius: 10)
        [00:01:21 -9257.816574] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:21 -9257.816574] Model parameter optimization (eps = 3.000000)
        [00:01:21 -9255.939818] FAST spr round 1 (radius: 5)
        [00:01:22 -9255.145561] FAST spr round 2 (radius: 5)
        [00:01:23 -9254.705666] FAST spr round 3 (radius: 5)
        [00:01:24 -9254.705322] Model parameter optimization (eps = 1.000000)
        [00:01:24 -9254.496399] SLOW spr round 1 (radius: 5)
        [00:01:25 -9254.490915] SLOW spr round 2 (radius: 10)
        [00:01:26 -9254.490858] Model parameter optimization (eps = 0.100000)

        [00:01:26] ML tree search #17, logLikelihood: -9254.468123

        [00:01:26 -10081.768397] Initial branch length optimization
        [00:01:26 -9281.302696] Model parameter optimization (eps = 10.000000)
        [00:01:26 -9258.771140] AUTODETECT spr round 1 (radius: 5)
        [00:01:27 -9256.554024] AUTODETECT spr round 2 (radius: 10)
        [00:01:27 -9256.550010] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:27 -9256.550010] Model parameter optimization (eps = 3.000000)
        [00:01:27 -9255.486015] FAST spr round 1 (radius: 5)
        [00:01:28 -9255.003340] FAST spr round 2 (radius: 5)
        [00:01:28 -9254.563501] FAST spr round 3 (radius: 5)
        [00:01:29 -9254.563197] Model parameter optimization (eps = 1.000000)
        [00:01:29 -9254.477456] SLOW spr round 1 (radius: 5)
        [00:01:31 -9254.475292] SLOW spr round 2 (radius: 10)
        [00:01:32 -9254.475264] Model parameter optimization (eps = 0.100000)

        [00:01:32] ML tree search #18, logLikelihood: -9254.466081

        [00:01:32 -10099.994375] Initial branch length optimization
        [00:01:32 -9294.661041] Model parameter optimization (eps = 10.000000)
        [00:01:32 -9274.043541] AUTODETECT spr round 1 (radius: 5)
        [00:01:32 -9256.028014] AUTODETECT spr round 2 (radius: 10)
        [00:01:32 -9256.025950] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:32 -9256.025950] Model parameter optimization (eps = 3.000000)
        [00:01:32 -9254.693210] FAST spr round 1 (radius: 5)
        [00:01:33 -9254.655145] Model parameter optimization (eps = 1.000000)
        [00:01:33 -9254.488768] SLOW spr round 1 (radius: 5)
        [00:01:35 -9254.484106] SLOW spr round 2 (radius: 10)
        [00:01:35 -9254.484068] Model parameter optimization (eps = 0.100000)

        [00:01:35] ML tree search #19, logLikelihood: -9254.467243

        [00:01:35 -10055.623974] Initial branch length optimization
        [00:01:36 -9293.528310] Model parameter optimization (eps = 10.000000)
        [00:01:36 -9271.345011] AUTODETECT spr round 1 (radius: 5)
        [00:01:36 -9256.948464] AUTODETECT spr round 2 (radius: 10)
        [00:01:36 -9256.947284] SPR radius for FAST iterations: 5 (autodetect)
        [00:01:36 -9256.947284] Model parameter optimization (eps = 3.000000)
        [00:01:36 -9255.544195] FAST spr round 1 (radius: 5)
        [00:01:37 -9254.605881] FAST spr round 2 (radius: 5)
        [00:01:38 -9254.605855] Model parameter optimization (eps = 1.000000)
        [00:01:38 -9254.482233] SLOW spr round 1 (radius: 5)
        [00:01:39 -9254.479713] SLOW spr round 2 (radius: 10)
        [00:01:40 -9254.479689] Model parameter optimization (eps = 0.100000)

        [00:01:40] ML tree search #20, logLikelihood: -9254.466599


        Optimized model parameters:

        Partition 0: noname
        Rate heterogeneity: GAMMA (4 cats, mean),  alpha: 2.259584 (ML),  weights&rates: (0.250000,0.322410) (0.250000,0.680054) (0.250000,1.075309) (0.250000,1.922227) 
        Base frequencies (model): 0.076748 0.051691 0.042645 0.051544 0.019803 0.040752 0.061830 0.073152 0.022944 0.053761 0.091904 0.058676 0.023826 0.040126 0.050901 0.068765 0.058565 0.014261 0.032102 0.066005 
        Substitution rates (model): 58.000000 54.000000 81.000000 56.000000 57.000000 105.000000 179.000000 27.000000 36.000000 30.000000 35.000000 54.000000 15.000000 194.000000 378.000000 475.000000 9.000000 11.000000 298.000000 45.000000 16.000000 113.000000 310.000000 29.000000 137.000000 328.000000 22.000000 38.000000 646.000000 44.000000 5.000000 74.000000 101.000000 64.000000 126.000000 20.000000 17.000000 528.000000 34.000000 86.000000 58.000000 81.000000 391.000000 47.000000 12.000000 263.000000 30.000000 10.000000 15.000000 503.000000 232.000000 8.000000 70.000000 16.000000 10.000000 49.000000 767.000000 130.000000 112.000000 11.000000 7.000000 26.000000 15.000000 4.000000 15.000000 59.000000 38.000000 4.000000 46.000000 31.000000 9.000000 5.000000 59.000000 69.000000 17.000000 23.000000 7.000000 31.000000 78.000000 14.000000 223.000000 42.000000 115.000000 209.000000 62.000000 323.000000 26.000000 597.000000 9.000000 72.000000 292.000000 43.000000 4.000000 164.000000 53.000000 51.000000 18.000000 24.000000 20.000000 119.000000 26.000000 12.000000 9.000000 181.000000 18.000000 5.000000 18.000000 30.000000 32.000000 10.000000 7.000000 45.000000 23.000000 6.000000 6.000000 27.000000 14.000000 5.000000 24.000000 201.000000 33.000000 55.000000 8.000000 47.000000 16.000000 56.000000 45.000000 33.000000 40.000000 115.000000 73.000000 46.000000 8.000000 573.000000 11.000000 229.000000 21.000000 479.000000 89.000000 10.000000 40.000000 245.000000 9.000000 32.000000 961.000000 14.000000 388.000000 248.000000 102.000000 59.000000 25.000000 52.000000 24.000000 180.000000 65.000000 4.000000 21.000000 47.000000 103.000000 10.000000 8.000000 14.000000 43.000000 16.000000 29.000000 226.000000 24.000000 18.000000 323.000000 17.000000 92.000000 12.000000 53.000000 536.000000 62.000000 285.000000 118.000000 6.000000 10.000000 23.000000 477.000000 35.000000 63.000000 38.000000 12.000000 21.000000 112.000000 71.000000 25.000000 16.000000 


        Final LogLikelihood: -9254.466015

        AIC score: 18560.932029 / AICc score: 18562.847445 / BIC score: 18681.398308
        Free parameters (model + branch lengths): 26

        Best ML tree saved to: /Users/melettedevore/Desktop/PP563/project/M3.raxml.bestTree
        All ML trees saved to: /Users/melettedevore/Desktop/PP563/project/M3.raxml.mlTrees
        Optimized model saved to: /Users/melettedevore/Desktop/PP563/project/M3.raxml.bestModel

        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/M3.raxml.log

        Analysis started: 03-May-2025 18:06:28 / finished: 03-May-2025 18:08:08

        Elapsed time: 100.707 seconds

    *Re-Ran with bootstrap analysis:

        (base) Melettes-MBP-6:project melettedevore$ raxml-ng --all --msa mafft-aligned-seqs.fasta --model JTT+G --prefix M4 --threads 2 --seed 2

        RAxML-NG v. 0.9.0 released on 20.05.2019 by The Exelixis Lab.
        Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
        Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
        Latest version: https://github.com/amkozlov/raxml-ng
        Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

        RAxML-NG was called at 04-May-2025 13:53:06 as follows:

        raxml-ng --all --msa mafft-aligned-seqs.fasta --model JTT+G --prefix M4 --threads 2 --seed 2

        Analysis options:
        run mode: ML tree search + bootstrapping (Felsenstein Bootstrap)
        start tree(s): random (10) + parsimony (10)
        bootstrap replicates: max: 1000 + bootstopping (autoMRE, cutoff: 0.030000)
        random seed: 2
        tip-inner: OFF
        pattern compression: ON
        per-rate scalers: OFF
        site repeats: ON
        branch lengths: proportional (ML estimate, algorithm: NR-FAST)
        SIMD kernels: AVX2
        parallelization: PTHREADS (2 threads), thread pinning: OFF

        [00:00:00] Reading alignment from file: mafft-aligned-seqs.fasta
        [00:00:00] Loaded alignment with 14 taxa and 760 sites

        Alignment comprises 1 partitions and 639 patterns

        Partition 0: noname
        Model: JTT+G4m
        Alignment sites / patterns: 760 / 639
        Gaps: 61.31 %
        Invariant sites: 27.11 %


        NOTE: Binary MSA file created: M4.raxml.rba

        [00:00:00] Generating 10 random starting tree(s) with 14 taxa
        [00:00:00] Generating 10 parsimony starting tree(s) with 14 taxa
        [00:00:00] Data distribution: max. partitions/sites/weight per thread: 1 / 320 / 25600

        Starting ML tree search with 20 distinct starting trees

        [00:00:11] ML tree search #1, logLikelihood: -9254.469284
        [00:00:19] ML tree search #2, logLikelihood: -9254.469948
        [00:00:26] ML tree search #3, logLikelihood: -9254.471990
        [00:00:38] ML tree search #4, logLikelihood: -9254.475545
        [00:00:47] ML tree search #5, logLikelihood: -9254.471678
        [00:01:01] ML tree search #6, logLikelihood: -9254.473591
        [00:01:09] ML tree search #7, logLikelihood: -9254.473402
        [00:01:17] ML tree search #8, logLikelihood: -9254.468644
        [00:01:24] ML tree search #9, logLikelihood: -9254.472174
        [00:01:34] ML tree search #10, logLikelihood: -9254.472887
        [00:01:43] ML tree search #11, logLikelihood: -9254.466824
        [00:01:48] ML tree search #12, logLikelihood: -9254.466603
        [00:01:52] ML tree search #13, logLikelihood: -9254.466015
        [00:01:56] ML tree search #14, logLikelihood: -9254.467841
        [00:02:00] ML tree search #15, logLikelihood: -9254.467177
        [00:02:04] ML tree search #16, logLikelihood: -9254.466375
        [00:02:09] ML tree search #17, logLikelihood: -9254.468123
        [00:02:14] ML tree search #18, logLikelihood: -9254.466081
        [00:02:18] ML tree search #19, logLikelihood: -9254.467243
        [00:02:22] ML tree search #20, logLikelihood: -9254.466599

        [00:02:22] ML tree search completed, best tree logLH: -9254.466015

        [00:02:22] Starting bootstrapping analysis with 1000 replicates.

        [00:02:25] Bootstrap tree #1, logLikelihood: -9502.318273
        [00:02:29] Bootstrap tree #2, logLikelihood: -9290.002611
        [00:02:32] Bootstrap tree #3, logLikelihood: -9208.351869
        [00:02:35] Bootstrap tree #4, logLikelihood: -9280.348053
        [00:02:39] Bootstrap tree #5, logLikelihood: -9434.734712
        [00:02:42] Bootstrap tree #6, logLikelihood: -9459.970273
        [00:02:45] Bootstrap tree #7, logLikelihood: -9389.296692
        [00:02:49] Bootstrap tree #8, logLikelihood: -9618.498326
        [00:02:54] Bootstrap tree #9, logLikelihood: -9011.852284
        [00:02:57] Bootstrap tree #10, logLikelihood: -9571.643863
        [00:03:01] Bootstrap tree #11, logLikelihood: -9214.554731
        [00:03:04] Bootstrap tree #12, logLikelihood: -9064.777325
        [00:03:07] Bootstrap tree #13, logLikelihood: -9325.459154
        [00:03:11] Bootstrap tree #14, logLikelihood: -9827.590602
        [00:03:13] Bootstrap tree #15, logLikelihood: -9032.087985
        [00:03:17] Bootstrap tree #16, logLikelihood: -9374.340095
        [00:03:22] Bootstrap tree #17, logLikelihood: -9354.555328
        [00:03:24] Bootstrap tree #18, logLikelihood: -9210.871839
        [00:03:28] Bootstrap tree #19, logLikelihood: -9373.940816
        [00:03:32] Bootstrap tree #20, logLikelihood: -9127.192344
        [00:03:34] Bootstrap tree #21, logLikelihood: -9400.965212
        [00:03:40] Bootstrap tree #22, logLikelihood: -9207.883761
        [00:03:44] Bootstrap tree #23, logLikelihood: -9512.551002
        [00:03:47] Bootstrap tree #24, logLikelihood: -9688.198144
        [00:03:51] Bootstrap tree #25, logLikelihood: -9398.481946
        [00:03:55] Bootstrap tree #26, logLikelihood: -9087.992768
        [00:03:59] Bootstrap tree #27, logLikelihood: -9174.901947
        [00:04:02] Bootstrap tree #28, logLikelihood: -9379.173218
        [00:04:06] Bootstrap tree #29, logLikelihood: -9282.971511
        [00:04:10] Bootstrap tree #30, logLikelihood: -9411.894852
        [00:04:13] Bootstrap tree #31, logLikelihood: -9435.361218
        [00:04:17] Bootstrap tree #32, logLikelihood: -9117.843381
        [00:04:21] Bootstrap tree #33, logLikelihood: -9314.479269
        [00:04:24] Bootstrap tree #34, logLikelihood: -9434.387382
        [00:04:28] Bootstrap tree #35, logLikelihood: -9077.039074
        [00:04:31] Bootstrap tree #36, logLikelihood: -9507.917710
        [00:04:36] Bootstrap tree #37, logLikelihood: -9362.618038
        [00:04:40] Bootstrap tree #38, logLikelihood: -8865.077434
        [00:04:44] Bootstrap tree #39, logLikelihood: -9404.820686
        [00:04:48] Bootstrap tree #40, logLikelihood: -9150.721897
        [00:04:53] Bootstrap tree #41, logLikelihood: -9331.309333
        [00:04:57] Bootstrap tree #42, logLikelihood: -8759.799737
        [00:05:02] Bootstrap tree #43, logLikelihood: -8927.334566
        [00:05:07] Bootstrap tree #44, logLikelihood: -9340.523190
        [00:05:11] Bootstrap tree #45, logLikelihood: -9251.165896
        [00:05:16] Bootstrap tree #46, logLikelihood: -8682.618777
        [00:05:20] Bootstrap tree #47, logLikelihood: -9413.101538
        [00:05:25] Bootstrap tree #48, logLikelihood: -9005.659298
        [00:05:28] Bootstrap tree #49, logLikelihood: -9216.717049
        [00:05:32] Bootstrap tree #50, logLikelihood: -9613.537834
        [00:05:36] Bootstrap tree #51, logLikelihood: -8658.396820
        [00:05:39] Bootstrap tree #52, logLikelihood: -9423.804835
        [00:05:44] Bootstrap tree #53, logLikelihood: -9328.681549
        [00:05:47] Bootstrap tree #54, logLikelihood: -9569.649941
        [00:05:50] Bootstrap tree #55, logLikelihood: -9264.704197
        [00:05:54] Bootstrap tree #56, logLikelihood: -9732.383603
        [00:05:57] Bootstrap tree #57, logLikelihood: -8635.550591
        [00:06:00] Bootstrap tree #58, logLikelihood: -8311.348241
        [00:06:03] Bootstrap tree #59, logLikelihood: -9312.088756
        [00:06:08] Bootstrap tree #60, logLikelihood: -9506.026223
        [00:06:12] Bootstrap tree #61, logLikelihood: -9436.607128
        [00:06:16] Bootstrap tree #62, logLikelihood: -9529.224510
        [00:06:19] Bootstrap tree #63, logLikelihood: -9489.231081
        [00:06:22] Bootstrap tree #64, logLikelihood: -9161.708045
        [00:06:26] Bootstrap tree #65, logLikelihood: -9383.259569
        [00:06:30] Bootstrap tree #66, logLikelihood: -9848.013691
        [00:06:34] Bootstrap tree #67, logLikelihood: -9775.688413
        [00:06:38] Bootstrap tree #68, logLikelihood: -9153.494477
        [00:06:41] Bootstrap tree #69, logLikelihood: -8907.518166
        [00:06:45] Bootstrap tree #70, logLikelihood: -8626.515202
        [00:06:48] Bootstrap tree #71, logLikelihood: -8957.324271
        [00:06:51] Bootstrap tree #72, logLikelihood: -9354.379564
        [00:06:55] Bootstrap tree #73, logLikelihood: -9406.138445
        [00:06:58] Bootstrap tree #74, logLikelihood: -9382.465919
        [00:07:01] Bootstrap tree #75, logLikelihood: -9451.151257
        [00:07:04] Bootstrap tree #76, logLikelihood: -8908.974560
        [00:07:07] Bootstrap tree #77, logLikelihood: -9323.017758
        [00:07:11] Bootstrap tree #78, logLikelihood: -9427.523482
        [00:07:14] Bootstrap tree #79, logLikelihood: -9190.785604
        [00:07:18] Bootstrap tree #80, logLikelihood: -9210.765587
        [00:07:22] Bootstrap tree #81, logLikelihood: -9303.430837
        [00:07:26] Bootstrap tree #82, logLikelihood: -9257.101837
        [00:07:30] Bootstrap tree #83, logLikelihood: -8925.269123
        [00:07:33] Bootstrap tree #84, logLikelihood: -8980.933358
        [00:07:36] Bootstrap tree #85, logLikelihood: -9257.259835
        [00:07:41] Bootstrap tree #86, logLikelihood: -8779.479271
        [00:07:46] Bootstrap tree #87, logLikelihood: -9650.630385
        [00:07:49] Bootstrap tree #88, logLikelihood: -9282.697182
        [00:07:52] Bootstrap tree #89, logLikelihood: -9212.856477
        [00:07:56] Bootstrap tree #90, logLikelihood: -8952.483330
        [00:07:58] Bootstrap tree #91, logLikelihood: -9149.130030
        [00:08:03] Bootstrap tree #92, logLikelihood: -9452.895229
        [00:08:07] Bootstrap tree #93, logLikelihood: -9272.491941
        [00:08:12] Bootstrap tree #94, logLikelihood: -8922.165973
        [00:08:17] Bootstrap tree #95, logLikelihood: -9231.821892
        [00:08:20] Bootstrap tree #96, logLikelihood: -9666.150795
        [00:08:24] Bootstrap tree #97, logLikelihood: -9063.769320
        [00:08:28] Bootstrap tree #98, logLikelihood: -9095.864881
        [00:08:34] Bootstrap tree #99, logLikelihood: -8996.746715
        [00:08:37] Bootstrap tree #100, logLikelihood: -9652.362819
        [00:08:40] Bootstrap tree #101, logLikelihood: -8912.552317
        [00:08:44] Bootstrap tree #102, logLikelihood: -9091.523278
        [00:08:47] Bootstrap tree #103, logLikelihood: -9233.262945
        [00:08:50] Bootstrap tree #104, logLikelihood: -9085.966022
        [00:08:54] Bootstrap tree #105, logLikelihood: -8780.065085
        [00:08:59] Bootstrap tree #106, logLikelihood: -8886.063580
        [00:09:03] Bootstrap tree #107, logLikelihood: -9189.204110
        [00:09:07] Bootstrap tree #108, logLikelihood: -9210.829345
        [00:09:11] Bootstrap tree #109, logLikelihood: -9042.002738
        [00:09:14] Bootstrap tree #110, logLikelihood: -9166.047993
        [00:09:19] Bootstrap tree #111, logLikelihood: -8984.506331
        [00:09:23] Bootstrap tree #112, logLikelihood: -9494.950162
        [00:09:27] Bootstrap tree #113, logLikelihood: -9416.419189
        [00:09:31] Bootstrap tree #114, logLikelihood: -8854.810142
        [00:09:34] Bootstrap tree #115, logLikelihood: -8861.231185
        [00:09:37] Bootstrap tree #116, logLikelihood: -9036.229486
        [00:09:41] Bootstrap tree #117, logLikelihood: -9596.513610
        [00:09:45] Bootstrap tree #118, logLikelihood: -9332.538535
        [00:09:50] Bootstrap tree #119, logLikelihood: -8932.400508
        [00:09:54] Bootstrap tree #120, logLikelihood: -9247.201142
        [00:09:58] Bootstrap tree #121, logLikelihood: -9229.817986
        [00:10:02] Bootstrap tree #122, logLikelihood: -9444.416973
        [00:10:07] Bootstrap tree #123, logLikelihood: -9354.443518
        [00:10:12] Bootstrap tree #124, logLikelihood: -9257.769601
        [00:10:16] Bootstrap tree #125, logLikelihood: -9389.733755
        [00:10:21] Bootstrap tree #126, logLikelihood: -8474.165924
        [00:10:25] Bootstrap tree #127, logLikelihood: -9171.745993
        [00:10:29] Bootstrap tree #128, logLikelihood: -9788.557990
        [00:10:32] Bootstrap tree #129, logLikelihood: -9119.466526
        [00:10:36] Bootstrap tree #130, logLikelihood: -8867.752762
        [00:10:40] Bootstrap tree #131, logLikelihood: -9179.894910
        [00:10:45] Bootstrap tree #132, logLikelihood: -8945.460082
        [00:10:48] Bootstrap tree #133, logLikelihood: -9003.649642
        [00:10:51] Bootstrap tree #134, logLikelihood: -9417.265448
        [00:10:56] Bootstrap tree #135, logLikelihood: -8865.613386
        [00:11:01] Bootstrap tree #136, logLikelihood: -9108.265377
        [00:11:05] Bootstrap tree #137, logLikelihood: -9426.548272
        [00:11:09] Bootstrap tree #138, logLikelihood: -9567.073566
        [00:11:12] Bootstrap tree #139, logLikelihood: -9478.292235
        [00:11:16] Bootstrap tree #140, logLikelihood: -8947.724446
        [00:11:20] Bootstrap tree #141, logLikelihood: -9442.889387
        [00:11:24] Bootstrap tree #142, logLikelihood: -9048.743435
        [00:11:27] Bootstrap tree #143, logLikelihood: -9148.837082
        [00:11:31] Bootstrap tree #144, logLikelihood: -9853.294822
        [00:11:37] Bootstrap tree #145, logLikelihood: -8897.903566
        [00:11:41] Bootstrap tree #146, logLikelihood: -9024.399639
        [00:11:45] Bootstrap tree #147, logLikelihood: -8935.223283
        [00:11:48] Bootstrap tree #148, logLikelihood: -9038.610453
        [00:11:51] Bootstrap tree #149, logLikelihood: -9149.430789
        [00:11:55] Bootstrap tree #150, logLikelihood: -9509.122415
        [00:11:58] Bootstrap tree #151, logLikelihood: -9191.044389
        [00:12:02] Bootstrap tree #152, logLikelihood: -9159.000436
        [00:12:07] Bootstrap tree #153, logLikelihood: -8960.001057
        [00:12:10] Bootstrap tree #154, logLikelihood: -9283.493995
        [00:12:14] Bootstrap tree #155, logLikelihood: -9473.315833
        [00:12:18] Bootstrap tree #156, logLikelihood: -9386.945901
        [00:12:21] Bootstrap tree #157, logLikelihood: -9310.892541
        [00:12:24] Bootstrap tree #158, logLikelihood: -8735.452400
        [00:12:29] Bootstrap tree #159, logLikelihood: -9203.634133
        [00:12:34] Bootstrap tree #160, logLikelihood: -9209.878386
        [00:12:37] Bootstrap tree #161, logLikelihood: -8877.980985
        [00:12:43] Bootstrap tree #162, logLikelihood: -9516.488028
        [00:12:46] Bootstrap tree #163, logLikelihood: -9375.213506
        [00:12:49] Bootstrap tree #164, logLikelihood: -9253.773447
        [00:12:52] Bootstrap tree #165, logLikelihood: -9260.171392
        [00:12:56] Bootstrap tree #166, logLikelihood: -9288.247590
        [00:13:01] Bootstrap tree #167, logLikelihood: -8965.240960
        [00:13:05] Bootstrap tree #168, logLikelihood: -9315.403491
        [00:13:10] Bootstrap tree #169, logLikelihood: -9278.792736
        [00:13:13] Bootstrap tree #170, logLikelihood: -9322.218585
        [00:13:16] Bootstrap tree #171, logLikelihood: -9242.759299
        [00:13:21] Bootstrap tree #172, logLikelihood: -9603.742270
        [00:13:26] Bootstrap tree #173, logLikelihood: -9042.850294
        [00:13:30] Bootstrap tree #174, logLikelihood: -9535.130878
        [00:13:34] Bootstrap tree #175, logLikelihood: -9211.245289
        [00:13:38] Bootstrap tree #176, logLikelihood: -9649.264692
        [00:13:41] Bootstrap tree #177, logLikelihood: -8867.266278
        [00:13:46] Bootstrap tree #178, logLikelihood: -9306.145366
        [00:13:50] Bootstrap tree #179, logLikelihood: -8866.753735
        [00:13:53] Bootstrap tree #180, logLikelihood: -8949.963666
        [00:13:57] Bootstrap tree #181, logLikelihood: -9653.692636
        [00:14:01] Bootstrap tree #182, logLikelihood: -9288.878039
        [00:14:04] Bootstrap tree #183, logLikelihood: -8976.951582
        [00:14:07] Bootstrap tree #184, logLikelihood: -9163.586633
        [00:14:11] Bootstrap tree #185, logLikelihood: -8924.121074
        [00:14:16] Bootstrap tree #186, logLikelihood: -9226.056937
        [00:14:20] Bootstrap tree #187, logLikelihood: -9227.675474
        [00:14:24] Bootstrap tree #188, logLikelihood: -9506.927599
        [00:14:29] Bootstrap tree #189, logLikelihood: -9324.899224
        [00:14:33] Bootstrap tree #190, logLikelihood: -9526.595737
        [00:14:38] Bootstrap tree #191, logLikelihood: -8887.072133
        [00:14:43] Bootstrap tree #192, logLikelihood: -9114.934940
        [00:14:48] Bootstrap tree #193, logLikelihood: -9176.296397
        [00:14:53] Bootstrap tree #194, logLikelihood: -8883.477728
        [00:14:57] Bootstrap tree #195, logLikelihood: -9201.830636
        [00:15:05] Bootstrap tree #196, logLikelihood: -8977.308831
        [00:15:09] Bootstrap tree #197, logLikelihood: -9075.410459
        [00:15:16] Bootstrap tree #198, logLikelihood: -8905.479835
        [00:15:25] Bootstrap tree #199, logLikelihood: -9007.466396
        [00:15:30] Bootstrap tree #200, logLikelihood: -9285.487873
        [00:15:34] Bootstrap tree #201, logLikelihood: -8953.217225
        [00:15:38] Bootstrap tree #202, logLikelihood: -8948.881810
        [00:15:42] Bootstrap tree #203, logLikelihood: -9178.734685
        [00:15:46] Bootstrap tree #204, logLikelihood: -9291.211299
        [00:15:49] Bootstrap tree #205, logLikelihood: -9384.988958
        [00:15:52] Bootstrap tree #206, logLikelihood: -9114.958357
        [00:15:56] Bootstrap tree #207, logLikelihood: -9570.685134
        [00:16:00] Bootstrap tree #208, logLikelihood: -8886.780682
        [00:16:03] Bootstrap tree #209, logLikelihood: -9246.250332
        [00:16:05] Bootstrap tree #210, logLikelihood: -8491.190837
        [00:16:10] Bootstrap tree #211, logLikelihood: -9261.252984
        [00:16:14] Bootstrap tree #212, logLikelihood: -9541.743007
        [00:16:17] Bootstrap tree #213, logLikelihood: -9115.790718
        [00:16:21] Bootstrap tree #214, logLikelihood: -9830.773157
        [00:16:25] Bootstrap tree #215, logLikelihood: -8894.494505
        [00:16:28] Bootstrap tree #216, logLikelihood: -8817.062285
        [00:16:31] Bootstrap tree #217, logLikelihood: -8911.299225
        [00:16:34] Bootstrap tree #218, logLikelihood: -9043.642264
        [00:16:37] Bootstrap tree #219, logLikelihood: -8701.114682
        [00:16:41] Bootstrap tree #220, logLikelihood: -9106.824938
        [00:16:44] Bootstrap tree #221, logLikelihood: -9397.339362
        [00:16:48] Bootstrap tree #222, logLikelihood: -9462.226720
        [00:16:52] Bootstrap tree #223, logLikelihood: -9696.714970
        [00:16:55] Bootstrap tree #224, logLikelihood: -9055.686087
        [00:16:58] Bootstrap tree #225, logLikelihood: -9151.052404
        [00:17:01] Bootstrap tree #226, logLikelihood: -8886.893301
        [00:17:04] Bootstrap tree #227, logLikelihood: -9282.487532
        [00:17:09] Bootstrap tree #228, logLikelihood: -9450.892959
        [00:17:13] Bootstrap tree #229, logLikelihood: -9881.246553
        [00:17:17] Bootstrap tree #230, logLikelihood: -9634.681224
        [00:17:20] Bootstrap tree #231, logLikelihood: -9167.133850
        [00:17:25] Bootstrap tree #232, logLikelihood: -9420.743509
        [00:17:28] Bootstrap tree #233, logLikelihood: -9408.768654
        [00:17:31] Bootstrap tree #234, logLikelihood: -8724.214559
        [00:17:34] Bootstrap tree #235, logLikelihood: -9097.440810
        [00:17:39] Bootstrap tree #236, logLikelihood: -9058.354148
        [00:17:43] Bootstrap tree #237, logLikelihood: -9147.511506
        [00:17:46] Bootstrap tree #238, logLikelihood: -8727.199770
        [00:17:49] Bootstrap tree #239, logLikelihood: -9079.244267
        [00:17:53] Bootstrap tree #240, logLikelihood: -9534.674532
        [00:17:57] Bootstrap tree #241, logLikelihood: -9421.791078
        [00:18:01] Bootstrap tree #242, logLikelihood: -9408.263976
        [00:18:05] Bootstrap tree #243, logLikelihood: -9089.353402
        [00:18:09] Bootstrap tree #244, logLikelihood: -8886.997377
        [00:18:13] Bootstrap tree #245, logLikelihood: -9202.259216
        [00:18:17] Bootstrap tree #246, logLikelihood: -8995.171966
        [00:18:21] Bootstrap tree #247, logLikelihood: -8723.420000
        [00:18:24] Bootstrap tree #248, logLikelihood: -9465.732142
        [00:18:26] Bootstrap tree #249, logLikelihood: -9076.303493
        [00:18:30] Bootstrap tree #250, logLikelihood: -9103.136225
        [00:18:33] Bootstrap tree #251, logLikelihood: -9101.822227
        [00:18:37] Bootstrap tree #252, logLikelihood: -9717.208720
        [00:18:40] Bootstrap tree #253, logLikelihood: -9529.145822
        [00:18:43] Bootstrap tree #254, logLikelihood: -9927.829788
        [00:18:46] Bootstrap tree #255, logLikelihood: -9149.096211
        [00:18:51] Bootstrap tree #256, logLikelihood: -9097.360369
        [00:18:53] Bootstrap tree #257, logLikelihood: -9432.735134
        [00:18:57] Bootstrap tree #258, logLikelihood: -9361.287009
        [00:19:01] Bootstrap tree #259, logLikelihood: -8909.534210
        [00:19:04] Bootstrap tree #260, logLikelihood: -9223.819996
        [00:19:08] Bootstrap tree #261, logLikelihood: -9323.321863
        [00:19:11] Bootstrap tree #262, logLikelihood: -9103.573140
        [00:19:15] Bootstrap tree #263, logLikelihood: -9434.228208
        [00:19:19] Bootstrap tree #264, logLikelihood: -9833.015962
        [00:19:23] Bootstrap tree #265, logLikelihood: -9607.305172
        [00:19:27] Bootstrap tree #266, logLikelihood: -9400.309514
        [00:19:31] Bootstrap tree #267, logLikelihood: -9047.668513
        [00:19:35] Bootstrap tree #268, logLikelihood: -9292.882400
        [00:19:39] Bootstrap tree #269, logLikelihood: -8886.920987
        [00:19:42] Bootstrap tree #270, logLikelihood: -8859.286561
        [00:19:46] Bootstrap tree #271, logLikelihood: -9180.712064
        [00:19:50] Bootstrap tree #272, logLikelihood: -9601.150734
        [00:19:53] Bootstrap tree #273, logLikelihood: -9021.279212
        [00:19:59] Bootstrap tree #274, logLikelihood: -9161.008750
        [00:20:03] Bootstrap tree #275, logLikelihood: -9559.730568
        [00:20:06] Bootstrap tree #276, logLikelihood: -9194.593598
        [00:20:11] Bootstrap tree #277, logLikelihood: -9552.532177
        [00:20:14] Bootstrap tree #278, logLikelihood: -8646.767989
        [00:20:18] Bootstrap tree #279, logLikelihood: -8765.922378
        [00:20:22] Bootstrap tree #280, logLikelihood: -8996.739430
        [00:20:26] Bootstrap tree #281, logLikelihood: -9393.223479
        [00:20:30] Bootstrap tree #282, logLikelihood: -8804.014042
        [00:20:33] Bootstrap tree #283, logLikelihood: -8774.069782
        [00:20:39] Bootstrap tree #284, logLikelihood: -9110.996286
        [00:20:42] Bootstrap tree #285, logLikelihood: -9445.080487
        [00:20:46] Bootstrap tree #286, logLikelihood: -9091.304693
        [00:20:50] Bootstrap tree #287, logLikelihood: -9040.827484
        [00:20:54] Bootstrap tree #288, logLikelihood: -8957.672573
        [00:20:57] Bootstrap tree #289, logLikelihood: -9387.990881
        [00:21:00] Bootstrap tree #290, logLikelihood: -9570.139874
        [00:21:04] Bootstrap tree #291, logLikelihood: -9090.021223
        [00:21:07] Bootstrap tree #292, logLikelihood: -8843.737805
        [00:21:10] Bootstrap tree #293, logLikelihood: -9248.071606
        [00:21:14] Bootstrap tree #294, logLikelihood: -9195.659748
        [00:21:18] Bootstrap tree #295, logLikelihood: -9382.426125
        [00:21:22] Bootstrap tree #296, logLikelihood: -9386.150029
        [00:21:25] Bootstrap tree #297, logLikelihood: -8958.476590
        [00:21:29] Bootstrap tree #298, logLikelihood: -9373.145034
        [00:21:32] Bootstrap tree #299, logLikelihood: -8869.083878
        [00:21:36] Bootstrap tree #300, logLikelihood: -8925.255459
        [00:21:39] Bootstrap tree #301, logLikelihood: -8873.074521
        [00:21:43] Bootstrap tree #302, logLikelihood: -9309.751118
        [00:21:47] Bootstrap tree #303, logLikelihood: -9304.862306
        [00:21:51] Bootstrap tree #304, logLikelihood: -9510.368297
        [00:21:54] Bootstrap tree #305, logLikelihood: -9161.402868
        [00:21:58] Bootstrap tree #306, logLikelihood: -9562.830560
        [00:22:02] Bootstrap tree #307, logLikelihood: -9612.513965
        [00:22:06] Bootstrap tree #308, logLikelihood: -8914.342850
        [00:22:09] Bootstrap tree #309, logLikelihood: -8857.691236
        [00:22:15] Bootstrap tree #310, logLikelihood: -9338.337190
        [00:22:18] Bootstrap tree #311, logLikelihood: -9491.826674
        [00:22:22] Bootstrap tree #312, logLikelihood: -9196.220981
        [00:22:25] Bootstrap tree #313, logLikelihood: -9053.840973
        [00:22:29] Bootstrap tree #314, logLikelihood: -8989.130859
        [00:22:33] Bootstrap tree #315, logLikelihood: -9435.841633
        [00:22:36] Bootstrap tree #316, logLikelihood: -8973.338549
        [00:22:39] Bootstrap tree #317, logLikelihood: -9327.651321
        [00:22:43] Bootstrap tree #318, logLikelihood: -9142.342728
        [00:22:47] Bootstrap tree #319, logLikelihood: -9514.571049
        [00:22:50] Bootstrap tree #320, logLikelihood: -8954.188196
        [00:22:54] Bootstrap tree #321, logLikelihood: -9128.444144
        [00:22:57] Bootstrap tree #322, logLikelihood: -9421.968470
        [00:23:01] Bootstrap tree #323, logLikelihood: -9700.405815
        [00:23:04] Bootstrap tree #324, logLikelihood: -9189.006257
        [00:23:08] Bootstrap tree #325, logLikelihood: -9275.934783
        [00:23:12] Bootstrap tree #326, logLikelihood: -8756.608926
        [00:23:16] Bootstrap tree #327, logLikelihood: -9039.083670
        [00:23:20] Bootstrap tree #328, logLikelihood: -9320.660642
        [00:23:25] Bootstrap tree #329, logLikelihood: -9013.478006
        [00:23:28] Bootstrap tree #330, logLikelihood: -9708.570074
        [00:23:31] Bootstrap tree #331, logLikelihood: -8953.393877
        [00:23:35] Bootstrap tree #332, logLikelihood: -9135.533171
        [00:23:39] Bootstrap tree #333, logLikelihood: -9333.648605
        [00:23:42] Bootstrap tree #334, logLikelihood: -9380.085769
        [00:23:45] Bootstrap tree #335, logLikelihood: -8827.389947
        [00:23:49] Bootstrap tree #336, logLikelihood: -8681.950240
        [00:23:52] Bootstrap tree #337, logLikelihood: -9341.951890
        [00:23:56] Bootstrap tree #338, logLikelihood: -9139.827370
        [00:23:59] Bootstrap tree #339, logLikelihood: -9297.609595
        [00:24:05] Bootstrap tree #340, logLikelihood: -9776.675242
        [00:24:09] Bootstrap tree #341, logLikelihood: -9028.740671
        [00:24:13] Bootstrap tree #342, logLikelihood: -9141.101966
        [00:24:17] Bootstrap tree #343, logLikelihood: -9036.347416
        [00:24:20] Bootstrap tree #344, logLikelihood: -8912.571250
        [00:24:24] Bootstrap tree #345, logLikelihood: -9488.886917
        [00:24:27] Bootstrap tree #346, logLikelihood: -8586.872621
        [00:24:33] Bootstrap tree #347, logLikelihood: -9525.741640
        [00:24:37] Bootstrap tree #348, logLikelihood: -9484.982088
        [00:24:40] Bootstrap tree #349, logLikelihood: -9668.489591
        [00:24:46] Bootstrap tree #350, logLikelihood: -8899.012759
        [00:24:49] Bootstrap tree #351, logLikelihood: -9079.240783
        [00:24:53] Bootstrap tree #352, logLikelihood: -9713.020256
        [00:24:58] Bootstrap tree #353, logLikelihood: -9408.037286
        [00:25:02] Bootstrap tree #354, logLikelihood: -9341.804411
        [00:25:07] Bootstrap tree #355, logLikelihood: -9235.745199
        [00:25:11] Bootstrap tree #356, logLikelihood: -9465.114899
        [00:25:15] Bootstrap tree #357, logLikelihood: -9629.981268
        [00:25:19] Bootstrap tree #358, logLikelihood: -9078.663961
        [00:25:23] Bootstrap tree #359, logLikelihood: -9099.814717
        [00:25:26] Bootstrap tree #360, logLikelihood: -8678.404562
        [00:25:29] Bootstrap tree #361, logLikelihood: -8825.208984
        [00:25:32] Bootstrap tree #362, logLikelihood: -9061.158321
        [00:25:35] Bootstrap tree #363, logLikelihood: -9228.407122
        [00:25:38] Bootstrap tree #364, logLikelihood: -8925.542807
        [00:25:42] Bootstrap tree #365, logLikelihood: -9801.069927
        [00:25:45] Bootstrap tree #366, logLikelihood: -9078.384018
        [00:25:48] Bootstrap tree #367, logLikelihood: -9278.901871
        [00:25:51] Bootstrap tree #368, logLikelihood: -9111.947870
        [00:25:55] Bootstrap tree #369, logLikelihood: -9056.505550
        [00:26:00] Bootstrap tree #370, logLikelihood: -9166.102049
        [00:26:03] Bootstrap tree #371, logLikelihood: -8859.990848
        [00:26:07] Bootstrap tree #372, logLikelihood: -9261.248759
        [00:26:11] Bootstrap tree #373, logLikelihood: -9661.776545
        [00:26:15] Bootstrap tree #374, logLikelihood: -8967.216687
        [00:26:18] Bootstrap tree #375, logLikelihood: -9612.684945
        [00:26:22] Bootstrap tree #376, logLikelihood: -9130.368293
        [00:26:26] Bootstrap tree #377, logLikelihood: -8954.338929
        [00:26:29] Bootstrap tree #378, logLikelihood: -8992.260985
        [00:26:33] Bootstrap tree #379, logLikelihood: -9253.704195
        [00:26:37] Bootstrap tree #380, logLikelihood: -9446.826635
        [00:26:41] Bootstrap tree #381, logLikelihood: -9354.525996
        [00:26:49] Bootstrap tree #382, logLikelihood: -9285.830637
        [00:26:52] Bootstrap tree #383, logLikelihood: -9616.838248
        [00:26:56] Bootstrap tree #384, logLikelihood: -9612.774740
        [00:26:59] Bootstrap tree #385, logLikelihood: -9414.013388
        [00:27:02] Bootstrap tree #386, logLikelihood: -9383.039738
        [00:27:06] Bootstrap tree #387, logLikelihood: -9195.446527
        [00:27:10] Bootstrap tree #388, logLikelihood: -9598.672870
        [00:27:13] Bootstrap tree #389, logLikelihood: -9699.949326
        [00:27:17] Bootstrap tree #390, logLikelihood: -9194.059019
        [00:27:20] Bootstrap tree #391, logLikelihood: -8756.168307
        [00:27:23] Bootstrap tree #392, logLikelihood: -8958.712745
        [00:27:27] Bootstrap tree #393, logLikelihood: -9300.235394
        [00:27:30] Bootstrap tree #394, logLikelihood: -9409.293137
        [00:27:34] Bootstrap tree #395, logLikelihood: -9006.847922
        [00:27:38] Bootstrap tree #396, logLikelihood: -9558.142948
        [00:27:42] Bootstrap tree #397, logLikelihood: -9417.633467
        [00:27:46] Bootstrap tree #398, logLikelihood: -9249.582368
        [00:27:50] Bootstrap tree #399, logLikelihood: -9433.731927
        [00:27:54] Bootstrap tree #400, logLikelihood: -9445.727372
        [00:27:58] Bootstrap tree #401, logLikelihood: -9269.653469
        [00:28:01] Bootstrap tree #402, logLikelihood: -9777.706135
        [00:28:04] Bootstrap tree #403, logLikelihood: -8880.443390
        [00:28:07] Bootstrap tree #404, logLikelihood: -9266.905029
        [00:28:10] Bootstrap tree #405, logLikelihood: -9343.049120
        [00:28:15] Bootstrap tree #406, logLikelihood: -9445.311111
        [00:28:19] Bootstrap tree #407, logLikelihood: -9249.489591
        [00:28:22] Bootstrap tree #408, logLikelihood: -9433.387359
        [00:28:25] Bootstrap tree #409, logLikelihood: -9342.078211
        [00:28:28] Bootstrap tree #410, logLikelihood: -9383.983700
        [00:28:32] Bootstrap tree #411, logLikelihood: -9570.889472
        [00:28:35] Bootstrap tree #412, logLikelihood: -9007.282879
        [00:28:40] Bootstrap tree #413, logLikelihood: -8880.978609
        [00:28:43] Bootstrap tree #414, logLikelihood: -9248.820824
        [00:28:47] Bootstrap tree #415, logLikelihood: -9311.414296
        [00:28:50] Bootstrap tree #416, logLikelihood: -9127.627223
        [00:28:53] Bootstrap tree #417, logLikelihood: -9144.599472
        [00:28:57] Bootstrap tree #418, logLikelihood: -9080.237452
        [00:29:00] Bootstrap tree #419, logLikelihood: -9402.309729
        [00:29:03] Bootstrap tree #420, logLikelihood: -9084.051143
        [00:29:07] Bootstrap tree #421, logLikelihood: -9342.457565
        [00:29:12] Bootstrap tree #422, logLikelihood: -8825.254005
        [00:29:16] Bootstrap tree #423, logLikelihood: -8908.204921
        [00:29:19] Bootstrap tree #424, logLikelihood: -9444.893132
        [00:29:23] Bootstrap tree #425, logLikelihood: -9235.095611
        [00:29:28] Bootstrap tree #426, logLikelihood: -9548.221339
        [00:29:32] Bootstrap tree #427, logLikelihood: -9138.065332
        [00:29:36] Bootstrap tree #428, logLikelihood: -8928.783040
        [00:29:39] Bootstrap tree #429, logLikelihood: -9683.373235
        [00:29:42] Bootstrap tree #430, logLikelihood: -8898.447438
        [00:29:47] Bootstrap tree #431, logLikelihood: -8930.431853
        [00:29:51] Bootstrap tree #432, logLikelihood: -9010.452630
        [00:29:54] Bootstrap tree #433, logLikelihood: -9240.662436
        [00:29:58] Bootstrap tree #434, logLikelihood: -9392.047576
        [00:30:02] Bootstrap tree #435, logLikelihood: -9515.733460
        [00:30:07] Bootstrap tree #436, logLikelihood: -9213.476288
        [00:30:11] Bootstrap tree #437, logLikelihood: -9370.652153
        [00:30:15] Bootstrap tree #438, logLikelihood: -8923.300904
        [00:30:19] Bootstrap tree #439, logLikelihood: -9377.231191
        [00:30:23] Bootstrap tree #440, logLikelihood: -9630.120504
        [00:30:27] Bootstrap tree #441, logLikelihood: -9072.011694
        [00:30:30] Bootstrap tree #442, logLikelihood: -9674.241806
        [00:30:34] Bootstrap tree #443, logLikelihood: -8962.598537
        [00:30:37] Bootstrap tree #444, logLikelihood: -8992.970543
        [00:30:43] Bootstrap tree #445, logLikelihood: -9138.516854
        [00:30:46] Bootstrap tree #446, logLikelihood: -9844.843812
        [00:30:49] Bootstrap tree #447, logLikelihood: -9345.323698
        [00:30:53] Bootstrap tree #448, logLikelihood: -9158.078130
        [00:30:57] Bootstrap tree #449, logLikelihood: -9683.278733
        [00:31:02] Bootstrap tree #450, logLikelihood: -9114.642169
        [00:31:02] Bootstrapping converged after 450 replicates.

        Optimized model parameters:

        Partition 0: noname
        Rate heterogeneity: GAMMA (4 cats, mean),  alpha: 2.259584 (ML),  weights&rates: (0.250000,0.322410) (0.250000,0.680054) (0.250000,1.075309) (0.250000,1.922227) 
        Base frequencies (model): 0.076748 0.051691 0.042645 0.051544 0.019803 0.040752 0.061830 0.073152 0.022944 0.053761 0.091904 0.058676 0.023826 0.040126 0.050901 0.068765 0.058565 0.014261 0.032102 0.066005 
        Substitution rates (model): 58.000000 54.000000 81.000000 56.000000 57.000000 105.000000 179.000000 27.000000 36.000000 30.000000 35.000000 54.000000 15.000000 194.000000 378.000000 475.000000 9.000000 11.000000 298.000000 45.000000 16.000000 113.000000 310.000000 29.000000 137.000000 328.000000 22.000000 38.000000 646.000000 44.000000 5.000000 74.000000 101.000000 64.000000 126.000000 20.000000 17.000000 528.000000 34.000000 86.000000 58.000000 81.000000 391.000000 47.000000 12.000000 263.000000 30.000000 10.000000 15.000000 503.000000 232.000000 8.000000 70.000000 16.000000 10.000000 49.000000 767.000000 130.000000 112.000000 11.000000 7.000000 26.000000 15.000000 4.000000 15.000000 59.000000 38.000000 4.000000 46.000000 31.000000 9.000000 5.000000 59.000000 69.000000 17.000000 23.000000 7.000000 31.000000 78.000000 14.000000 223.000000 42.000000 115.000000 209.000000 62.000000 323.000000 26.000000 597.000000 9.000000 72.000000 292.000000 43.000000 4.000000 164.000000 53.000000 51.000000 18.000000 24.000000 20.000000 119.000000 26.000000 12.000000 9.000000 181.000000 18.000000 5.000000 18.000000 30.000000 32.000000 10.000000 7.000000 45.000000 23.000000 6.000000 6.000000 27.000000 14.000000 5.000000 24.000000 201.000000 33.000000 55.000000 8.000000 47.000000 16.000000 56.000000 45.000000 33.000000 40.000000 115.000000 73.000000 46.000000 8.000000 573.000000 11.000000 229.000000 21.000000 479.000000 89.000000 10.000000 40.000000 245.000000 9.000000 32.000000 961.000000 14.000000 388.000000 248.000000 102.000000 59.000000 25.000000 52.000000 24.000000 180.000000 65.000000 4.000000 21.000000 47.000000 103.000000 10.000000 8.000000 14.000000 43.000000 16.000000 29.000000 226.000000 24.000000 18.000000 323.000000 17.000000 92.000000 12.000000 53.000000 536.000000 62.000000 285.000000 118.000000 6.000000 10.000000 23.000000 477.000000 35.000000 63.000000 38.000000 12.000000 21.000000 112.000000 71.000000 25.000000 16.000000 


        Final LogLikelihood: -9254.466015

        AIC score: 18560.932029 / AICc score: 18562.847445 / BIC score: 18681.398308
        Free parameters (model + branch lengths): 26

        Best ML tree saved to: /Users/melettedevore/Desktop/PP563/project/M4.raxml.bestTree
        All ML trees saved to: /Users/melettedevore/Desktop/PP563/project/M4.raxml.mlTrees
        Best ML tree with Felsenstein bootstrap (FBP) support values saved to: /Users/melettedevore/Desktop/PP563/project/M4.raxml.support
        Optimized model saved to: /Users/melettedevore/Desktop/PP563/project/M4.raxml.bestModel
        Bootstrap trees saved to: /Users/melettedevore/Desktop/PP563/project/M4.raxml.bootstraps

        Execution log saved to: /Users/melettedevore/Desktop/PP563/project/M4.raxml.log

        Analysis started: 04-May-2025 13:53:06 / finished: 04-May-2025 14:24:08

        Elapsed time: 1862.532 seconds

Visualizing with FigTree

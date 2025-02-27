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
        
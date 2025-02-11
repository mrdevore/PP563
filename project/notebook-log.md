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
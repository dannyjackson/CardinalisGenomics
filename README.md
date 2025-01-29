# Project Overview
This repository contains scripts and code used to analyze medium-coverage whole-genome sequence data from two Cardinalis species, each sampled in both urban and rural contexts.

### Central Question
How does urbanization affect two closely related species that differentially respond to urbanization?

Northern cardinals (*Cardinalis cardinalis*) appear slightly more tolerant of urbanization than pyrrhuloxia (*Cardinalis sinuatus*), despite existing in very similar ecological niches in the absence of hunam disturbance. We modeled this in Maxent as the foundational justification for this project in a paper that is still in revision (see [Spatial_Github](https://github.com/dannyjackson/Spatial_Github) for relevant code and latest draft). These analyses aim to identify signatures of selection that have acted on each species in an urban context. 

### Abbreviations:
   - **NOCA / noca**: northern cardinal (*Cardinalis cardinalis*)
   - **PYRR / pyrr**: pyrrhuloxia (*Cardinalis sinuatus*)

### Pipeline
To implement the pipeline, follow the markdown files in sequential order. Each markdown filename follows a structured format:
**[Stage Letter][Step Number]_[Description]**
   - The **letter** represents the pipeline stage.
   - The **number** indicates the order in which the scripts should be excecuted within that stage.
   - The **description** briefly summarizes the purpose of the script.
Each markdown file calls various scripts within this repository, all labeled using a similar syntax.

## Pipeline Stages

A. **Preprocessing**  
We follow the guidance for mitigating batch effects in low-coverage genomic sequence data from [Lou and Therkildsen 2022](https://doi.org/10.1111/1755-0998.13559).
   - Trims sequence reads
   - Aligns and filters the genome
   - Clips overlapping read pairs
   - Indel realignment
   - Removes individuals with low depth of coverage (<3x)
   - Generates a list of SNPs
   - Methods
      - [Trimmomatic](https://github.com/timflutre/trimmomatic)
      - [bwa](https://github.com/lh3/bwa)
      - [BamUtil clipOverlap](https://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap)
      - [gatk](https://gatk.broadinstitute.org/hc/en-us) (requires version 3.7-0 or earlier)
      - [samtools](https://www.htslib.org/)
      - [angsd](https://www.popgen.dk/angsd/)


B. **Population Structure Analysis**  
   - Identifies species and population clusters
   - Also identifies related individuals for filtering from downstream analysis
   - Methods
      - [PCAngsd](https://github.com/Rosemeis/pcangsd/)
         - *Principal component analysis (PCA)*
         - *Admixture modeling*

C. **Selection Analysis**  
   - Detects signatures of natural selection through:
     - **Direct comparisons of urban vs rural populations** 
       - Implemented in [angsd](https://www.popgen.dk/angsd/)
       - *F<sub>ST</sub>*
       - *dxy*
     - **Indirect comparisons of urban vs rural populations:**
       - Tests for selection signals independently in each population (urban and rural
       - Signals exclusive to the urban population are considered indicative of selection due to urbanization.
        - Methods
          - *Tajima's D* (implemented in [angsd](https://www.popgen.dk/angsd/))
          - *[RAiSD](https://github.com/alachins/raisd)*
  
D. **Gene Investigations**  
   - Diving deep into the context of interesting genes identified in Stage C.
     


### Acknowledgements
This project is the 4th chapter of my dissertation, which I defended in September 2023. It is a collaboration between myself, my PhD advisor **Kevin McGraw** (Michigan State University), **Scott Taylor** (University of Colorado), and my postdoctoral advisor **Sabrina McNew** (University of Arizona). I am very grateful for the support and feedback of the other members of my dissertation committee: **Gro Amdam** and **Karen Sweazea** (Arizona State University). 

All sequence data wille be available on the **Sequence Read Archive** upon publicaton.

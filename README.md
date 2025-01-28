# G2P-BIS Genomics To Proteins - Bio Integration Suite
[![DOI](https://img.shields.io/badge/DOI-10.1038/s41592--024--02409--0-blue.svg)](https://doi.org/10.1038/s41592-024-02409-0)
[![PubMed](https://img.shields.io/badge/PubMed-39294369-blue.svg)](https://pubmed.ncbi.nlm.nih.gov/39294369/)

A suite of tools and algorithms used by [Genomics 2 Proteins portal](https://g2p.broadinstitute.org) to integrate and align genomics and protein data.

## How to cite
If you use Genomics 2 Protein portal, please cite [the paper](https://www.nature.com/articles/s41592-024-02409-0):

- *Kwon, S., Safer, J., Nguyen, D.T. et al. Genomics 2 Proteins portal: a resource and discovery tool for linking genetic screening outputs to protein sequences and structures. Nat Methods 21, 1947–1957 (2024), [https://doi.org/10.1101/2024.01.02.573913](https://doi.org/10.1101/2024.01.02.573913)*

or, if you prefer the `BibTeX` format:

```
@article{kwon_genomics_2024,
	author = {Kwon, Seulki and Safer, Jordan and Nguyen, Duyen T. and Hoksza, David and May, Patrick and Arbesfeld, Jeremy A. and Rubin, Alan F. and Campbell, Arthur J. and Burgin, Alex and Iqbal, Sumaiya},
	title = {Genomics 2 Proteins portal: a resource and discovery tool for linking genetic screening outputs to protein sequences and structures},
	journal = {Nature Methods},
	volume = {21},
	pages = {1947--1957},
	year = {2024},
	month = oct,
	issn = {1548-7105},
	doi = {10.1038/s41592-024-02409-0},
	url = {https://doi.org/10.1038/s41592-024-02409-0},
}
```


## Protein-Protein Interactions

### Introduction

Proteins are complex molecules that rely on various types of interactions to maintain their structure and function. Genomics 2 Proteins portal provide intra- and inter-molecular protein interactions at each residue positions. Intra-protein interactions are calculated from both AlphaFold and monomeric PDB structures, and inter-protein interactions are calculated from available multimeric PDB structures. 


1. **Hydrogen Bonds**
1. **Nonbonding Interactions**
1. **Disulfide Bonds**
1. **Salt Bridges**

These interactions play critical roles in stabilizing protein structures and facilitating biological functions. Hydrogen bonds, nonbonding interactions, and disulfide bonds are calculated from HBPLUSv.3.06 software (https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/). Salt-bridges are calculated from the included python script (src/PPI_salt_bridges.py)

### Interaction Types


#### 1. Hydrogen Bonds

**Hydrogen bonds** are interactions between a hydrogen atom attached to an electronegative atom (donor) and another electronegative atom (acceptor). These are calculated using the software HBPLUS, which identifies potential hydrogen bonds based on geometric criteria. HBPLUS considers an interaction to qualify as a hydrogen bond if the donor-acceptor distance is less than 3.9 Ångströms and the hydrogen-acceptor distance is less than 2.5 Ångströms.


#### 2. Nonbonding Interactions

**Nonbonding interactions** include various types of interactions that do not involve covalent bonds or hydrogen bonds, such as van der Waals forces. HBPLUS is also used to calculate these interactions using its default settings where the distance between two residues is less than 3.9 Ångströms.

#### 3. Disulfide Bonds

**Disulfide bonds** are covalent bonds formed between the sulfur atoms of two cysteine residues. These bonds are crucial for stabilizing the three-dimensional structure of proteins. HBPLUS identifies these bonds based on the proximity (3.0 Ångströms) of sulfur atoms within the cysteine residues.

#### 4. Salt Bridges

A **salt bridge** is a type of non-covalent interaction between oppositely charged residues. In the provided script, salt bridges are calculated by identifying pairs of charged residues (Aspartic Acid (ASP), Glutamic Acid (GLU), Arginine (ARG), and Lysine (LYS)) and measuring the distance between their charged atoms. A salt bridge is considered present if the distance between a pair of oppositely charged atoms is less than 3.2 Ångströms.

## G2P APIs
Genomics To Proteins portal makes available two APIs, the G2P3D API for high throughput and scalable access to Gene-Protein Isoform-Transcript-Structure mapping and the Protein Feature API for high throughput and scalable access to protein feature data.

Swagger page: https://g2p.broadinsitute.org/api-docs/

Client Library: https://github.com/broadinstitute/g2papi

## Pairwise Protein Isoform Sequence Alignment
For pairwise protein isoform sequence alignment, Genomics To Proteins portal uses Gotoh's optimization for the Needleman-Wunsch algorithm. The implementation uses a gap open penalty of -12 and gap extend penalty of -1, as described in biopython based on Durbin et al's description. For isoform alignment, the implementation uses a mismatch penalty of -20 and a match reward of +1, as opposed to the BLASTP matrix commonly used for other peptides. The reason for this is that in the case of an isoform alignment for genetic interpretation, the purpose of the alignment is to identify identical sequences resulting from translation of the same exon shared between isoforms. 

Gotoh O. An improved algorithm for matching biological sequences. J Mol Biol. 1982 Dec 15;162(3):705-8. doi: 10.1016/0022-2836(82)90398-9. PMID: 7166760.

Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol. 1970 Mar;48(3):443-53. doi: 10.1016/0022-2836(70)90057-4. PMID: 5420325.

Durbin R, Eddy SR, Krogh A, Mitchison G. Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids. Cambridge University Press; 1998.

Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878

## Pocket annotations
For pocket detection and assigment, pockets are detected on AlphaFold structures using [fpocket](https://github.com/Discngine/fpocket). Pockets are ranked according the the computed fpocket druggability score, and if a residue participates in multiple pockets, it is annotated with the pocket with the highest druggability score. fpocket computes the list of residues which form the boundaries of a pocket as part of its computation. Metadata computed by fpocket is shown in the portal for each pocket. The metadata shown is: druggability score, pocket volume in cubic angstroms and the mean pLDDT of residues making up the pocket.

Le Guilloux, Vincent, Peter Schmidtke, and Pierre Tuffery. "Fpocket: an open source platform for ligand pocket detection." BMC bioinformatics 10 (2009): 1-11.

## Usage

```javascript
// Import the Needleman-Wunsch function from alignment.js
const { needlemanWunsch } = require('./core/alignment.js');

// Define your protein sequences
let sequence1 = 'MALSASPCANGAGGAEGGAAAAGGAUGA';
let sequence2 = 'MLAVSAAAGAAGGAAAAGGAUGA';

// Define your scoring parameters
let gapOpen = -12;
let gapExtend = -1;

// Call the Needleman-Wunsch function with the specified arguments
let result = needlemanWunsch(sequence1, sequence2, { gapOpen, gapExtend });

console.log(result);
```


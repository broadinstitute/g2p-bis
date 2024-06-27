# G2P-BIS Genomics To Proteins - Bio Integration Suite
A suite of tools and algorithms used by Genomics 2 Proteins portal to integrate and align genomics and protein data.

## Pairwise Isoform Sequence Alignment
For pairwise isoform sequence alignment, Genomics To Proteins portal uses Gotoh's optimization for the Needleman-Wunsch algorithm. The implementation uses a gap open penalty of -12 and gap extend penalty of -1, as described in biopython based on Durbin et al's description. For isoform alignment, the implementation uses a mismatch penalty of -20 and a match reward of +1, as opposed to the BLASTP matrix commonly used for other peptides. The reason for this is that in the case of an isoform alignment for genetic interpretation, the purpose of the alignment is to identify identical sequences resulting from translation of the same exon shared between isoforms. 

Gotoh O. An improved algorithm for matching biological sequences. J Mol Biol. 1982 Dec 15;162(3):705-8. doi: 10.1016/0022-2836(82)90398-9. PMID: 7166760.

Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol. 1970 Mar;48(3):443-53. doi: 10.1016/0022-2836(70)90057-4. PMID: 5420325.

Durbin R, Eddy SR, Krogh A, Mitchison G. Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids. Cambridge University Press; 1998.

Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878

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


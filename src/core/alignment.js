import BLASTP from "../constants/blastp.js";

function isoform_match_score_fn(residue1, residue2) {
  if (residue1 === residue2) {
    return 1;
  }
  return -20;
}

/**
 * Needleman-Wunsch sequence alignment with affine gap penalties.
 * Using the affine gap method described in Gotoh's 1982 paper.
 * https://doi.org/10.1016/0022-2836(82)90398-9
 *
 * Gap penalties and BLASTP matrix derived from Biopython's implementation.
 * https://github.com/biopython/biopython/blob/master/Bio/Align/__init__.py#L3894
 * Which is based on Durbin et al. 1998.
 *
 * Gap penalties are based on biopython's implementation
 * BLASTP matrix is biopython derivation of BLOSUM62
 **/

function initializeMatricesAffine(n, m, gapOpen, gapExtend) {
  const M = Array.from(Array(n + 1), () => Array(m + 1).fill(0));
  const Ix = Array.from(Array(n + 1), () => Array(m + 1).fill(0));
  const Iy = Array.from(Array(n + 1), () => Array(m + 1).fill(0));

  // Initialize the first row and column of Ix, Iy, and M
  for (let i = 1; i <= n; i++) {
    // Start from 1 to leave (0,0) as 0
    M[i][0] = gapOpen + (i - 1) * gapExtend; // Adjusted to start from gapOpen for the first gap
    Ix[i][0] = M[i][0]; // Ix should match M here since it's a gap opening in seq1
    Iy[i][0] = -Infinity; // Iy should be -inf for the first column
  }
  for (let j = 1; j <= m; j++) {
    // Start from 1 to leave (0,0) as 0
    M[0][j] = gapOpen + (j - 1) * gapExtend; // Adjusted to start from gapOpen for the first gap
    Iy[0][j] = M[0][j]; // Iy should match M here since it's a gap opening in seq2
    Ix[0][j] = -Infinity; // Ix should be -inf for the first row
  }

  // Set the top-left corner to 0 explicitly to clarify no score has been accumulated
  M[0][0] = 0;
  Ix[0][0] = 0; // Ensuring that the starting point for Ix is also 0
  Iy[0][0] = 0; // Ensuring that the starting point for Iy is also 0

  return { M, Ix, Iy };
}

/**
 * Fills the scoring matrices for sequence alignment with affine gap penalties.
 * @param {string} seq1 The first sequence.
 * @param {string} seq2 The second sequence.
 * @param {Object} matrices The scoring matrices (M, Ix, Iy) and traceback matrix.
 * @param {number} gapOpen The penalty for opening a gap.
 * @param {number} gapExtend The penalty for extending a gap.
 */
function fillScoreMatricesAffine(seq1, seq2, { M, Ix, Iy }, gapOpen, gapExtend) {
  const n = seq1.length;
  const m = seq2.length;

  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      // Use this for general alignments
      // const scoreMatchMismatch = BLASTP[seq1[i - 1]][seq2[j - 1]];
      const scoreMatchMismatch = isoform_match_score_fn(seq1[i - 1], seq2[j - 1]);

      // Calculate scores for M, Ix, and Iy matrices
      Ix[i][j] = Math.max(M[i - 1][j] + gapOpen, Ix[i - 1][j] + gapExtend);
      Iy[i][j] = Math.max(M[i][j - 1] + gapOpen, Iy[i][j - 1] + gapExtend);
      M[i][j] = Math.max(M[i - 1][j - 1] + scoreMatchMismatch, Ix[i][j], Iy[i][j]);
    }
  }
}

/**
 * Performs traceback to reconstruct the optimal alignment from the traceback matrix.
 * @param {string} seq1 The first sequence.
 * @param {string} seq2 The second sequence.
 * @param {Array} tracebackMatrix The traceback matrix.
 * @returns {Array} An array containing the aligned sequences.
 */
function tracebackAffine(seq1, seq2, M, Ix, Iy, gapOpen) {
  let i = seq1.length;
  let j = seq2.length;
  const alignedSeq1 = [];
  const alignedSeq2 = [];

  let currentMatrix = "M";

  while (i > 0 || j > 0) {
    if (currentMatrix === "M") {
      const currentCell = M[i][j];
      const upCell = Ix[i][j];
      const leftCell = Iy[i][j];
      if (currentCell === upCell) {
        currentMatrix = "Ix";
      } else if (currentCell === leftCell) {
        currentMatrix = "Iy";
      } else {
        alignedSeq1.unshift(seq1[i - 1]);
        alignedSeq2.unshift(seq2[j - 1]);
        i--;
        j--;
      }
    } else if (currentMatrix === "Ix") {
      const currentCell = Ix[i][j];
      const openCell = M[i - 1][j] + gapOpen;
      if (currentCell === openCell) {
        currentMatrix = "M";
      }
      i--;
      alignedSeq1.unshift(seq1[i]);
      alignedSeq2.unshift("-");
    } else if (currentMatrix === "Iy") {
      const currentCell = Iy[i][j];
      const openCell = M[i][j - 1] + gapOpen;
      if (currentCell === openCell) {
        currentMatrix = "M";
      }
      j--;
      alignedSeq1.unshift("-");
      alignedSeq2.unshift(seq2[j]);
    }
  }

  return [alignedSeq1.join(""), alignedSeq2.join("")];
}

/**
 * The main function to perform Needleman-Wunsch sequence alignment with affine gap penalties.
 * @param {string} seq1 The first sequence.
 * @param {string} seq2 The second sequence.
 * @param {Object} options An object containing scoring options (gapOpen, gapExtend, defaults from).
 * @returns {Array} An array containing the aligned sequences and their alignment score.
 */
function needlemanWunsch(seq1, seq2, { gapOpen = -12, gapExtend = -1 } = {}) {
  const { M, Ix, Iy } = initializeMatricesAffine(seq1.length, seq2.length, gapOpen, gapExtend);

  fillScoreMatricesAffine(seq1, seq2, { M, Ix, Iy }, gapOpen, gapExtend);

  const [alignedSeq1, alignedSeq2] = tracebackAffine(seq1, seq2, M, Ix, Iy, gapOpen);
  const score = Math.max(
    M[seq1.length][seq2.length],
    Ix[seq1.length][seq2.length],
    Iy[seq1.length][seq2.length]
  );

  return [alignedSeq1, alignedSeq2, score];
}

const generateConservationLine = (seq1Chunk, seq2Chunk) => {
  let conservationLine = "";
  for (let i = 0; i < seq1Chunk.length; i++) {
    if (seq1Chunk[i] === seq2Chunk[i]) {
      conservationLine += "*"; // Exact match
    } else {
      // Check if the amino acids are similar, dissimilar, or no match
      if (seq1Chunk[i] === "-" || seq2Chunk[i] === "-") {
        conservationLine += " "; // No match
      } else {
        const score = BLASTP[seq1Chunk[i]][seq2Chunk[i]];
        if (score > 0) {
          conservationLine += ":"; // Similar
        } else {
          conservationLine += "."; // Dissimilar
        }
      }
    }
  }
  return conservationLine;
};

function toClustalWFormat(header1, seq1, header2, seq2) {
  const clustalWHeader = "CLUSTAL W alignment\n\n";
  // Break sequences into chunks for formatting, assuming 60 characters per line
  const chunkSize = 60;
  let formattedText = clustalWHeader;
  for (let i = 0; i < seq1.length; i += chunkSize) {
    const seq1Chunk = seq1.substring(i, i + chunkSize);
    const seq2Chunk = seq2.substring(i, i + chunkSize);
    const headerSpace = " ".repeat(header1.length);
    const margin = " ".repeat(6);
    formattedText += header1 + margin + seq1Chunk + "\n";
    formattedText += header2 + margin + seq2Chunk + "\n";
    formattedText += headerSpace + margin + generateConservationLine(seq1Chunk, seq2Chunk) + "\n\n";
  }
  return formattedText;
}

function clustalWDownloadURL(header1, seq1, header2, seq2) {
  const clustalWContent = toClustalWFormat(header1, seq1, header2, seq2);
  const blob = new Blob([clustalWContent], { type: "text/plain" });
  return {
    url: URL.createObjectURL(blob),
    filename: `G2P_${header1}_${header2}_alignment.clustalw.txt`,
  };
}

export { needlemanWunsch, clustalWDownloadURL };

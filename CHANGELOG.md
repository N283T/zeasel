# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- **Alphabet enhancements** — degeneracy bitmask tables, `avgScore`, `expectScore`, `count`, `match`, `dealign`, `dsqrlen`, `validateSeq`, `isResidue` (#34)
- **Sequence enhancements** — `SequenceBlock` for threaded pipelines, metadata setters (`setName`, `setAccession`, `setDescription`, `setSource`), `crcChecksum`, `countResidues`, `convertDegen2X`, `guessAlphabet`, `eql`, `isValid` (#35)
- **MSA enhancements** — `sequenceSubset`, `minimGaps`, `noGaps`, `reverseComplement`, `convertDegen2X`, `checkUniqueNames`, `hash`, `checksum`, `validate`, `compare`, `reasonableRF`, `guessAlphabet`, annotation management (`addGF`/`addGS`/`appendGC`/`appendGR`), per-sequence metadata (#36)
- **Score matrix enhancements** — BLOSUM45/80, PAM30/70/120/250, `getByName` lookup, `isSymmetric`, `maxScore`, `minScore`, `expectedScore`, `relativeEntropy` (#37)
- **SSI enhancements** — alias/secondary keys, `findNumber`, multi-file indexing, `findSubseq` (#38)
- **Tree enhancements** — `readNewick` parser, `singleLinkage`, `completeLinkage`, `wpgma` clustering (#39)
- **Matrix enhancements** — `luDecompose`, `invert`, `exp(tQ)`, `sum`, `min`, `max`, `frobeniusNorm` (#40)
- **Stats enhancements** — `stats/functions.zig` with `logGamma`, `psi`, `trigamma`, `erfc`, `chiSquaredTest`, `gTest`, `linearRegression`, `meanVariance` (#41)
- **Threading** — `threads.zig` with `ThreadPool` and `WorkQueue` (producer-consumer pattern) (#42)
- **dsqdata threading** — `PrefetchReader` with background loader thread for parallel I/O (#43)
- **Evolutionary models** — `ratematrix.zig` (Q matrix, normalization, P=exp(tQ)), `paml.zig` (PAML format reader), `mixdchlet.zig` (mixture Dirichlet priors) (#44)
- **Clustering** — `cluster.zig` with `singleLinkageFromDist` (#45)
- **I/O formats** — PHYLIP (`io/phylip.zig`), A2M (`io/a2m.zig`), PSI-BLAST (`io/psiblast.zig`), SELEX (`io/selex.zig`) (#46)
- **Data structures** — `varint.zig` (base-128 encoding), `graph.zig` (bipartite matching), `red_black.zig` (red-black tree), `huffman.zig` (Huffman coding) (#47)
- **Utilities** — `cpu.zig` (SIMD detection), `iset.zig` (independent set splitting), `recorder.zig` (I/O recorder with rewind) (#48)
- **Composition** — `composition.zig` with BL62, WAG, SW34, SW50 background frequencies (#37)
- SSI (Sequence/Subsequence Index) for fast random access to sequences (#14)
- Clustal and aligned FASTA format parsers (#10)
- Phylogenetic tree construction (UPGMA) and Newick I/O (#13)
- Matrix operations module (#12)
- RNA secondary structure WUSS notation parsing and validation (#15)
- dsqdata binary sequence database format (#16)
- Stockholm markup round-trip (GF/GC/GS/GR tags) (#17)
- Normal and Weibull statistical distributions (#8)
- MSA clustering (single/complete/average linkage) and manipulation (#11)
- Genetic code module with NCBI translation tables and codon translation (#9)
- Histogram with descriptive statistics and chi-squared goodness-of-fit (#7)
- Numerical optimization: conjugate gradient minimizer and Brent's root finder (#6)
- Sequence distance calculations (pairwise identity, Jukes-Cantor) (#5)
- Dirichlet and mixture Dirichlet priors (#4)
- MSA position-based and BLOSUM sequence weighting (#3)
- Random number generator and random sequence generation/shuffling (#2)
- Score matrix module with BLOSUM62 (#1)
- GenBank and EMBL format parsers and writers
- CLI tools: `zeasel-seqstat`, `zeasel-reformat`, `zeasel-seqfetch`
- Gumbel, exponential, and gamma statistical distributions
- SIMD vector operations with comptime generics and `@Vector`
- Stockholm format parser and writer
- MSA struct with `fromText`, `extractSeq`, `selectColumns`
- Unified `Reader`/`Writer` I/O interface with format auto-detection
- FASTA parser and writer
- `Sequence` struct with `init`, `fromText`, `clone`, `toText`, `reverseComplement`, `subseq`
- `Alphabet` module with DNA, RNA, and amino acid alphabets (comptime-generated)
- `digitize`, `textize`, `reverseComplement`, and alphabet type guessing

### Changed

- Renamed project from `zeacel` to `zeasel` (#32)

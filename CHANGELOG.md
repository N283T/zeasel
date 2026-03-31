# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

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

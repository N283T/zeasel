# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- **SSI dual-format support** — Read both zeasel native and Easel SSI v3.0 index formats, auto-detected by magic bytes. Easel SSI uses disk-based binary search for O(log n) lookup without loading keys into memory (#63)
- **Comprehensive Easel review** — 58 issues identified and fixed via systematic comparison with Easel C source (#61-#72, PRs #73-#101)

### Fixed (Easel review — CRITICAL)

- Fix Polak-Ribière beta denominator bug in conjugate gradient minimizer (#61)

### Fixed (Easel review — HIGH)

- Fix amino acid U/O degeneracy: U→C, O→K instead of all 20 residues (#64)
- Replace sequence checksum with Jenkins one-at-a-time hash for Easel compatibility (#64)
- Improve `guessType` with robust heuristics: 98% threshold, RNA/DNA distinction, minimum residues (#64)
- Fix `reverseComplement` to swap start/end coordinates and invalidate secondary structure (#64)
- Fix `countResidues` to distribute degenerate symbols across canonical residues (#64)
- Register phylip/a2m/psiblast/selex in Reader/Writer Format dispatch (#65)
- Rewrite A2M parser to correctly handle insert/consensus column distinction (#65)
- Add Stockholm multi-MSA parsing, PHYLIP sequential format, gzip/stdin support (#65)
- Add SSI random-access fetch/fetchSubseq integration in Reader (#65)
- Add `WindowReader` for streaming chromosome-scale sequence reading (#65)
- Add DDBJ/Pfam format variants, CLUSTAL-like detection, GenBank header detection (#65)
- Fix FASTA parser to reject empty sequences (#65)
- Fix PB weights: per-sequence-length normalization, degenerate residue handling (#66)
- Fix `percentIdentity` to use MIN(len1,len2) denominator matching Easel (#66)
- Add GSC (Gerstein-Sonnhammer-Chothia) weighting via UPGMA guide tree (#66)
- Add `idFilter` for identity-based sequence filtering (#66)
- Add bootstrap resampling, vertical shuffle, FlushLeftInserts, broken base-pair repair to msa_ops (#66)
- Replace erf approximation (1.5e-7) with C library erfc for full double precision (#67)
- Fix normal CDF/surv to use erfc avoiding catastrophic cancellation in tails (#67)
- Fix Gumbel surv() precision loss; add bisection fallback to fitComplete (#67)
- Add Gumbel fitCompleteLoc, fitCensored, fitCensoredLoc (#67)
- Add MLE fitting for exponential, gamma (Minka), and Weibull distributions (#67)
- Add log-space functions (logPdf/logCdf/logSurv) for all distributions (#67)
- Add gamma inverse CDF via bisection (#67)
- Replace minimizer backtracking with bracket+Brent line search (#67)
- Add histogram dynamic resizing, setTail, setExpectedTail (#67)
- Add Dirichlet sampling (sample, sampleUniform) (#67)
- Add CTG/TTG as standard code start codons (#68)
- Add all 18 NCBI translation tables with `byId()` lookup (#68)
- Add ambiguity-aware codon translation via degen bitmasks (#68)
- Add ORF finding pipeline with 6-frame search (#68)
- Add gencode `writeTable()` for human-readable output (#68)
- Fix `expectedScore` to take separate fi/fj for query/target backgrounds (#69)
- Widen ScoreMatrix from i8 to i16 with Kp×Kp (28×28) degenerate residue support (#69)
- Add BLAST format score matrix Read/Write (#69)
- Add `relativeEntropy` probability validation (#69)
- Fix PAML reader to permute from PAML alphabet order to zeasel order (#70)
- Generalize Jukes-Cantor to accept alphabet_size parameter (DNA K=4, protein K=20) (#70)
- Fix Newick parser to return error on multifurcation instead of silent data loss (#70)
- Add JukesCantor, Kimura, WAG rate matrix constructors (#70)
- Add distance variance estimation and average pairwise identity with subsampling (#70)
- Add tree compare (topological equivalence) and validate (internal consistency) (#70)
- Add random tree simulation (#70)
- Fix ratematrix `probMatrix` unused allocator parameter (#70)
- Add P-matrix validation, ratematrix relativeEntropy/expectedScore (#70)
- Add MixDirichlet Read/Write I/O, validate, compare (#70)
- Add DualWorkQueue with reader/worker buffer-recycling pattern (#71)
- Replace WorkQueue ArrayList with O(1) circular buffer (#71)
- Add ThreadPool StartBarrier, WorkerInfo, per-worker data, comprehensive tests (#71)
- Add Huffman encode/decode operations (#72)
- Add WUSS noPseudo (pseudoknot removal) and reverse (complement) (#72)
- Implement Altschul-Erickson dinucleotide-preserving shuffle (#72)
- Add Cobalt/Blue independent set algorithms (Petti & Eddy 2022) with LinkFn callback (#72)
- Add VectorOps: logNorm, relativeEntropy, sort, validate, log/exp element-wise (#72)
- Add varint exponential-Golomb and Golomb-Rice encodings (#72)
- Add red-black tree toSortedSlice (in-order traversal) (#72)
- Add random gamma deviate (Marsaglia-Tsang) and uniformPositive sampling (#72)
- Add recorder addLine, markBlock, getBlock for stream integration (#72)
- Add matrix Jacobi eigenvalue decomposition (#72)
- Fix Matrix.scale() to return new matrix (immutable convention) (#72)
- Add sequence appendDigital for incremental building (#72)
- Add runtime CPU feature detection via CPUID (x86_64) and NEON (aarch64) (#72)

### Previously added

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

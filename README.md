# zeasel

A Zig reimplementation of [Easel](https://github.com/EddyRivasLab/easel), the C library for biological sequence analysis from the Eddy/Rivas laboratory. Easel underpins [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/). zeasel serves as the foundation for a Zig reimplementation of HMMER.

## Features

- **Biological alphabets** — DNA, RNA, amino acid with comptime lookup tables, degeneracy resolution, scoring/counting, robust alphabet guessing
- **Sequence I/O** — FASTA, Stockholm (multi-MSA), GenBank, EMBL, DDBJ, Clustal, aligned FASTA, PHYLIP (interleaved+sequential), A2M, PSI-BLAST, SELEX, Pfam; gzip auto-decompression, stdin support, SSI random-access fetch
- **Streaming I/O** — `WindowReader` for chromosome-scale sequences with overlapping windows
- **Multiple sequence alignment** — MSA operations, weighting (PB/BLOSUM/GSC), clustering, gap manipulation, bootstrap resampling, vertical shuffle, FlushLeftInserts, identity filtering (IDFilter), broken base-pair repair
- **SIMD vector operations** — comptime-generic `@Vector`-based ops with automatic ISA selection; logNorm, relativeEntropy, sort, validate, log/exp
- **Statistical distributions** — Gumbel, gamma, exponential, normal, Weibull, Dirichlet, mixture Dirichlet; full PDF/CDF/surv/logPdf/logCdf/logSurv/invcdf; MLE fitting for all distributions; censored/truncated Gumbel fitting
- **Statistical functions** — logGamma, digamma, chi-squared test, G-test, linear regression, high-precision erfc (C library)
- **Numerical optimization** — conjugate gradient minimizer with bracket+Brent line search, Brent's root finder
- **Score matrices** — BLOSUM45/62/80, PAM30/70/120/250; i16 scores with Kp×Kp degenerate residue support; BLAST format I/O; named lookup, probability analysis
- **Evolutionary models** — rate matrix Q with JukesCantor/Kimura/WAG constructors, probability matrix P=exp(tQ), PAML reader with alphabet permutation, mixture Dirichlet I/O
- **Genetic code** — all 18 NCBI translation tables, ambiguity-aware codon translation, ORF finding pipeline (6-frame)
- **Sequence distance** — pairwise identity (MIN denominator), generalized Jukes-Cantor (any K), Kimura, variance estimation, average pairwise identity with subsampling
- **Phylogenetic trees** — UPGMA, single/complete linkage, WPGMA, Newick read/write (multifurcation detection), tree compare/validate, random tree simulation
- **SSI index** — sequence/subsequence index with aliases, multi-file, subseq positioning, Reader fetch integration
- **dsqdata** — binary packed sequence database with threaded prefetch reader, DualWorkQueue buffer recycling
- **Threading** — thread pool with start barrier and per-worker data, WorkQueue (O(1) circular buffer), DualWorkQueue (reader/worker pattern), runtime CPU detection (CPUID/NEON)
- **Data structures** — red-black tree (with sorted slice export), Huffman coding (encode/decode), varint (Google + exp-Golomb + Rice), bipartite matching, independent set splitting (Cobalt/Blue algorithms)
- **Matrix** — dense matrix with LU decomposition, inverse, exp(tQ), Jacobi eigenvalue decomposition
- **WUSS notation** — RNA secondary structure parsing, validation, pseudoknot removal, reverse complement

## Requirements

- [Zig](https://ziglang.org/) 0.15.0 or later

## Build

```bash
# Build the library
zig build

# Run all tests
zig build test

# Build CLI tools
zig build tools
```

## CLI Tools

| Tool | Description | Easel equivalent |
|------|-------------|------------------|
| `zeasel-seqstat` | Report sequence file statistics | `esl-seqstat` |
| `zeasel-reformat` | Convert between sequence formats | `esl-reformat` |
| `zeasel-seqfetch` | Fetch sequences by name | `esl-seqfetch` |

```bash
# After `zig build tools`, binaries are in zig-out/bin/
./zig-out/bin/zeasel-seqstat sequences.fasta
./zig-out/bin/zeasel-reformat stockholm input.sto
./zig-out/bin/zeasel-seqfetch sequences.fasta seq1 seq2
```

## Usage as a Library

Add zeasel as a Zig dependency, then import:

```zig
const zeasel = @import("zeasel");

// Digitize a DNA sequence
const abc = &zeasel.alphabet.dna;
const dsq = try abc.digitize(allocator, "ACGTACGT");
defer allocator.free(dsq);

// Read sequences from a FASTA file
var reader = try zeasel.io.Reader.fromMemory(allocator, abc, file_data, .fasta);
defer reader.deinit();
const sequences = try reader.readAll();
```

## Project Structure

```
src/
├── root.zig             # Public API re-exports
├── alphabet.zig         # Biological alphabets + degeneracy tables
├── sequence.zig         # Sequence type + SequenceBlock
├── io.zig               # I/O module hub
│   └── io/              # Format-specific parsers/writers
│       ├── reader.zig, writer.zig, window_reader.zig
│       ├── fasta.zig, stockholm.zig, genbank.zig
│       ├── clustal.zig, afa.zig, phylip.zig
│       ├── a2m.zig, psiblast.zig, selex.zig
├── msa.zig              # Multiple sequence alignment
├── msa_ops.zig          # MSA clustering and manipulation
├── msa_weight.zig       # MSA sequence weighting
├── score_matrix.zig     # Substitution matrices (BLOSUM, PAM)
├── composition.zig      # Background residue frequencies
├── distance.zig         # Sequence distance calculations
├── tree.zig             # Phylogenetic trees (UPGMA, Newick)
├── matrix.zig           # Matrix operations + LU/inverse/exp
├── ratematrix.zig       # Evolutionary rate matrices
├── paml.zig             # PAML model reader
├── mixdchlet.zig        # Mixture Dirichlet priors
├── ssi.zig              # Sequence/Subsequence Index
├── dsqdata.zig          # Binary sequence database + prefetch
├── gencode.zig          # Genetic code / translation
├── wuss.zig             # RNA secondary structure (WUSS)
├── simd.zig             # SIMD vector operations
│   └── simd/vector_ops.zig
├── stats.zig            # Statistical distributions hub
│   └── stats/           # Distribution modules + functions
├── threads.zig          # Thread pool + work queue
├── cpu.zig              # CPU/SIMD feature detection
├── cluster.zig          # Generalized clustering
├── red_black.zig        # Red-black tree
├── huffman.zig          # Huffman coding
├── varint.zig           # Variable-length integer encoding
├── graph.zig            # Bipartite matching
├── iset.zig             # Independent set splitting
├── recorder.zig         # I/O recorder with rewind
├── util.zig             # Utilities (RNG, random sequences)
└── tools/               # CLI tools
    ├── seqstat.zig
    ├── reformat.zig
    └── seqfetch.zig
```

## Easel Module Coverage

Full comparison of Easel's C modules and their zeasel status.

### Implemented in zeasel

| Easel module | zeasel module | Notes |
|---|---|---|
| `esl_alphabet` | `alphabet.zig` | Comptime alphabets + degeneracy tables, scoring/counting functions |
| `esl_sq` | `sequence.zig` | Digital sequences + SequenceBlock for threaded pipelines |
| `esl_msa` | `msa.zig` | Full MSA operations: subset, gaps, RC, metadata, validation |
| `esl_msaweight` | `msa_weight.zig` | PB weights (consensus columns, per-seq normalization), BLOSUM weights, GSC weights, IDFilter |
| `esl_msacluster` | `msa_ops.zig` | Single-linkage clustering by % identity |
| `esl_msashuffle` | `msa_ops.zig` | Column/vertical shuffling, bootstrap resampling, FlushLeftInserts, broken base-pair repair |
| `esl_scorematrix` | `score_matrix.zig` | BLOSUM45/62/80, PAM30/70/120/250; i16 Kp×Kp with degenerate support; BLAST format I/O |
| `esl_composition` | `composition.zig` | BL62, WAG, SW34, SW50 background frequencies |
| `esl_distance` | `distance.zig` | Pairwise identity (MIN denominator), generalized JC, Kimura, variance, average identity |
| `esl_tree` | `tree.zig` | UPGMA, single/complete/WPGMA, Newick I/O, compare, validate, simulate |
| `esl_dmatrix` | `matrix.zig` | Dense matrix + LU decomposition, inverse, exp(tQ), Jacobi eigenvalue decomposition |
| `esl_ssi` | `ssi.zig` | SSI index with aliases, multi-file, subseq positioning |
| `esl_dsqdata` | `dsqdata.zig` | Binary sequence DB + threaded PrefetchReader |
| `esl_gencode` | `gencode.zig` | All 18 NCBI tables, ambiguity-aware translation, ORF finding pipeline |
| `esl_wuss` | `wuss.zig` | WUSS notation: parse, validate, noPseudo, reverse |
| `esl_cluster` | `cluster.zig` | Generalized single-linkage clustering |
| `esl_vectorops` | `simd/vector_ops.zig` | Comptime-generic SIMD; logNorm, relativeEntropy, sort, validate |
| `esl_gumbel` | `stats/gumbel.zig` | Full PDF/CDF/surv/log-space; complete + censored MLE fitting |
| `esl_exponential` | `stats/exponential.zig` | Full PDF/CDF/surv/log-space; MLE fitting |
| `esl_gamma` | `stats/gamma.zig` | Full PDF/CDF/surv/log-space/invcdf; Minka MLE fitting |
| `esl_normal` | `stats/normal.zig` | High-precision erfc-based CDF/surv |
| `esl_weibull` | `stats/weibull.zig` | Full PDF/CDF/surv/log-space; CG minimizer MLE fitting |
| `esl_dirichlet` | `stats/dirichlet.zig` | Dirichlet distribution + sampling |
| `esl_histogram` | `stats/histogram.zig` | Dynamic resizing, tail fitting, goodness-of-fit |
| `esl_minimizer` | `stats/minimizer.zig` | Conjugate gradient with bracket+Brent line search |
| `esl_rootfinder` | `stats/rootfinder.zig` | Brent's method root finder |
| `esl_stats` | `stats/functions.zig` | logGamma, psi, trigamma, erfc, chi-squared, G-test, linear regression |
| `esl_mixdchlet` | `mixdchlet.zig` | Mixture Dirichlet priors; Read/Write I/O, validate, compare |
| `esl_ratematrix` | `ratematrix.zig` | Rate matrix Q with JC/Kimura/WAG constructors, P=exp(tQ), validation |
| `esl_paml` | `paml.zig` | PAML reader with PAML→zeasel alphabet permutation |
| `esl_red_black` | `red_black.zig` | Red-black tree (HMMER sparse DP) |
| `esl_huffman` | `huffman.zig` | Huffman coding with encode/decode operations |
| `esl_varint` | `varint.zig` | Google varint, exponential-Golomb, Golomb-Rice encoding |
| `esl_graph` | `graph.zig` | Maximum bipartite matching |
| `esl_iset` | `iset.zig` | Cobalt/Blue independent set splitting (Petti & Eddy 2022) |
| `esl_recorder` | `recorder.zig` | Line-based I/O with rewind/replay, stream integration, block marking |
| `esl_cpu` | `cpu.zig` | Comptime + runtime (CPUID/NEON) feature detection |
| `esl_threads` | `threads.zig` | ThreadPool with start barrier, WorkQueue (circular), DualWorkQueue |
| `esl_random` | `util/random.zig` | Mersenne Twister RNG + gamma deviate + uniformPositive |
| `esl_randomseq` | `util/random_seq.zig` | Random sequence generation, Altschul-Erickson dinucleotide shuffle |
| `esl_sqio_ascii` | `io/fasta.zig` | FASTA/EMBL sequence I/O |
| `esl_msafile_stockholm` | `io/stockholm.zig` | Stockholm format |
| `esl_msafile_clustal` | `io/clustal.zig` | Clustal format |
| `esl_msafile_afa` | `io/afa.zig` | Aligned FASTA |
| `esl_msafile_phylip` | `io/phylip.zig` | PHYLIP interleaved + sequential auto-detection |
| `esl_msafile_a2m` | `io/a2m.zig` | A2M (UCSC SAM) format |
| `esl_msafile_psiblast` | `io/psiblast.zig` | PSI-BLAST format |
| `esl_msafile_selex` | `io/selex.zig` | SELEX format |

### Not needed (replaced by Zig standard library or language features)

| Easel module | Zig replacement | Why not needed |
|---|---|---|
| `esl_alloc` | `std.mem.Allocator` | Zig allocator interface supports alignment natively |
| `esl_arr2`, `esl_arr3` | Native slices | Zig multi-dimensional arrays are first-class |
| `esl_avx`, `esl_sse`, `esl_neon`, `esl_vmx`, `esl_avx512` | `@Vector` built-in | Zig comptime SIMD replaces per-ISA C files |
| `esl_bitfield` | `std.StaticBitSet` | Standard library bit set |
| `esl_buffer` | `[]const u8` slices | Zig slice-based parsers don't need buffered I/O abstraction |
| `esl_fileparser` | `std.mem.tokenizeAny` | Standard library tokenizer |
| `esl_getopts` | `std.process.argsAlloc` | Zig arg parsing or third-party libraries |
| `esl_heap` | `std.PriorityQueue` | Standard library priority queue |
| `esl_json` | `std.json` | Standard library JSON parser |
| `esl_keyhash` | `std.StringHashMap` | Standard library hash map |
| `esl_mem` | Zig slice operations | `std.mem.eql`, `std.mem.indexOf`, etc. |
| `esl_quicksort` | `std.mem.sort` | Standard library sort |
| `esl_rand64` | `std.Random` | Standard library RNG |
| `esl_regexp` | `std.mem` string ops | Pattern matching via standard library |
| `esl_stack` | `std.ArrayList` | Standard library dynamic array as stack |
| `esl_stopwatch` | `std.time.Timer` | Standard library timer |
| `esl_subcmd` | Direct dispatch | Zig switch-based subcommand dispatch |
| `esl_workqueue` | `threads.WorkQueue` | Integrated into threads module |

### Not applicable

| Easel module | Reason |
|---|---|
| `esl_mpi` | MPI cluster parallelism — out of scope |
| `esl_msafile2` | Legacy Pfam Stockholm reader — superseded |
| `esl_sqio_ncbi` | NCBI BLAST database format — specialized |
| `esl_swat` | Smith-Waterman — marked UNFINISHED in Easel |
| `esl_hmm` | General discrete HMMs — not profile HMMs (HMMER has its own) |
| `esl_gev`, `esl_mixgev` | Generalized extreme value — research only |
| `esl_hyperexp`, `esl_stretchexp`, `esl_lognormal` | Niche distributions — not used by HMMER core |
| `esl_matrixops` | C array helpers — replaced by `matrix.zig` |
| `interface_gsl`, `interface_lapack` | External C library wrappers — pure Zig implementations instead |

## Design Decisions

| Aspect | Easel (C) | zeasel (Zig) |
|--------|-----------|--------------|
| Digital sequences | `ESL_DSQ` (uint8_t) + sentinel bytes | `[]u8` slice, no sentinels |
| Format dispatch | Function pointer vtable | Tagged union + switch |
| Memory management | malloc/realloc + manual free | `std.mem.Allocator` interface |
| Error handling | Macros (`ESL_FAIL`) | Zig error unions (`!T`) |
| Type-specific functions | Separate `DSet/FSet/ISet` | `comptime` generics |
| SIMD | Separate files per ISA | `@Vector` + `suggestVectorLength` |
| Hash tables | Custom `ESL_KEYHASH` | `std.StringHashMap` |

## License

See [LICENSE](LICENSE).

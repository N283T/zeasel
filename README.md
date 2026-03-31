# zeasel

A Zig reimplementation of [Easel](https://github.com/EddyRivasLab/easel), the C library for biological sequence analysis from the Eddy/Rivas laboratory. Easel underpins [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/). zeasel serves as the foundation for a Zig reimplementation of HMMER.

## Features

- **Biological alphabets** — DNA, RNA, and amino acid with comptime-generated lookup tables
- **Sequence I/O** — FASTA, Stockholm, GenBank, EMBL, Clustal, aligned FASTA
- **Multiple sequence alignment** — MSA operations, weighting (PB/BLOSUM), clustering, gap manipulation
- **SIMD vector operations** — comptime-generic `@Vector`-based ops with automatic ISA selection
- **Statistical distributions** — Gumbel, gamma, exponential, normal, Weibull, Dirichlet
- **Numerical optimization** — minimizer (conjugate gradient), root finder (Brent's method)
- **Score matrices** — BLOSUM62
- **Genetic code** — codon translation with NCBI genetic code tables
- **Sequence distance** — pairwise identity and Jukes-Cantor correction
- **Phylogenetic trees** — UPGMA construction, Newick I/O
- **SSI index** — sequence/subsequence index for fast random access
- **dsqdata** — binary packed sequence database for high-throughput search
- **WUSS notation** — RNA secondary structure parsing and validation

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
├── root.zig           # Public API re-exports
├── alphabet.zig       # Biological alphabets (DNA/RNA/Amino)
├── sequence.zig       # Sequence type
├── io.zig             # I/O module hub
│   └── io/            # Format-specific parsers/writers
│       ├── reader.zig, writer.zig
│       ├── fasta.zig, stockholm.zig
│       ├── genbank.zig, clustal.zig, afa.zig
├── msa.zig            # Multiple sequence alignment
├── msa_ops.zig        # MSA clustering and manipulation
├── msa_weight.zig     # MSA sequence weighting
├── simd.zig           # SIMD vector operations
├── stats.zig          # Statistical distributions
├── score_matrix.zig   # Substitution matrices (BLOSUM62)
├── gencode.zig        # Genetic code / translation
├── distance.zig       # Sequence distance calculations
├── tree.zig           # Phylogenetic trees (UPGMA, Newick)
├── matrix.zig         # Matrix operations
├── ssi.zig            # Sequence/Subsequence Index
├── dsqdata.zig        # Binary sequence database
├── wuss.zig           # RNA secondary structure (WUSS)
├── util.zig           # Utilities (RNG, random sequences)
└── tools/             # CLI tools
    ├── seqstat.zig
    ├── reformat.zig
    └── seqfetch.zig
```

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

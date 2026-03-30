# zeacel — Zig Reimplementation of Easel

## Overview

zeacel is a Zig reimplementation of [Easel](https://github.com/EddyRivasLab/easel), a C library for biological sequence analysis developed by the Eddy/Rivas laboratory at Harvard. Easel underpins HMMER and Infernal. zeacel serves as the foundation library for a Zig reimplementation of HMMER.

**Goals:**
- Feature parity with Easel, achieved incrementally
- Idiomatic Zig design (comptime generics, error unions, `@Vector` SIMD, allocator interface)
- No C API compatibility required — clean Zig-native interfaces
- No circular dependencies

**Non-goals:**
- Drop-in replacement for Easel's C API
- Backward compatibility with Easel's file formats for internal state (SSI index files etc.)

## Project Structure

```
zeacel/
├── build.zig
├── build.zig.zon
├── src/
│   ├── root.zig              # Public API re-exports
│   ├── alphabet.zig          # Biological alphabets (DNA/RNA/Amino)
│   ├── sequence.zig          # Single sequence object
│   ├── digital.zig           # Digital sequence (u8 encoded) utilities
│   ├── io/
│   │   ├── reader.zig        # Unified sequence reader (format dispatch)
│   │   ├── writer.zig        # Sequence writer
│   │   ├── fasta.zig         # FASTA format parser/writer
│   │   ├── stockholm.zig     # Stockholm format parser/writer
│   │   ├── genbank.zig       # GenBank/EMBL/DDBJ parser
│   │   └── buffer.zig        # I/O buffer abstraction
│   ├── msa.zig               # Multiple sequence alignment
│   ├── stats/
│   │   ├── gumbel.zig        # Gumbel distribution
│   │   ├── gamma.zig         # Gamma distribution
│   │   └── ...               # Other distributions
│   ├── simd/
│   │   ├── vector_ops.zig    # SIMD vector operations via @Vector
│   │   └── search.zig        # SIMD-accelerated sequence search
│   └── util/
│       └── random.zig        # Random number generator
└── test/
    ├── alphabet_test.zig
    ├── sequence_test.zig
    └── ...
```

## Dependency Graph (no cycles)

```
alphabet
    ↓
digital ← sequence
    ↓        ↓
  io/*     msa
    ↓
  simd/*   stats/*
```

## Key Design Decisions vs Easel

| Aspect | Easel (C) | zeacel (Zig) |
|--------|-----------|--------------|
| Digital sequences | `ESL_DSQ` (uint8_t) + sentinel bytes | `[]u8` slice, no sentinels |
| Text/Digital mode | Dual-mode coexistence | Digital-only default, text via conversion |
| Format dispatch | Function pointer vtable | Tagged union (`Format`) + switch |
| Memory management | malloc/realloc + manual free | `std.mem.Allocator` interface |
| Error handling | Macros (ESL_FAIL, ESL_XFAIL) | Zig error unions (`!T`) |
| Type-specific functions | DSet/FSet/ISet per type | `comptime` generics |
| SIMD | Separate files per ISA (sse, avx, neon) | `@Vector` + `suggestVectorLength` auto-dispatch |
| Object reuse | `esl_sq_Reuse()` pattern | Fresh alloc + deinit; arena/pool allocator externally |
| Hash tables | Custom `ESL_KEYHASH` | `std.StringHashMap` from stdlib |
| Data structures | Custom stack, red-black tree | Zig stdlib equivalents |

## Module Designs

### 1. Alphabet (`alphabet.zig`)

Biological sequence alphabets with comptime-generated lookup tables.

```zig
pub const AlphabetType = enum { dna, rna, amino };

pub const Alphabet = struct {
    type: AlphabetType,
    k: u8,                    // Canonical residue count (4 or 20)
    kp: u8,                   // Total symbol count (canonical + degenerate + gap + missing)
    symbols: []const u8,      // Symbol chars indexed by digital code
    encode_map: [128]u8,      // ASCII -> digital code (0xFF = invalid)
    degen: []const []const bool,  // Degeneracy matrix
    complement: ?[]const u8,      // Complement table (DNA/RNA only)

    pub fn encode(self: Alphabet, char: u8) !u8;
    pub fn decode(self: Alphabet, code: u8) u8;
    pub fn digitize(self: Alphabet, allocator: Allocator, text: []const u8) ![]u8;
    pub fn textize(self: Alphabet, allocator: Allocator, dsq: []const u8) ![]u8;
    pub fn reverseComplement(self: Alphabet, dsq: []u8) !void;
    pub fn isCanonical(self: Alphabet, code: u8) bool;
    pub fn isGap(self: Alphabet, code: u8) bool;
};

pub const dna = comptime initDna();
pub const rna = comptime initRna();
pub const amino = comptime initAmino();
```

Standard alphabets are constructed at comptime with zero runtime allocation.

### 2. Sequence (`sequence.zig`)

Digital-only sequence object. Text representation available via conversion.

```zig
pub const Sequence = struct {
    name: []const u8,
    accession: ?[]const u8 = null,
    description: ?[]const u8 = null,
    taxonomy_id: ?i32 = null,
    dsq: []u8,
    secondary_structure: ?[]const u8 = null,
    source: ?Source = null,
    abc: *const Alphabet,
    allocator: Allocator,

    pub const Source = struct {
        name: []const u8,
        start: i64,       // 1-indexed; start > end = reverse strand
        end: i64,
        full_length: i64,
    };

    pub fn init(allocator: Allocator, abc: *const Alphabet, name: []const u8, dsq: []const u8) !Sequence;
    pub fn fromText(allocator: Allocator, abc: *const Alphabet, name: []const u8, text: []const u8) !Sequence;
    pub fn clone(self: Sequence) !Sequence;
    pub fn toText(self: Sequence, allocator: Allocator) ![]u8;
    pub fn reverseComplement(self: *Sequence) !void;
    pub fn subseq(self: Sequence, allocator: Allocator, start: usize, end: usize) !Sequence;
    pub fn len(self: Sequence) usize;
    pub fn deinit(self: *Sequence) void;
};
```

Disk offsets and I/O state are kept in the reader, not in Sequence.

### 3. I/O (`io/`)

Unified reader with tagged union format dispatch and iterator pattern.

```zig
// io/reader.zig
pub const Format = enum {
    fasta, stockholm, genbank, embl, uniprot,
    pub fn detect(header: []const u8) ?Format;
};

pub const Reader = struct {
    format: Format,
    source: Source,
    abc: *const Alphabet,
    allocator: Allocator,

    pub const Source = union(enum) { file: FileSource, memory: MemorySource, stream: StreamSource };

    pub fn open(allocator: Allocator, abc: *const Alphabet, path: []const u8, format: ?Format) !Reader;
    pub fn fromMemory(allocator: Allocator, abc: *const Alphabet, data: []const u8, format: Format) !Reader;
    pub fn next(self: *Reader) !?Sequence;
    pub fn deinit(self: *Reader) void;
};

// io/writer.zig
pub const Writer = struct {
    format: Format,
    dest: std.io.AnyWriter,

    pub fn init(dest: std.io.AnyWriter, format: Format) Writer;
    pub fn write(self: *Writer, seq: Sequence) !void;
    pub fn writeAll(self: *Writer, seqs: []const Sequence) !void;
};
```

Usage pattern:
```zig
var reader = try Reader.open(allocator, &alphabet.dna, "sequences.fasta", null);
defer reader.deinit();
while (try reader.next()) |seq| {
    defer seq.deinit();
    // process...
}
```

Format-specific parsing lives in `io/fasta.zig`, `io/stockholm.zig`, etc.

SSI-equivalent indexing, compressed file support, and caching are deferred to later phases.

### 4. MSA (`msa.zig`)

Multiple sequence alignment with immutable operations and Stockholm markup preservation.

```zig
pub const Msa = struct {
    names: [][]const u8,
    seqs: [][]u8,               // Digital sequences [nseq][alen]
    weights: ?[]f64 = null,
    alen: usize,
    abc: *const Alphabet,
    allocator: Allocator,

    // Optional metadata
    name: ?[]const u8 = null,
    accession: ?[]const u8 = null,
    description: ?[]const u8 = null,
    consensus_ss: ?[]const u8 = null,
    reference: ?[]const u8 = null,
    seq_accessions: ?[]?[]const u8 = null,
    seq_descriptions: ?[]?[]const u8 = null,
    seq_ss: ?[]?[]const u8 = null,

    // Stockholm generic markup (preserved for round-trip fidelity)
    gf_tags: std.StringArrayHashMap([]const u8),
    gc_tags: std.StringArrayHashMap([]const u8),
    gs_tags: std.StringArrayHashMap([][]const u8),
    gr_tags: std.StringArrayHashMap([][]const u8),

    pub fn nseq(self: Msa) usize;
    pub fn extractSeq(self: Msa, allocator: Allocator, idx: usize) !Sequence;
    pub fn minimizeGaps(self: *Msa) !void;
    pub fn removeAllGaps(self: *Msa) !void;
    pub fn selectSeqs(self: Msa, allocator: Allocator, indices: []const usize) !Msa;
    pub fn selectColumns(self: Msa, allocator: Allocator, mask: []const bool) !Msa;
    pub fn fromText(allocator: Allocator, abc: *const Alphabet, names: [][]const u8, text_seqs: [][]const u8) !Msa;
    pub fn checksum(self: Msa) u32;
    pub fn deinit(self: *Msa) void;
};
```

Operations like `selectSeqs` and `selectColumns` return new Msa instances (immutable pattern).
Parser builds alignment incrementally with ArrayList, then freezes into Msa.

### 5. SIMD / Vector Operations (`simd/vector_ops.zig`)

Comptime generic vector operations with automatic SIMD via `@Vector`.

```zig
pub fn VectorOps(comptime T: type) type {
    return struct {
        // Element-wise
        pub fn set(dst: []T, value: T) void;
        pub fn scale(dst: []T, factor: T) void;
        pub fn add(dst: []T, src: []const T) void;
        pub fn addScaled(dst: []T, src: []const T, a: T) void;

        // Reductions
        pub fn sum(v: []const T) T;
        pub fn dot(a: []const T, b: []const T) T;
        pub fn max(v: []const T) T;
        pub fn min(v: []const T) T;
        pub fn argmax(v: []const T) usize;
        pub fn argmin(v: []const T) usize;

        // Normalization (float types)
        pub fn normalize(v: []T) void;
        pub fn logNormalize(v: []T) void;

        // Log-space (float types)
        pub fn logSum(v: []const T) T;

        // Information theory (float types)
        pub fn entropy(p: []const T) T;
        pub fn relativeEntropy(p: []const T, q: []const T) T;

        // Utilities
        pub fn copy(dst: []T, src: []const T) void;
        pub fn reverse(v: []T) void;
        pub fn sort(v: []T, order: enum { ascending, descending }) void;
    };
}

pub const VecF64 = VectorOps(f64);
pub const VecF32 = VectorOps(f32);
pub const VecI32 = VectorOps(i32);
pub const VecU8  = VectorOps(u8);
```

SIMD implementation uses `@Vector` with `std.simd.suggestVectorLength` for automatic ISA selection.
Replaces Easel's ~2,200 lines of type-duplicated code with a single comptime-generic module.

### 6. Utilities

- **Random**: Wrap `std.Random` or implement Mersenne Twister for Easel compatibility
- **Hash tables**: Use `std.StringHashMap` directly (replaces `ESL_KEYHASH`)
- **Stack, trees**: Use Zig stdlib equivalents (replaces `esl_stack`, `esl_red_black`)

No custom data structure wrappers unless stdlib proves insufficient.

## Testing Strategy

- Inline `test` blocks in each module, run via `zig test`
- Reuse Easel's test data files (sequence files from `testsuite/`, `esl_msa_testfiles/`)
- Key test categories:
  - **Round-trip**: encode→decode, digitize→textize, read→write→read
  - **Edge cases**: empty sequences, single residue, maximum length
  - **Error paths**: invalid characters, malformed files, truncated input
  - **SIMD correctness**: scalar vs SIMD result equivalence
  - **Cross-alphabet**: verify DNA/RNA/Amino specific behavior

## Implementation Phases

| Phase | Modules | Deliverable | Easel Reference |
|-------|---------|-------------|-----------------|
| 1 | `alphabet`, `digital` | DNA/RNA/Amino alphabets, encode/decode | `esl_alphabet.*` |
| 2 | `sequence` | Sequence type, fromText, clone, subseq | `esl_sq.*` |
| 3 | `io/buffer`, `io/fasta` | FASTA read/write | `esl_buffer.*`, `esl_sqio_ascii.*` |
| 4 | `io/stockholm`, `msa` | Stockholm read/write, MSA operations | `esl_msa.*`, `esl_msafile_stockholm.*` |
| 5 | `simd/vector_ops` | Comptime generic vector ops + SIMD | `esl_vectorops.*` |
| 6 | `stats/*` | Probability distributions (Gumbel, Gamma, etc.) | `esl_gumbel.*`, `esl_gamma.*`, etc. |
| 7 | `io/genbank` + other formats | GenBank/EMBL/DDBJ/UniProt parsers | `esl_sqio_ascii.*` |
| 8 | Remaining utilities + miniapps | CLI tools | `miniapps/*` |

Each phase delivers tested, working code committed independently.

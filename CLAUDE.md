# CLAUDE.md — zeasel

## Project Overview

zeasel is a Zig reimplementation of [Easel](https://github.com/EddyRivasLab/easel), the biological sequence analysis library from the Eddy/Rivas lab. It serves as the foundation for a Zig reimplementation of HMMER.

## Build & Test

```bash
zig build          # Build library
zig build test     # Run all tests
zig build tools    # Build CLI tools (zeasel-seqstat, zeasel-reformat, zeasel-seqfetch)
```

Minimum Zig version: **0.15.0**

## Architecture

### Module dependency graph (no cycles)

```
alphabet
   ↓
sequence ← gencode, score_matrix
   ↓
io (fasta, stockholm, genbank, clustal, afa)
   ↓
msa ← msa_weight, msa_ops, distance
   ↓
ssi, dsqdata, tree, matrix, wuss
   ↓
simd (vector_ops)
   ↓
stats (gumbel, gamma, exponential, normal, weibull, dirichlet, histogram, minimizer, rootfinder)
   ↓
util (random, random_seq)
```

### Key patterns

- **Comptime alphabets**: `alphabet.dna`, `alphabet.rna`, `alphabet.amino` are fully constructed at comptime. No runtime allocation for standard alphabets.
- **Digital-only sequences**: Sequences store digitally encoded `[]u8`. Text representation via `alphabet.textize()`.
- **Tagged union format dispatch**: `io.Format` enum drives reader/writer selection via switch, not vtables.
- **Immutable MSA operations**: `selectSeqs`, `selectColumns` return new `Msa` instances.
- **Comptime generics for SIMD**: `simd.VectorOps(T)` replaces Easel's per-type duplicated code.
- **Allocator interface**: All allocating functions take `std.mem.Allocator`. No global state.

### Public API

All public modules are re-exported from `src/root.zig`:

```zig
const zeasel = @import("zeasel");
// zeasel.alphabet, zeasel.sequence, zeasel.io, zeasel.msa, etc.
```

## Conventions

- Tests are inline `test` blocks within each module, not in separate files.
- Each module is a single `.zig` file. Sub-modules use directories (e.g., `io/`, `stats/`, `simd/`, `util/`).
- CLI tools live in `src/tools/` and are registered in `build.zig`.
- No sentinels on digital sequences (unlike Easel's `ESL_DSQ`).
- Error handling uses Zig error unions (`!T`), not status codes.

## Adding a New Module

1. Create `src/<module>.zig`
2. Add `pub const <module> = @import("<module>.zig");` to `src/root.zig`
3. Add `_ = <module>;` to the `test` block in `src/root.zig`
4. If it's a CLI tool, add it to the `tool_defs` array in `build.zig`

## Adding a New I/O Format

1. Create `src/io/<format>.zig` with `parse()` and optionally `write()` functions
2. Add the format variant to `io.Format` enum in `src/io/reader.zig`
3. Add the dispatch case in `Reader` and `Writer`
4. Add `Format.detect()` signature recognition if the format has a distinctive header

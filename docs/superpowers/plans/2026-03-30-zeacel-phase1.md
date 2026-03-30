# zeacel Phase 1: Project Scaffold + Alphabet Module

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Set up the zeacel Zig project and implement the Alphabet module — the foundational building block for biological sequence digitization.

**Architecture:** The Alphabet module encodes biological sequences (DNA, RNA, amino acids) as compact u8 digital codes using comptime-generated lookup tables. Standard alphabets (dna, rna, amino) are fully constructed at compile time with zero runtime allocation. The encoding follows Easel's internal ordering: canonical residues (0..K-1), gap (K), degeneracies (K+1..Kp-4), any (Kp-3), nonresidue (Kp-2), missing (Kp-1).

**Tech Stack:** Zig 0.15.2, `zig build` / `zig test`

**Reference:** Easel source at `/Users/nagaet/ghq/github.com/EddyRivasLab/easel/` — specifically `esl_alphabet.h` and `esl_alphabet.c`.

---

### Task 1: Initialize Zig Project

**Files:**
- Create: `build.zig`
- Create: `build.zig.zon`
- Create: `src/root.zig`

- [ ] **Step 1: Create build.zig.zon**

```zig
.{
    .name = .@"zeacel",
    .version = .@"0.1.0",
    .fingerprint = .@"zig-default-fingerprint",
    .minimum_zig_version = .@"0.15.0",
    .paths = .{
        "build.zig",
        "build.zig.zon",
        "src",
    },
}
```

- [ ] **Step 2: Create build.zig**

```zig
const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const lib_mod = b.addModule("zeacel", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const lib = b.addLibrary(.{
        .linkage = .static,
        .name = "zeacel",
        .root_module = lib_mod,
    });
    b.installArtifact(lib);

    const lib_unit_tests = b.addTest(.{
        .root_module = lib_mod,
    });
    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);
}
```

- [ ] **Step 3: Create src/root.zig**

```zig
pub const alphabet = @import("alphabet.zig");
```

This will fail to compile until alphabet.zig exists — that is expected and confirms the module wiring.

- [ ] **Step 4: Verify build setup compiles (expect failure)**

Run: `cd /Users/nagaet/zeacel && zig build 2>&1`
Expected: Compilation error about missing `alphabet.zig`. This confirms `build.zig` wiring is correct.

- [ ] **Step 5: Create empty alphabet.zig to unblock build**

Create `src/alphabet.zig` with a placeholder:

```zig
// Biological sequence alphabets (DNA, RNA, Amino Acid).
```

- [ ] **Step 6: Verify clean build**

Run: `cd /Users/nagaet/zeacel && zig build`
Expected: Success with no errors.

- [ ] **Step 7: Initialize git and commit**

```bash
cd /Users/nagaet/zeacel
git init
git add build.zig build.zig.zon src/root.zig src/alphabet.zig
git commit -m "chore: initialize zeacel Zig project"
```

---

### Task 2: Alphabet Core Types and Constants

**Files:**
- Modify: `src/alphabet.zig`

This task defines the core data types that all other tasks depend on. It follows Easel's internal symbol ordering:
- Positions `0..K-1`: Canonical residues (e.g., A=0, C=1, G=2, T=3 for DNA)
- Position `K`: Gap symbol (`-`)
- Positions `K+1..Kp-4`: Degeneracy symbols (e.g., R, Y, M, K, S, W, H, B, V, D for DNA)
- Position `Kp-3`: Any/unknown symbol (N for DNA, X for amino)
- Position `Kp-2`: Nonresidue symbol (`*`)
- Position `Kp-1`: Missing data symbol (`~`)

- [ ] **Step 1: Write tests for core type properties**

Add to `src/alphabet.zig`:

```zig
const std = @import("std");

pub const AlphabetType = enum {
    dna,
    rna,
    amino,
};

/// Sentinel value for invalid/unmapped ASCII characters in encode_map.
pub const INVALID_CODE: u8 = 0xFF;

pub const Alphabet = struct {
    kind: AlphabetType,
    /// Number of canonical residues (4 for DNA/RNA, 20 for amino).
    k: u8,
    /// Total symbol count: canonical + gap + degeneracies + any + nonresidue + missing.
    kp: u8,
    /// Symbol characters indexed by digital code.
    /// Example for DNA: "ACGT-RYMKSWHBVDN*~"
    symbols: []const u8,
    /// ASCII (0..127) -> digital code lookup table. INVALID_CODE for unmapped characters.
    encode_map: [128]u8,
    /// Complement table for DNA/RNA. complement[code] = complementary code.
    /// null for amino acid alphabets.
    complement: ?[]const u8,

    /// Encode a single ASCII character to its digital code.
    pub fn encode(self: *const Alphabet, char: u8) error{InvalidCharacter}!u8 {
        if (char > 127) return error.InvalidCharacter;
        const code = self.encode_map[char];
        if (code == INVALID_CODE) return error.InvalidCharacter;
        return code;
    }

    /// Decode a digital code back to its ASCII character.
    pub fn decode(self: *const Alphabet, code: u8) u8 {
        return self.symbols[code];
    }

    /// Check if a digital code represents a canonical residue (not gap, not degenerate).
    pub fn isCanonical(self: *const Alphabet, code: u8) bool {
        return code < self.k;
    }

    /// Check if a digital code is the gap symbol.
    pub fn isGap(self: *const Alphabet, code: u8) bool {
        return code == self.k;
    }

    /// Check if a digital code is a degeneracy symbol (excluding the "any" symbol).
    pub fn isDegenerate(self: *const Alphabet, code: u8) bool {
        return code > self.k and code < self.kp - 2;
    }

    /// Check if a digital code is the "any" / unknown symbol (N for DNA, X for amino).
    pub fn isUnknown(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 3;
    }

    /// Check if a digital code is the nonresidue symbol (*).
    pub fn isNonresidue(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 2;
    }

    /// Check if a digital code is the missing data symbol (~).
    pub fn isMissing(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 1;
    }

    /// Get the digital code for the gap symbol.
    pub fn gapCode(self: *const Alphabet) u8 {
        return self.k;
    }

    /// Get the digital code for the "any" / unknown symbol.
    pub fn unknownCode(self: *const Alphabet) u8 {
        return self.kp - 3;
    }
};

test "Alphabet struct size constants" {
    // DNA: K=4, Kp=18 (ACGT + gap + 10 degeneracies + any + nonresidue + missing)
    const a = Alphabet{
        .kind = .dna,
        .k = 4,
        .kp = 18,
        .symbols = "ACGT-RYMKSWHBVDN*~",
        .encode_map = [_]u8{INVALID_CODE} ** 128,
        .complement = null,
    };
    try std.testing.expectEqual(@as(u8, 4), a.k);
    try std.testing.expectEqual(@as(u8, 18), a.kp);
}

test "Alphabet classification methods" {
    const a = Alphabet{
        .kind = .dna,
        .k = 4,
        .kp = 18,
        .symbols = "ACGT-RYMKSWHBVDN*~",
        .encode_map = [_]u8{INVALID_CODE} ** 128,
        .complement = null,
    };

    // Canonical: 0..3 (A, C, G, T)
    try std.testing.expect(a.isCanonical(0));
    try std.testing.expect(a.isCanonical(3));
    try std.testing.expect(!a.isCanonical(4));

    // Gap: position 4
    try std.testing.expect(a.isGap(4));
    try std.testing.expect(!a.isGap(0));

    // Degenerate: positions 5..13 (R, Y, M, K, S, W, H, B, V, D)
    try std.testing.expect(a.isDegenerate(5));
    try std.testing.expect(a.isDegenerate(13));
    try std.testing.expect(!a.isDegenerate(4));  // gap
    try std.testing.expect(!a.isDegenerate(15)); // any (N)

    // Unknown/Any: position Kp-3 = 15 (N)
    try std.testing.expect(a.isUnknown(15));
    try std.testing.expect(!a.isUnknown(0));

    // Nonresidue: position Kp-2 = 16 (*)
    try std.testing.expect(a.isNonresidue(16));

    // Missing: position Kp-1 = 17 (~)
    try std.testing.expect(a.isMissing(17));
}
```

- [ ] **Step 2: Run tests**

Run: `cd /Users/nagaet/zeacel && zig build test`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
cd /Users/nagaet/zeacel
git add src/alphabet.zig
git commit -m "feat: add Alphabet core types with classification methods"
```

---

### Task 3: Comptime DNA Alphabet Construction

**Files:**
- Modify: `src/alphabet.zig`

This task builds the DNA alphabet entirely at comptime. The encode_map maps ASCII characters to digital codes, handling case-insensitivity and synonyms (U→T, X→N, _→-, .→-). The complement table maps each code to its Watson-Crick complement.

- [ ] **Step 1: Write tests for DNA alphabet**

Add to `src/alphabet.zig`:

```zig
/// Build an encode_map at comptime from a symbol string.
/// Maps each character in symbols to its index, plus case-insensitive variants.
/// Additional synonyms can be specified as pairs: {from, to_symbol}.
fn buildEncodeMap(
    comptime symbols: []const u8,
    comptime synonyms: []const [2]u8,
) [128]u8 {
    var map = [_]u8{INVALID_CODE} ** 128;

    // Map each symbol to its index
    for (symbols, 0..) |sym, i| {
        map[sym] = @intCast(i);
    }

    // Case-insensitive: map lowercase of each uppercase letter
    for (symbols, 0..) |sym, i| {
        if (sym >= 'A' and sym <= 'Z') {
            map[sym + 32] = @intCast(i); // lowercase
        }
    }

    // Apply synonyms
    for (synonyms) |pair| {
        const from = pair[0];
        const to_sym = pair[1];
        // Find the code for the target symbol
        const code = map[to_sym];
        map[from] = code;
        // Also map lowercase of synonym
        if (from >= 'A' and from <= 'Z') {
            map[from + 32] = code;
        }
    }

    return map;
}

const dna_symbols = "ACGT-RYMKSWHBVDN*~";
const dna_complement_table = [18]u8{
    3,  // A(0) -> T(3)
    2,  // C(1) -> G(2)
    1,  // G(2) -> C(1)
    0,  // T(3) -> A(0)
    4,  // -(4) -> -(4)
    6,  // R(5) -> Y(6)   (AG -> CT)
    5,  // Y(6) -> R(5)   (CT -> AG)
    8,  // M(7) -> K(8)   (AC -> GT)
    7,  // K(8) -> M(7)   (GT -> AC)
    9,  // S(9) -> S(9)   (CG -> CG)
    10, // W(10)-> W(10)  (AT -> AT)
    14, // H(11)-> D(14)  (ACT -> AGT)
    13, // B(12)-> V(13)  (CGT -> ACG)
    12, // V(13)-> B(12)  (ACG -> CGT)
    11, // D(14)-> H(11)  (AGT -> ACT)
    15, // N(15)-> N(15)
    16, // *(16)-> *(16)
    17, // ~(17)-> ~(17)
};

pub const dna = Alphabet{
    .kind = .dna,
    .k = 4,
    .kp = 18,
    .symbols = dna_symbols,
    .encode_map = comptime buildEncodeMap(dna_symbols, &.{
        .{ 'U', 'T' }, // read U as T
        .{ 'X', 'N' }, // read X as N
        .{ 'I', 'A' }, // Inosine as A
        .{ '_', '-' }, // underscore as gap
        .{ '.', '-' }, // dot as gap
    }),
    .complement = &dna_complement_table,
};

test "DNA alphabet encode canonical residues" {
    try std.testing.expectEqual(@as(u8, 0), try dna.encode('A'));
    try std.testing.expectEqual(@as(u8, 1), try dna.encode('C'));
    try std.testing.expectEqual(@as(u8, 2), try dna.encode('G'));
    try std.testing.expectEqual(@as(u8, 3), try dna.encode('T'));
}

test "DNA alphabet encode is case-insensitive" {
    try std.testing.expectEqual(@as(u8, 0), try dna.encode('a'));
    try std.testing.expectEqual(@as(u8, 1), try dna.encode('c'));
    try std.testing.expectEqual(@as(u8, 2), try dna.encode('g'));
    try std.testing.expectEqual(@as(u8, 3), try dna.encode('t'));
}

test "DNA alphabet encode synonyms" {
    // U -> T
    try std.testing.expectEqual(@as(u8, 3), try dna.encode('U'));
    try std.testing.expectEqual(@as(u8, 3), try dna.encode('u'));
    // X -> N
    try std.testing.expectEqual(@as(u8, 15), try dna.encode('X'));
    try std.testing.expectEqual(@as(u8, 15), try dna.encode('x'));
    // Gap synonyms
    try std.testing.expectEqual(@as(u8, 4), try dna.encode('_'));
    try std.testing.expectEqual(@as(u8, 4), try dna.encode('.'));
}

test "DNA alphabet encode degeneracies" {
    try std.testing.expectEqual(@as(u8, 5), try dna.encode('R'));  // purine (AG)
    try std.testing.expectEqual(@as(u8, 6), try dna.encode('Y'));  // pyrimidine (CT)
    try std.testing.expectEqual(@as(u8, 15), try dna.encode('N')); // any
}

test "DNA alphabet encode invalid character" {
    try std.testing.expectError(error.InvalidCharacter, dna.encode('!'));
    try std.testing.expectError(error.InvalidCharacter, dna.encode('Z'));
    try std.testing.expectError(error.InvalidCharacter, dna.encode(0));
}

test "DNA alphabet decode" {
    try std.testing.expectEqual(@as(u8, 'A'), dna.decode(0));
    try std.testing.expectEqual(@as(u8, 'C'), dna.decode(1));
    try std.testing.expectEqual(@as(u8, 'G'), dna.decode(2));
    try std.testing.expectEqual(@as(u8, 'T'), dna.decode(3));
    try std.testing.expectEqual(@as(u8, '-'), dna.decode(4));
    try std.testing.expectEqual(@as(u8, 'N'), dna.decode(15));
}

test "DNA alphabet complement" {
    const comp = dna.complement.?;
    try std.testing.expectEqual(@as(u8, 3), comp[0]); // A -> T
    try std.testing.expectEqual(@as(u8, 2), comp[1]); // C -> G
    try std.testing.expectEqual(@as(u8, 1), comp[2]); // G -> C
    try std.testing.expectEqual(@as(u8, 0), comp[3]); // T -> A
    try std.testing.expectEqual(@as(u8, 6), comp[5]); // R -> Y
    try std.testing.expectEqual(@as(u8, 5), comp[6]); // Y -> R
}
```

- [ ] **Step 2: Run tests**

Run: `cd /Users/nagaet/zeacel && zig build test`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
cd /Users/nagaet/zeacel
git add src/alphabet.zig
git commit -m "feat: add comptime DNA alphabet with encode/decode/complement"
```

---

### Task 4: RNA and Amino Acid Alphabets

**Files:**
- Modify: `src/alphabet.zig`

RNA is nearly identical to DNA (U instead of T, same degeneracies adapted for U). Amino acid has K=20, Kp=29. Symbols: `ACDEFGHIKLMNPQRSTVWY-BJZOUX*~`. Degeneracies: B=ND, J=IL, Z=QE, U=C (selenocysteine), O=K (pyrrolysine).

- [ ] **Step 1: Write RNA and amino acid alphabets with tests**

Add to `src/alphabet.zig`:

```zig
const rna_symbols = "ACGU-RYMKSWHBVDN*~";
const rna_complement_table = [18]u8{
    3,  // A(0) -> U(3)
    2,  // C(1) -> G(2)
    1,  // G(2) -> C(1)
    0,  // U(3) -> A(0)
    4,  // -(4) -> -(4)
    6,  // R(5) -> Y(6)
    5,  // Y(6) -> R(5)
    8,  // M(7) -> K(8)
    7,  // K(8) -> M(7)
    9,  // S(9) -> S(9)
    10, // W(10)-> W(10)
    14, // H(11)-> D(14)
    13, // B(12)-> V(13)
    12, // V(13)-> B(12)
    11, // D(14)-> H(11)
    15, // N(15)-> N(15)
    16, // *(16)-> *(16)
    17, // ~(17)-> ~(17)
};

pub const rna = Alphabet{
    .kind = .rna,
    .k = 4,
    .kp = 18,
    .symbols = rna_symbols,
    .encode_map = comptime buildEncodeMap(rna_symbols, &.{
        .{ 'T', 'U' }, // read T as U
        .{ 'X', 'N' }, // read X as N
        .{ 'I', 'A' }, // Inosine as A
        .{ '_', '-' }, // underscore as gap
        .{ '.', '-' }, // dot as gap
    }),
    .complement = &rna_complement_table,
};

const amino_symbols = "ACDEFGHIKLMNPQRSTVWY-BJZOUX*~";

pub const amino = Alphabet{
    .kind = .amino,
    .k = 20,
    .kp = 29,
    .symbols = amino_symbols,
    .encode_map = comptime buildEncodeMap(amino_symbols, &.{
        .{ '_', '-' }, // underscore as gap
        .{ '.', '-' }, // dot as gap
    }),
    .complement = null,
};

test "RNA alphabet encode" {
    try std.testing.expectEqual(@as(u8, 0), try rna.encode('A'));
    try std.testing.expectEqual(@as(u8, 3), try rna.encode('U'));
    // T -> U synonym
    try std.testing.expectEqual(@as(u8, 3), try rna.encode('T'));
    try std.testing.expectEqual(@as(u8, 3), try rna.encode('t'));
    // Case-insensitive
    try std.testing.expectEqual(@as(u8, 0), try rna.encode('a'));
}

test "RNA alphabet complement" {
    const comp = rna.complement.?;
    try std.testing.expectEqual(@as(u8, 3), comp[0]); // A -> U
    try std.testing.expectEqual(@as(u8, 0), comp[3]); // U -> A
}

test "Amino acid alphabet encode" {
    // First residue: A=0
    try std.testing.expectEqual(@as(u8, 0), try amino.encode('A'));
    // Last canonical residue: Y=19
    try std.testing.expectEqual(@as(u8, 19), try amino.encode('Y'));
    // Case insensitive
    try std.testing.expectEqual(@as(u8, 0), try amino.encode('a'));
    // Gap
    try std.testing.expectEqual(@as(u8, 20), try amino.encode('-'));
    try std.testing.expectEqual(@as(u8, 20), try amino.encode('_'));
    try std.testing.expectEqual(@as(u8, 20), try amino.encode('.'));
}

test "Amino acid alphabet sizes" {
    try std.testing.expectEqual(@as(u8, 20), amino.k);
    try std.testing.expectEqual(@as(u8, 29), amino.kp);
    // X is the "any" symbol at position Kp-3 = 26
    try std.testing.expectEqual(@as(u8, 26), amino.kp - 3);
    try std.testing.expect(amino.isUnknown(26));
    try std.testing.expectEqual(@as(u8, 'X'), amino.decode(26));
}

test "Amino acid alphabet has no complement" {
    try std.testing.expect(amino.complement == null);
}

test "Amino acid alphabet degeneracy codes" {
    // B (Asx) = position 21
    try std.testing.expectEqual(@as(u8, 21), try amino.encode('B'));
    try std.testing.expect(amino.isDegenerate(21));
    // J (Xle) = position 22
    try std.testing.expectEqual(@as(u8, 22), try amino.encode('J'));
    // Z (Glx) = position 23
    try std.testing.expectEqual(@as(u8, 23), try amino.encode('Z'));
    // O (pyrrolysine) = position 24
    try std.testing.expectEqual(@as(u8, 24), try amino.encode('O'));
    // U (selenocysteine) = position 25
    try std.testing.expectEqual(@as(u8, 25), try amino.encode('U'));
}
```

- [ ] **Step 2: Run tests**

Run: `cd /Users/nagaet/zeacel && zig build test`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
cd /Users/nagaet/zeacel
git add src/alphabet.zig
git commit -m "feat: add RNA and amino acid alphabets"
```

---

### Task 5: Digitize and Textize Functions

**Files:**
- Modify: `src/alphabet.zig`

These functions convert between text sequences (e.g., "ACGTACGT") and digital sequences (e.g., [0,1,2,3,0,1,2,3]). They allocate memory via the provided allocator.

- [ ] **Step 1: Write tests for digitize and textize**

Add to `src/alphabet.zig`:

```zig
    /// Digitize an entire text sequence.
    /// Caller owns the returned slice and must free it with the same allocator.
    pub fn digitize(self: *const Alphabet, allocator: std.mem.Allocator, text: []const u8) ![]u8 {
        const dsq = try allocator.alloc(u8, text.len);
        for (text, 0..) |char, i| {
            dsq[i] = self.encode(char) catch return error.InvalidCharacter;
        }
        return dsq;
    }

    /// Convert a digital sequence back to text.
    /// Caller owns the returned slice and must free it with the same allocator.
    pub fn textize(self: *const Alphabet, allocator: std.mem.Allocator, dsq: []const u8) ![]u8 {
        const text = try allocator.alloc(u8, dsq.len);
        for (dsq, 0..) |code, i| {
            text[i] = self.decode(code);
        }
        return text;
    }

    /// Reverse complement a digital sequence in place.
    /// Returns error if the alphabet has no complement table (e.g., amino acids).
    pub fn reverseComplement(self: *const Alphabet, dsq: []u8) error{NoComplement}!void {
        const comp = self.complement orelse return error.NoComplement;
        // Reverse and complement in a single pass
        var i: usize = 0;
        var j: usize = dsq.len;
        while (i < j) {
            j -= 1;
            const tmp = comp[dsq[i]];
            dsq[i] = comp[dsq[j]];
            dsq[j] = tmp;
            i += 1;
        }
        // If odd length, the middle element needs complementing (was swapped with itself)
        if (dsq.len % 2 == 1) {
            const mid = dsq.len / 2;
            dsq[mid] = comp[dsq[mid]];
        }
    }
```

Add tests:

```zig
test "digitize DNA sequence" {
    const allocator = std.testing.allocator;
    const dsq = try dna.digitize(allocator, "ACGT");
    defer allocator.free(dsq);

    try std.testing.expectEqual(@as(u8, 0), dsq[0]); // A
    try std.testing.expectEqual(@as(u8, 1), dsq[1]); // C
    try std.testing.expectEqual(@as(u8, 2), dsq[2]); // G
    try std.testing.expectEqual(@as(u8, 3), dsq[3]); // T
}

test "digitize handles lowercase" {
    const allocator = std.testing.allocator;
    const dsq = try dna.digitize(allocator, "acgt");
    defer allocator.free(dsq);

    try std.testing.expectEqual(@as(u8, 0), dsq[0]);
    try std.testing.expectEqual(@as(u8, 1), dsq[1]);
    try std.testing.expectEqual(@as(u8, 2), dsq[2]);
    try std.testing.expectEqual(@as(u8, 3), dsq[3]);
}

test "digitize with degeneracies" {
    const allocator = std.testing.allocator;
    const dsq = try dna.digitize(allocator, "ACNR");
    defer allocator.free(dsq);

    try std.testing.expectEqual(@as(u8, 0), dsq[0]);  // A
    try std.testing.expectEqual(@as(u8, 1), dsq[1]);  // C
    try std.testing.expectEqual(@as(u8, 15), dsq[2]); // N (any)
    try std.testing.expectEqual(@as(u8, 5), dsq[3]);  // R (purine)
}

test "digitize invalid character returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(error.InvalidCharacter, dna.digitize(allocator, "ACGZ"));
}

test "digitize empty sequence" {
    const allocator = std.testing.allocator;
    const dsq = try dna.digitize(allocator, "");
    defer allocator.free(dsq);
    try std.testing.expectEqual(@as(usize, 0), dsq.len);
}

test "textize round-trip" {
    const allocator = std.testing.allocator;
    const original = "ACGTRYNSWHBVD";
    const dsq = try dna.digitize(allocator, original);
    defer allocator.free(dsq);
    const text = try dna.textize(allocator, dsq);
    defer allocator.free(text);

    try std.testing.expectEqualStrings(original, text);
}

test "textize amino acid round-trip" {
    const allocator = std.testing.allocator;
    const original = "ACDEFGHIKLMNPQRSTVWY";
    const dsq = try amino.digitize(allocator, original);
    defer allocator.free(dsq);
    const text = try amino.textize(allocator, dsq);
    defer allocator.free(text);

    try std.testing.expectEqualStrings(original, text);
}

test "reverseComplement DNA" {
    const allocator = std.testing.allocator;
    // ACGT -> reverse: TGCA -> complement: ACGT
    var dsq = try dna.digitize(allocator, "ACGT");
    defer allocator.free(dsq);

    try dna.reverseComplement(dsq);
    const text = try dna.textize(allocator, dsq);
    defer allocator.free(text);

    try std.testing.expectEqualStrings("ACGT", text);
}

test "reverseComplement DNA asymmetric" {
    const allocator = std.testing.allocator;
    // AACG -> reverse: GCAA -> complement: CGTT
    var dsq = try dna.digitize(allocator, "AACG");
    defer allocator.free(dsq);

    try dna.reverseComplement(dsq);
    const text = try dna.textize(allocator, dsq);
    defer allocator.free(text);

    try std.testing.expectEqualStrings("CGTT", text);
}

test "reverseComplement odd length" {
    const allocator = std.testing.allocator;
    // ACG -> reverse: GCA -> complement: CGT
    var dsq = try dna.digitize(allocator, "ACG");
    defer allocator.free(dsq);

    try dna.reverseComplement(dsq);
    const text = try dna.textize(allocator, dsq);
    defer allocator.free(text);

    try std.testing.expectEqualStrings("CGT", text);
}

test "reverseComplement single base" {
    const allocator = std.testing.allocator;
    var dsq = try dna.digitize(allocator, "A");
    defer allocator.free(dsq);

    try dna.reverseComplement(dsq);
    const text = try dna.textize(allocator, dsq);
    defer allocator.free(text);

    try std.testing.expectEqualStrings("T", text);
}

test "reverseComplement empty" {
    const allocator = std.testing.allocator;
    var dsq = try dna.digitize(allocator, "");
    defer allocator.free(dsq);

    try dna.reverseComplement(dsq);
    try std.testing.expectEqual(@as(usize, 0), dsq.len);
}

test "reverseComplement amino acid returns error" {
    const allocator = std.testing.allocator;
    var dsq = try amino.digitize(allocator, "ACD");
    defer allocator.free(dsq);

    try std.testing.expectError(error.NoComplement, amino.reverseComplement(dsq));
}

test "reverseComplement degeneracies" {
    const allocator = std.testing.allocator;
    // R(AG) -> complement Y(CT), and vice versa
    var dsq = try dna.digitize(allocator, "RY");
    defer allocator.free(dsq);

    try dna.reverseComplement(dsq);
    const text = try dna.textize(allocator, dsq);
    defer allocator.free(text);

    try std.testing.expectEqualStrings("RY", text);
}
```

- [ ] **Step 2: Run tests**

Run: `cd /Users/nagaet/zeacel && zig build test`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
cd /Users/nagaet/zeacel
git add src/alphabet.zig
git commit -m "feat: add digitize, textize, and reverseComplement"
```

---

### Task 6: Alphabet Guessing from Sequence Content

**Files:**
- Modify: `src/alphabet.zig`

Easel's `esl_abc_GuessAlphabet` counts character frequencies to decide if a sequence is DNA/RNA or amino acid. This is useful for auto-detecting file content.

- [ ] **Step 1: Write test and implementation for guessType**

Add to `src/alphabet.zig` inside the `Alphabet` struct definition (as a standalone function, since it does not operate on an existing alphabet):

```zig
/// Guess the alphabet type from sequence content.
/// Counts character frequencies: if >80% of residues are in {A,C,G,T,U,N},
/// guess nucleic acid (DNA). Otherwise guess amino acid.
/// Returns null if the sequence is empty or contains no valid residues.
pub fn guessType(text: []const u8) ?AlphabetType {
    var nuc_count: usize = 0;
    var amino_only_count: usize = 0;
    var total: usize = 0;

    for (text) |c| {
        const upper = if (c >= 'a' and c <= 'z') c - 32 else c;
        switch (upper) {
            'A', 'C', 'G', 'T', 'U', 'N' => {
                nuc_count += 1;
                total += 1;
            },
            'D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'P', 'Q', 'R', 'S', 'V', 'W', 'Y' => {
                // These are valid in both amino and nucleic alphabets (some as degens),
                // but D, E, F, I, L, P, Q are amino-only canonical residues.
                switch (upper) {
                    'E', 'F', 'I', 'L', 'P', 'Q' => amino_only_count += 1,
                    else => {},
                }
                total += 1;
            },
            '-', '.', '_', '*', '~' => {}, // skip gap/special chars
            else => {},
        }
    }

    if (total == 0) return null;

    // If any amino-only canonical residues found, it's protein
    if (amino_only_count > 0) return .amino;

    // If > 80% nucleic characters, guess DNA
    if (nuc_count * 100 / total > 80) return .dna;

    return .amino;
}
```

Add tests:

```zig
test "guessType pure DNA" {
    try std.testing.expectEqual(AlphabetType.dna, Alphabet.guessType("ACGTACGTACGT").?);
}

test "guessType DNA with some degeneracies" {
    try std.testing.expectEqual(AlphabetType.dna, Alphabet.guessType("ACGTNNNNACGT").?);
}

test "guessType amino acid" {
    try std.testing.expectEqual(AlphabetType.amino, Alphabet.guessType("ACDEFGHIKLMNPQRSTVWY").?);
}

test "guessType amino acid detected by unique residues" {
    // Contains E, F, L which are amino-only
    try std.testing.expectEqual(AlphabetType.amino, Alphabet.guessType("ACGTEF").?);
}

test "guessType empty returns null" {
    try std.testing.expect(Alphabet.guessType("") == null);
}

test "guessType case-insensitive" {
    try std.testing.expectEqual(AlphabetType.dna, Alphabet.guessType("acgtacgt").?);
}

test "guessType gaps only returns null" {
    try std.testing.expect(Alphabet.guessType("---...") == null);
}
```

- [ ] **Step 2: Run tests**

Run: `cd /Users/nagaet/zeacel && zig build test`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
cd /Users/nagaet/zeacel
git add src/alphabet.zig
git commit -m "feat: add alphabet type guessing from sequence content"
```

---

### Task 7: Final Cleanup and Verification

**Files:**
- Modify: `src/root.zig` (if needed)

- [ ] **Step 1: Run full test suite**

Run: `cd /Users/nagaet/zeacel && zig build test 2>&1`
Expected: All tests pass with no warnings.

- [ ] **Step 2: Verify library builds as static lib**

Run: `cd /Users/nagaet/zeacel && zig build`
Expected: Successful build, artifact in `zig-out/lib/`.

- [ ] **Step 3: Verify clean git status**

Run: `cd /Users/nagaet/zeacel && git log --oneline && echo "---" && git status`
Expected: Clean working tree with commits:
1. `chore: initialize zeacel Zig project`
2. `feat: add Alphabet core types with classification methods`
3. `feat: add comptime DNA alphabet with encode/decode/complement`
4. `feat: add RNA and amino acid alphabets`
5. `feat: add digitize, textize, and reverseComplement`
6. `feat: add alphabet type guessing from sequence content`

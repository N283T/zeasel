// Biological sequence alphabets (DNA, RNA, Amino Acid).

const std = @import("std");

pub const AlphabetType = enum {
    dna,
    rna,
    amino,
};

pub const INVALID_CODE: u8 = 0xFF;

pub const Alphabet = struct {
    kind: AlphabetType,
    k: u8, // Canonical residue count (4 for DNA/RNA, 20 for amino)
    kp: u8, // Total symbol count
    symbols: []const u8,
    encode_map: [128]u8, // ASCII -> digital code. INVALID_CODE for unmapped.
    complement: ?[]const u8, // Complement table (DNA/RNA only)
    // Degeneracy bitmasks: for code x, degen[x] has bit y set if canonical residue y
    // is part of the degeneracy. Canonical codes map to themselves. Max 32 symbols.
    degen: [32]u32,
    ndegen: [32]u8, // Number of canonical residues per code

    /// True if code represents a residue (canonical, degenerate, or unknown).
    pub fn isResidue(self: *const Alphabet, code: u8) bool {
        return self.isCanonical(code) or self.isDegenerate(code) or self.isUnknown(code);
    }

    pub fn encode(self: *const Alphabet, char: u8) error{InvalidCharacter}!u8 {
        if (char > 127) return error.InvalidCharacter;
        const code = self.encode_map[char];
        if (code == INVALID_CODE) return error.InvalidCharacter;
        return code;
    }

    pub fn decode(self: *const Alphabet, code: u8) u8 {
        std.debug.assert(code < self.symbols.len);
        return self.symbols[code];
    }

    pub fn isCanonical(self: *const Alphabet, code: u8) bool {
        return code < self.k;
    }

    pub fn isGap(self: *const Alphabet, code: u8) bool {
        return code == self.k;
    }

    pub fn isDegenerate(self: *const Alphabet, code: u8) bool {
        return code > self.k and code < self.kp - 3;
    }

    pub fn isUnknown(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 3;
    }

    pub fn isNonresidue(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 2;
    }

    pub fn isMissing(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 1;
    }

    pub fn gapCode(self: *const Alphabet) u8 {
        return self.k;
    }

    pub fn unknownCode(self: *const Alphabet) u8 {
        return self.kp - 3;
    }

    pub fn digitize(self: *const Alphabet, allocator: std.mem.Allocator, text: []const u8) ![]u8 {
        const dsq = try allocator.alloc(u8, text.len);
        errdefer allocator.free(dsq);
        for (text, 0..) |char, i| {
            dsq[i] = self.encode(char) catch return error.InvalidCharacter;
        }
        return dsq;
    }

    pub fn textize(self: *const Alphabet, allocator: std.mem.Allocator, dsq: []const u8) ![]u8 {
        const text = try allocator.alloc(u8, dsq.len);
        for (dsq, 0..) |code, i| {
            text[i] = self.decode(code);
        }
        return text;
    }

    pub fn reverseComplement(self: *const Alphabet, dsq: []u8) error{NoComplement}!void {
        const comp = self.complement orelse return error.NoComplement;
        var i: usize = 0;
        var j: usize = dsq.len;
        while (i < j) {
            j -= 1;
            const tmp = comp[dsq[i]];
            dsq[i] = comp[dsq[j]];
            dsq[j] = tmp;
            i += 1;
        }
    }

    /// Average score for a (possibly degenerate) residue code, given a
    /// score vector sc[0..K-1] for canonical residues.
    /// For canonical residues, returns sc[x] directly.
    /// For degenerate codes, returns the mean of scores for constituent residues.
    /// For gap/nonresidue/missing, returns 0.
    pub fn avgScore(self: *const Alphabet, x: u8, sc: []const f64) f64 {
        if (x >= self.kp) return 0;
        if (self.isCanonical(x)) return sc[x];
        if (!self.isResidue(x)) return 0;
        const mask = self.degen[x];
        const n = self.ndegen[x];
        if (n == 0) return 0;
        var sum: f64 = 0;
        for (0..self.k) |i| {
            if (mask & (@as(u32, 1) << @intCast(i)) != 0) {
                sum += sc[i];
            }
        }
        return sum / @as(f64, @floatFromInt(n));
    }

    /// Fill out a score vector for degenerate codes using avgScore.
    /// sc must be length kp. Canonical scores (0..K-1) must be pre-filled.
    /// This fills in degenerate codes (K+1..Kp-3).
    pub fn avgScVec(self: *const Alphabet, sc: []f64) void {
        var x: u8 = self.k + 1;
        while (x <= self.kp -| 3) : (x += 1) {
            if (x < self.kp) {
                sc[x] = self.avgScore(x, sc[0..self.k]);
            }
        }
    }

    /// Expected score for a (possibly degenerate) residue code, given
    /// a score vector sc[0..K-1] and background probabilities p[0..K-1].
    /// For canonical residues, returns sc[x].
    /// For degenerate codes, returns weighted average: sum(sc[i]*p[i]) / sum(p[i])
    /// for constituent residues i.
    pub fn expectScore(self: *const Alphabet, x: u8, sc: []const f64, p: []const f64) f64 {
        if (x >= self.kp) return 0;
        if (self.isCanonical(x)) return sc[x];
        if (!self.isResidue(x)) return 0;
        const mask = self.degen[x];
        var result: f64 = 0;
        var denom: f64 = 0;
        for (0..self.k) |i| {
            if (mask & (@as(u32, 1) << @intCast(i)) != 0) {
                result += sc[i] * p[i];
                denom += p[i];
            }
        }
        if (denom == 0) return 0;
        return result / denom;
    }

    /// Fill out a score vector for degenerate codes using expectScore.
    pub fn expectScVec(self: *const Alphabet, sc: []f64, p: []const f64) void {
        var x: u8 = self.k + 1;
        while (x <= self.kp -| 3) : (x += 1) {
            if (x < self.kp) {
                sc[x] = self.expectScore(x, sc[0..self.k], p);
            }
        }
    }

    /// Count a (possibly degenerate) symbol into a count vector ct[0..K-1].
    /// Canonical residues add wt to ct[x]. Degenerate codes distribute wt
    /// equally across constituent canonical residues. Gap adds to ct[K] if
    /// ct is long enough. Missing/nonresidue are ignored.
    pub fn count(self: *const Alphabet, ct: []f64, x: u8, wt: f64) void {
        if (x >= self.kp) return;
        if (self.isCanonical(x)) {
            ct[x] += wt;
            return;
        }
        if (self.isGap(x)) {
            if (ct.len > self.k) ct[self.k] += wt;
            return;
        }
        if (self.isMissing(x) or self.isNonresidue(x)) return;
        // Degenerate or unknown: distribute weight equally
        const mask = self.degen[x];
        const n = self.ndegen[x];
        if (n == 0) return;
        const share = wt / @as(f64, @floatFromInt(n));
        for (0..self.k) |i| {
            if (mask & (@as(u32, 1) << @intCast(i)) != 0) {
                ct[i] += share;
            }
        }
    }

    /// Match probability between two residue codes, accounting for degeneracy.
    /// If p is null, uses uniform (1/K) probabilities.
    /// Returns 0 for non-residue comparisons (gap, missing, etc).
    pub fn match(self: *const Alphabet, x: u8, y: u8, p: ?[]const f64) f64 {
        // Easy case: both canonical
        if (self.isCanonical(x) and self.isCanonical(y)) {
            return if (x == y) 1.0 else 0.0;
        }
        // Non-residue: return 0
        if (!self.isResidue(x) or !self.isResidue(y)) return 0;

        const mx = self.degen[x];
        const my = self.degen[y];

        var prob: f64 = 0;
        var sx: f64 = 0;
        var sy: f64 = 0;

        for (0..self.k) |i| {
            const bit = @as(u32, 1) << @intCast(i);
            const in_x = (mx & bit) != 0;
            const in_y = (my & bit) != 0;
            const pi: f64 = if (p) |bg| bg[i] else 1.0 / @as(f64, @floatFromInt(self.k));
            if (in_x) sx += pi;
            if (in_y) sy += pi;
            if (in_x and in_y) prob += pi * pi;
        }
        if (sx == 0 or sy == 0) return 0;
        return prob / (sx * sy);
    }

    /// Count non-gap residues in a digital sequence.
    pub fn dsqrlen(self: *const Alphabet, dsq: []const u8) usize {
        var n: usize = 0;
        for (dsq) |code| {
            if (self.isResidue(code)) n += 1;
        }
        return n;
    }

    /// Remove gap characters from an aligned digital sequence.
    /// Returns a new slice with gaps removed. Caller owns the returned memory.
    pub fn dealign(self: *const Alphabet, allocator: std.mem.Allocator, dsq: []const u8) ![]u8 {
        const rlen = self.dsqrlen(dsq);
        const result = try allocator.alloc(u8, rlen);
        var j: usize = 0;
        for (dsq) |code| {
            if (self.isResidue(code)) {
                result[j] = code;
                j += 1;
            }
        }
        return result;
    }

    /// Validate that a text sequence contains only valid characters for this alphabet.
    pub fn validateSeq(self: *const Alphabet, text: []const u8) bool {
        for (text) |ch| {
            if (ch > 127) return false;
            if (self.encode_map[ch] == INVALID_CODE) return false;
        }
        return true;
    }
};

// --- Task 3: Comptime DNA alphabet construction ---

/// Build an encode map at comptime from a symbols string and synonym pairs.
///
/// - Each char in `symbols` maps to its index.
/// - Lowercase versions of each symbol map to the same code.
/// - Each synonym pair `[2]u8{src, dst}` means encode_map[src] = encode_map[dst].
///   Lowercase of src is also mapped to the same code.
fn buildEncodeMap(
    comptime symbols: []const u8,
    comptime synonyms: []const [2]u8,
) [128]u8 {
    var map = [_]u8{INVALID_CODE} ** 128;

    // Map each symbol and its lowercase to its index.
    for (symbols, 0..) |sym, i| {
        map[sym] = @intCast(i);
        const lower = std.ascii.toLower(sym);
        if (lower != sym) {
            map[lower] = @intCast(i);
        }
    }

    // Apply synonym pairs: encode_map[src] = encode_map[dst]
    // Also map lowercase of src.
    for (synonyms) |pair| {
        const src = pair[0];
        const dst = pair[1];
        map[src] = map[dst];
        const lower_src = std.ascii.toLower(src);
        if (lower_src != src) {
            map[lower_src] = map[dst];
        }
    }

    return map;
}

/// Build degeneracy bitmask and count tables at comptime.
/// k = number of canonical residues, kp = total symbol count.
/// entries is an array of {code, bitmask} pairs for degenerate codes.
fn buildDegenTable(
    comptime k: u8,
    comptime kp: u8,
    comptime entries: []const [2]u32,
) struct { masks: [32]u32, counts: [32]u8 } {
    var masks = [_]u32{0} ** 32;
    var counts = [_]u8{0} ** 32;

    // Canonical residues map to themselves
    for (0..k) |i| {
        masks[i] = @as(u32, 1) << @intCast(i);
        counts[i] = 1;
    }
    // Gap, nonresidue, missing: no residues
    _ = kp;

    // Fill in degenerate entries
    for (entries) |entry| {
        const code: usize = entry[0];
        const mask: u32 = entry[1];
        masks[code] = mask;
        counts[code] = @popCount(mask);
    }

    return .{ .masks = masks, .counts = counts };
}

// DNA symbol string. Positions:
//   0..3  : canonical residues A, C, G, T
//   4     : gap (-)
//   5..14 : degeneracy symbols R, Y, M, K, S, W, H, B, V, D
//   15    : unknown (N)
//   16    : nonresidue (*)
//   17    : missing (~)
pub const dna_symbols: []const u8 = "ACGT-RYMKSWHBVDN*~";

// Complement table: indexed by digital code, value is complement code.
// A(0)<->T(3), C(1)<->G(2), -(4)->-(4),
// R(5)<->Y(6), M(7)<->K(8), S(9)->S(9), W(10)->W(10),
// H(11)<->D(14), B(12)<->V(13),
// N(15)->N(15), *(16)->*(16), ~(17)->~(17)
pub const dna_complement_table: [18]u8 = .{
    3, 2, 1, 0, // A->T, C->G, G->C, T->A
    4, // gap->gap
    6, 5, // R->Y, Y->R
    8, 7, // M->K, K->M
    9, 10, // S->S, W->W
    14, 13, 12, 11, // H->D, B->V, V->B, D->H
    15, 16, 17, // N->N, *->*, ~->~
};

// Synonym pairs for DNA: U->T, X->N, I->A, _->-, .->-
const dna_synonyms: []const [2]u8 = &.{
    .{ 'U', 'T' },
    .{ 'X', 'N' },
    .{ 'I', 'A' },
    .{ '_', '-' },
    .{ '.', '-' },
};

// DNA degeneracy bitmasks (bit 0=A, 1=C, 2=G, 3=T):
//   R=A|G, Y=C|T, M=A|C, K=G|T, S=C|G, W=A|T,
//   H=A|C|T, B=C|G|T, V=A|C|G, D=A|G|T, N=A|C|G|T
const dna_degen = buildDegenTable(4, 18, &.{
    .{ 5,  0b0101 }, // R = A|G
    .{ 6,  0b1010 }, // Y = C|T
    .{ 7,  0b0011 }, // M = A|C
    .{ 8,  0b1100 }, // K = G|T
    .{ 9,  0b0110 }, // S = C|G
    .{ 10, 0b1001 }, // W = A|T
    .{ 11, 0b1011 }, // H = A|C|T
    .{ 12, 0b1110 }, // B = C|G|T
    .{ 13, 0b0111 }, // V = A|C|G
    .{ 14, 0b1101 }, // D = A|G|T
    .{ 15, 0b1111 }, // N = any
});

pub const dna: Alphabet = .{
    .kind = .dna,
    .k = 4,
    .kp = 18,
    .symbols = dna_symbols,
    .encode_map = buildEncodeMap(dna_symbols, dna_synonyms),
    .complement = &dna_complement_table,
    .degen = dna_degen.masks,
    .ndegen = dna_degen.counts,
};

// RNA symbol string. Same structure as DNA, U instead of T.
//   0..3  : canonical residues A, C, G, U
//   4     : gap (-)
//   5..14 : degeneracy symbols R, Y, M, K, S, W, H, B, V, D
//   15    : unknown (N)
//   16    : nonresidue (*)
//   17    : missing (~)
pub const rna_symbols: []const u8 = "ACGU-RYMKSWHBVDN*~";

// Complement table for RNA: A(0)<->U(3), C(1)<->G(2), same degeneracy pattern as DNA.
pub const rna_complement_table: [18]u8 = .{
    3, 2, 1, 0, // A->U, C->G, G->C, U->A
    4, // gap->gap
    6, 5, // R->Y, Y->R
    8, 7, // M->K, K->M
    9, 10, // S->S, W->W
    14, 13, 12, 11, // H->D, B->V, V->B, D->H
    15, 16, 17, // N->N, *->*, ~->~
};

// Synonym pairs for RNA: T->U, X->N, I->A, _->-, .->-
const rna_synonyms: []const [2]u8 = &.{
    .{ 'T', 'U' },
    .{ 'X', 'N' },
    .{ 'I', 'A' },
    .{ '_', '-' },
    .{ '.', '-' },
};

// RNA uses the same degeneracy pattern as DNA
const rna_degen = buildDegenTable(4, 18, &.{
    .{ 5,  0b0101 }, // R = A|G
    .{ 6,  0b1010 }, // Y = C|U
    .{ 7,  0b0011 }, // M = A|C
    .{ 8,  0b1100 }, // K = G|U
    .{ 9,  0b0110 }, // S = C|G
    .{ 10, 0b1001 }, // W = A|U
    .{ 11, 0b1011 }, // H = A|C|U
    .{ 12, 0b1110 }, // B = C|G|U
    .{ 13, 0b0111 }, // V = A|C|G
    .{ 14, 0b1101 }, // D = A|G|U
    .{ 15, 0b1111 }, // N = any
});

pub const rna: Alphabet = .{
    .kind = .rna,
    .k = 4,
    .kp = 18,
    .symbols = rna_symbols,
    .encode_map = buildEncodeMap(rna_symbols, rna_synonyms),
    .complement = &rna_complement_table,
    .degen = rna_degen.masks,
    .ndegen = rna_degen.counts,
};

// Amino acid symbol string. Positions:
//   0..19 : canonical residues A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
//   20    : gap (-)
//   21..25: degeneracy symbols B, J, Z, O, U
//   26    : unknown (X)
//   27    : nonresidue (*)
//   28    : missing (~)
pub const amino_symbols: []const u8 = "ACDEFGHIKLMNPQRSTVWY-BJZOUX*~";

// No complement table for amino acids.
// Synonym pairs: _->-, .->-
const amino_synonyms: []const [2]u8 = &.{
    .{ '_', '-' },
    .{ '.', '-' },
};

// Amino acid degeneracies:
//   B(21) = D|N, J(22) = I|L, Z(23) = E|Q,
//   O(24) = pyrrolysine (maps to K), U(25) = selenocysteine (maps to C),
//   X(26) = any (all 20)
// Alphabet order: ACDEFGHIKLMNPQRSTVWY (indices 0-19)
//   C=1, K=9
const amino_degen = buildDegenTable(20, 29, &.{
    .{ 21, (1 << 2) | (1 << 11) },        // B = D|N
    .{ 22, (1 << 7) | (1 << 9) },         // J = I|L
    .{ 23, (1 << 3) | (1 << 13) },        // Z = E|Q
    .{ 24, 1 << 9 },                      // O = K (pyrrolysine)
    .{ 25, 1 << 1 },                      // U = C (selenocysteine)
    .{ 26, (1 << 20) - 1 },               // X = all 20
});

pub const amino: Alphabet = .{
    .kind = .amino,
    .k = 20,
    .kp = 29,
    .symbols = amino_symbols,
    .encode_map = buildEncodeMap(amino_symbols, amino_synonyms),
    .complement = null,
    .degen = amino_degen.masks,
    .ndegen = amino_degen.counts,
};

// --- Task 6: Guess alphabet type from sequence content ---

/// Guess the alphabet type from sequence content.
/// If any amino-only canonical residues (E, F, I, L, P, Q) are found, returns .amino.
/// If >80% of recognized residues are nucleic (A,C,G,T,U,N), returns .dna.
/// Otherwise returns .amino.
/// Returns null if the sequence is empty or has no recognized residues.
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
                switch (upper) {
                    'E', 'F', 'I', 'L', 'P', 'Q' => amino_only_count += 1,
                    else => {},
                }
                total += 1;
            },
            '-', '.', '_', '*', '~' => {},
            else => {},
        }
    }

    if (total == 0) return null;
    if (amino_only_count > 0) return .amino;
    if (nuc_count * 100 / total > 80) return .dna;
    return .amino;
}

// --- Tests for Task 2 ---

test "Alphabet: struct size constants for DNA" {
    const alpha = dna;

    try std.testing.expectEqual(@as(u8, 4), alpha.k);
    try std.testing.expectEqual(@as(u8, 18), alpha.kp);
    try std.testing.expectEqual(@as(u8, 4), alpha.gapCode()); // gap at index k=4
    try std.testing.expectEqual(@as(u8, 15), alpha.unknownCode()); // kp-3 = 18-3 = 15
}

test "Alphabet: classification methods" {
    const alpha = dna;

    // Canonical: codes 0..3 (A, C, G, T)
    try std.testing.expect(alpha.isCanonical(0)); // A
    try std.testing.expect(alpha.isCanonical(3)); // T
    try std.testing.expect(!alpha.isCanonical(4)); // gap

    // Gap: code == k == 4
    try std.testing.expect(alpha.isGap(4));
    try std.testing.expect(!alpha.isGap(0));
    try std.testing.expect(!alpha.isGap(5));

    // Degenerate: codes > k and < kp-2, i.e., 5..14 (R,Y,M,K,S,W,H,B,V,D)
    try std.testing.expect(alpha.isDegenerate(5)); // R
    try std.testing.expect(alpha.isDegenerate(14)); // D
    try std.testing.expect(!alpha.isDegenerate(4)); // gap
    try std.testing.expect(!alpha.isDegenerate(15)); // N (unknown)
    try std.testing.expect(!alpha.isDegenerate(16)); // * (nonresidue)

    // Unknown: code == kp-3 == 15 (N)
    try std.testing.expect(alpha.isUnknown(15)); // N
    try std.testing.expect(!alpha.isUnknown(14));
    try std.testing.expect(!alpha.isUnknown(16));

    // Nonresidue: code == kp-2 == 16 (*)
    try std.testing.expect(alpha.isNonresidue(16)); // *
    try std.testing.expect(!alpha.isNonresidue(15));
    try std.testing.expect(!alpha.isNonresidue(17));

    // Missing: code == kp-1 == 17 (~)
    try std.testing.expect(alpha.isMissing(17)); // ~
    try std.testing.expect(!alpha.isMissing(16));
    try std.testing.expect(!alpha.isMissing(0));
}

// --- Tests for Task 3 ---

test "DNA encode canonical residues" {
    try std.testing.expectEqual(@as(u8, 0), try dna.encode('A'));
    try std.testing.expectEqual(@as(u8, 1), try dna.encode('C'));
    try std.testing.expectEqual(@as(u8, 2), try dna.encode('G'));
    try std.testing.expectEqual(@as(u8, 3), try dna.encode('T'));
}

test "DNA encode case-insensitive" {
    try std.testing.expectEqual(@as(u8, 0), try dna.encode('a'));
    try std.testing.expectEqual(@as(u8, 1), try dna.encode('c'));
    try std.testing.expectEqual(@as(u8, 2), try dna.encode('g'));
    try std.testing.expectEqual(@as(u8, 3), try dna.encode('t'));
}

test "DNA encode synonyms" {
    // U -> T (3)
    try std.testing.expectEqual(@as(u8, 3), try dna.encode('U'));
    try std.testing.expectEqual(@as(u8, 3), try dna.encode('u'));
    // X -> N (15)
    try std.testing.expectEqual(@as(u8, 15), try dna.encode('X'));
    try std.testing.expectEqual(@as(u8, 15), try dna.encode('x'));
    // _ -> - (4)
    try std.testing.expectEqual(@as(u8, 4), try dna.encode('_'));
    // . -> - (4)
    try std.testing.expectEqual(@as(u8, 4), try dna.encode('.'));
}

test "DNA encode degeneracies" {
    try std.testing.expectEqual(@as(u8, 5), try dna.encode('R'));
    try std.testing.expectEqual(@as(u8, 6), try dna.encode('Y'));
    try std.testing.expectEqual(@as(u8, 15), try dna.encode('N'));
}

test "DNA encode invalid characters" {
    try std.testing.expectError(error.InvalidCharacter, dna.encode('!'));
    try std.testing.expectError(error.InvalidCharacter, dna.encode('Z'));
    try std.testing.expectError(error.InvalidCharacter, dna.encode(0));
}

test "DNA decode" {
    try std.testing.expectEqual(@as(u8, 'A'), dna.decode(0));
    try std.testing.expectEqual(@as(u8, 'C'), dna.decode(1));
    try std.testing.expectEqual(@as(u8, 'G'), dna.decode(2));
    try std.testing.expectEqual(@as(u8, 'T'), dna.decode(3));
    try std.testing.expectEqual(@as(u8, '-'), dna.decode(4));
    try std.testing.expectEqual(@as(u8, 'N'), dna.decode(15));
}

test "DNA complement table" {
    // A(0) <-> T(3)
    try std.testing.expectEqual(@as(u8, 3), dna_complement_table[0]);
    try std.testing.expectEqual(@as(u8, 0), dna_complement_table[3]);
    // C(1) <-> G(2)
    try std.testing.expectEqual(@as(u8, 2), dna_complement_table[1]);
    try std.testing.expectEqual(@as(u8, 1), dna_complement_table[2]);
    // R(5) <-> Y(6)
    try std.testing.expectEqual(@as(u8, 6), dna_complement_table[5]);
    try std.testing.expectEqual(@as(u8, 5), dna_complement_table[6]);
}

// --- Tests for Task 4: RNA and Amino Acid Alphabets ---

test "RNA encode canonical residues" {
    try std.testing.expectEqual(@as(u8, 0), try rna.encode('A'));
    try std.testing.expectEqual(@as(u8, 3), try rna.encode('U'));
}

test "RNA encode synonym T->U" {
    try std.testing.expectEqual(@as(u8, 3), try rna.encode('T'));
    try std.testing.expectEqual(@as(u8, 3), try rna.encode('t'));
}

test "RNA encode case-insensitive" {
    try std.testing.expectEqual(@as(u8, 0), try rna.encode('a'));
}

test "RNA complement" {
    const comp = rna.complement.?;
    // A(0) <-> U(3)
    try std.testing.expectEqual(@as(u8, 3), comp[0]);
    try std.testing.expectEqual(@as(u8, 0), comp[3]);
}

test "Amino encode canonical residues" {
    try std.testing.expectEqual(@as(u8, 0), try amino.encode('A'));
    try std.testing.expectEqual(@as(u8, 19), try amino.encode('Y'));
}

test "Amino encode case-insensitive" {
    try std.testing.expectEqual(@as(u8, 0), try amino.encode('a'));
}

test "Amino gap codes" {
    try std.testing.expectEqual(@as(u8, 20), try amino.encode('-'));
    try std.testing.expectEqual(@as(u8, 20), try amino.encode('_'));
    try std.testing.expectEqual(@as(u8, 20), try amino.encode('.'));
}

test "Amino sizes" {
    try std.testing.expectEqual(@as(u8, 20), amino.k);
    try std.testing.expectEqual(@as(u8, 29), amino.kp);
}

test "Amino unknown X at kp-3" {
    const x_code = amino.kp - 3; // 29 - 3 = 26
    try std.testing.expectEqual(@as(u8, 26), try amino.encode('X'));
    try std.testing.expectEqual(@as(u8, 'X'), amino.decode(26));
    try std.testing.expectEqual(x_code, try amino.encode('X'));
}

test "Amino no complement" {
    try std.testing.expectEqual(@as(?[]const u8, null), amino.complement);
}

test "Amino degeneracy codes" {
    try std.testing.expectEqual(@as(u8, 21), try amino.encode('B'));
    try std.testing.expectEqual(@as(u8, 22), try amino.encode('J'));
    try std.testing.expectEqual(@as(u8, 23), try amino.encode('Z'));
    try std.testing.expectEqual(@as(u8, 24), try amino.encode('O'));
    try std.testing.expectEqual(@as(u8, 25), try amino.encode('U'));
}

// --- Tests for Task 5: digitize, textize, reverseComplement ---

test "digitize ACGT" {
    const result = try dna.digitize(std.testing.allocator, "ACGT");
    defer std.testing.allocator.free(result);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, result);
}

test "digitize lowercase acgt" {
    const result = try dna.digitize(std.testing.allocator, "acgt");
    defer std.testing.allocator.free(result);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, result);
}

test "digitize with degeneracies ACNR" {
    const result = try dna.digitize(std.testing.allocator, "ACNR");
    defer std.testing.allocator.free(result);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 15, 5 }, result);
}

test "digitize invalid character returns error" {
    try std.testing.expectError(error.InvalidCharacter, dna.digitize(std.testing.allocator, "ACGZ"));
}

test "digitize empty string" {
    const result = try dna.digitize(std.testing.allocator, "");
    defer std.testing.allocator.free(result);
    try std.testing.expectEqual(@as(usize, 0), result.len);
}

test "textize round-trip DNA" {
    const input = "ACGTRYNSWHBVD";
    const dsq = try dna.digitize(std.testing.allocator, input);
    defer std.testing.allocator.free(dsq);
    const text = try dna.textize(std.testing.allocator, dsq);
    defer std.testing.allocator.free(text);
    try std.testing.expectEqualSlices(u8, input, text);
}

test "textize round-trip amino acid" {
    const input = "ACDEFGHIKLMNPQRSTVWY";
    const dsq = try amino.digitize(std.testing.allocator, input);
    defer std.testing.allocator.free(dsq);
    const text = try amino.textize(std.testing.allocator, dsq);
    defer std.testing.allocator.free(text);
    try std.testing.expectEqualSlices(u8, input, text);
}

test "reverseComplement palindrome ACGT" {
    var seq = [_]u8{ 0, 1, 2, 3 }; // A, C, G, T
    try dna.reverseComplement(&seq);
    // reverse of [T,G,C,A] complement = [A,C,G,T] -> palindrome
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, &seq);
}

test "reverseComplement AACG -> CGTT" {
    var seq = [_]u8{ 0, 0, 1, 2 }; // A, A, C, G
    try dna.reverseComplement(&seq);
    // reverse of [G,C,A,A] complement = [C,G,T,T]
    try std.testing.expectEqualSlices(u8, &[_]u8{ 1, 2, 3, 3 }, &seq);
}

test "reverseComplement odd length ACG -> CGT" {
    var seq = [_]u8{ 0, 1, 2 }; // A, C, G
    try dna.reverseComplement(&seq);
    // reverse = [G,C,A], complement = [C,G,T]
    try std.testing.expectEqualSlices(u8, &[_]u8{ 1, 2, 3 }, &seq);
}

test "reverseComplement single A -> T" {
    var seq = [_]u8{0}; // A
    try dna.reverseComplement(&seq);
    try std.testing.expectEqualSlices(u8, &[_]u8{3}, &seq);
}

test "reverseComplement empty" {
    var seq = [_]u8{};
    try dna.reverseComplement(&seq);
    try std.testing.expectEqual(@as(usize, 0), seq.len);
}

test "reverseComplement amino returns error" {
    var seq = [_]u8{ 0, 1, 2 };
    try std.testing.expectError(error.NoComplement, amino.reverseComplement(&seq));
}

test "reverseComplement degeneracies RY -> RY" {
    // R=5, Y=6. reverse of [Y,R] complement = [R,Y] -> same
    var seq = [_]u8{ 5, 6 }; // R, Y
    try dna.reverseComplement(&seq);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 5, 6 }, &seq);
}

// --- Tests for Task 6: guessType ---

test "guessType: pure DNA returns .dna" {
    try std.testing.expectEqual(@as(?AlphabetType, .dna), guessType("ACGTACGTACGT"));
}

test "guessType: DNA with N's returns .dna" {
    try std.testing.expectEqual(@as(?AlphabetType, .dna), guessType("ACGTNNNNACGT"));
}

test "guessType: all canonical amino acids returns .amino" {
    try std.testing.expectEqual(@as(?AlphabetType, .amino), guessType("ACDEFGHIKLMNPQRSTVWY"));
}

test "guessType: DNA mixed with amino-only residues returns .amino" {
    try std.testing.expectEqual(@as(?AlphabetType, .amino), guessType("ACGTEF"));
}

test "guessType: empty string returns null" {
    try std.testing.expectEqual(@as(?AlphabetType, null), guessType(""));
}

test "guessType: lowercase DNA returns .dna" {
    try std.testing.expectEqual(@as(?AlphabetType, .dna), guessType("acgtacgt"));
}

test "guessType: gaps only returns null" {
    try std.testing.expectEqual(@as(?AlphabetType, null), guessType("---..."));
}

// --- Tests for degeneracy functions ---

test "isResidue: canonical and degenerate are residues" {
    try std.testing.expect(dna.isResidue(0)); // A
    try std.testing.expect(dna.isResidue(5)); // R (degenerate)
    try std.testing.expect(dna.isResidue(15)); // N (unknown)
    try std.testing.expect(!dna.isResidue(4)); // gap
    try std.testing.expect(!dna.isResidue(16)); // nonresidue
    try std.testing.expect(!dna.isResidue(17)); // missing
}

test "avgScore: canonical returns direct score" {
    const sc = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), dna.avgScore(0, &sc), 1e-9); // A
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), dna.avgScore(3, &sc), 1e-9); // T
}

test "avgScore: R (A|G) returns average of A and G scores" {
    const sc = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    // R = A|G, avg = (1+3)/2 = 2.0
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), dna.avgScore(5, &sc), 1e-9);
}

test "avgScore: N (all 4) returns mean of all" {
    const sc = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    // N = all 4, avg = (1+2+3+4)/4 = 2.5
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), dna.avgScore(15, &sc), 1e-9);
}

test "avgScore: gap returns 0" {
    const sc = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), dna.avgScore(4, &sc), 1e-9);
}

test "expectScore: canonical returns direct score" {
    const sc = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    const p = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), dna.expectScore(0, &sc, &p), 1e-9);
}

test "expectScore: R with non-uniform background" {
    const sc = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    const p = [_]f64{ 0.3, 0.2, 0.3, 0.2 };
    // R = A|G, expected = (1.0*0.3 + 3.0*0.3) / (0.3+0.3) = 1.2/0.6 = 2.0
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), dna.expectScore(5, &sc, &p), 1e-9);
}

test "count: canonical residue" {
    var ct = [_]f64{ 0, 0, 0, 0, 0 };
    dna.count(&ct, 0, 1.0); // A
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), ct[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), ct[1], 1e-9);
}

test "count: degenerate R distributes equally to A and G" {
    var ct = [_]f64{ 0, 0, 0, 0, 0 };
    dna.count(&ct, 5, 1.0); // R = A|G
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), ct[0], 1e-9); // A
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), ct[1], 1e-9); // C
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), ct[2], 1e-9); // G
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), ct[3], 1e-9); // T
}

test "count: N distributes to all 4" {
    var ct = [_]f64{ 0, 0, 0, 0, 0 };
    dna.count(&ct, 15, 2.0); // N = all, wt=2
    for (0..4) |i| {
        try std.testing.expectApproxEqAbs(@as(f64, 0.5), ct[i], 1e-9);
    }
}

test "count: gap counts into ct[K]" {
    var ct = [_]f64{ 0, 0, 0, 0, 0 };
    dna.count(&ct, 4, 1.0); // gap
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), ct[4], 1e-9);
}

test "match: identical canonical = 1.0" {
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), dna.match(0, 0, null), 1e-9);
}

test "match: different canonical = 0.0" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), dna.match(0, 1, null), 1e-9);
}

test "match: R vs A = 0.5 (uniform)" {
    // R = A|G. A matches R with prob 0.5
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), dna.match(5, 0, null), 1e-9);
}

test "match: N vs R = 0.5 (uniform)" {
    // N = ACGT (4), R = AG (2). Overlap = AG (2).
    // prob = sum(p_i^2 for i in overlap) / (sx * sy)
    // = (0.25^2 + 0.25^2) / (1.0 * 0.5) = 0.125 / 0.5 = 0.25
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), dna.match(15, 5, null), 1e-9);
}

test "dsqrlen: counts non-gap residues" {
    const dsq = [_]u8{ 0, 4, 1, 2, 4, 3 }; // A, gap, C, G, gap, T
    try std.testing.expectEqual(@as(usize, 4), dna.dsqrlen(&dsq));
}

test "dealign: removes gaps" {
    const allocator = std.testing.allocator;
    const dsq = [_]u8{ 0, 4, 1, 2, 4, 3 }; // A, gap, C, G, gap, T
    const result = try dna.dealign(allocator, &dsq);
    defer allocator.free(result);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, result);
}

test "validateSeq: valid DNA" {
    try std.testing.expect(dna.validateSeq("ACGTacgtNRY-"));
}

test "validateSeq: invalid character" {
    try std.testing.expect(!dna.validateSeq("ACGZ"));
}

test "amino degeneracy: B = D|N" {
    const sc = [_]f64{ 0, 0, 5.0, 0, 0, 0, 0, 0, 0, 0, 0, 3.0, 0, 0, 0, 0, 0, 0, 0, 0 };
    // B(21) = D(2)|N(11), avg = (5+3)/2 = 4.0
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), amino.avgScore(21, &sc), 1e-9);
}

test "avgScVec: fills degenerate scores" {
    var sc: [18]f64 = undefined;
    sc[0] = 1.0; // A
    sc[1] = 2.0; // C
    sc[2] = 3.0; // G
    sc[3] = 4.0; // T
    dna.avgScVec(&sc);
    // R(5) = (A+G)/2 = 2.0
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), sc[5], 1e-9);
    // N(15) = (1+2+3+4)/4 = 2.5
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), sc[15], 1e-9);
}

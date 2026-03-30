// Substitution score matrices for biological sequences.
//
// Canonical residue order (amino acid digital codes 0-19):
//   A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
//   M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19

const std = @import("std");
const alphabet_mod = @import("alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;

pub const ScoreMatrix = struct {
    /// Raw score table: data[i][j] for canonical residue digital codes i, j (0..K-1).
    data: *const [20][20]i8,
    /// Human-readable name (e.g. "BLOSUM62").
    name: []const u8,
    /// Associated alphabet.
    abc: *const Alphabet,

    /// Return score for two digital codes (canonical residues only).
    pub fn score(self: ScoreMatrix, i: u8, j: u8) i8 {
        std.debug.assert(i < 20 and j < 20);
        return self.data[i][j];
    }

    /// Return score for two ASCII residue characters.
    /// Returns error.InvalidCharacter if either character is not a canonical residue.
    pub fn scoreByChar(self: ScoreMatrix, a: u8, b: u8) !i8 {
        const i = try self.abc.encode(a);
        const j = try self.abc.encode(b);
        if (i >= self.abc.k or j >= self.abc.k) return error.InvalidCharacter;
        return self.data[i][j];
    }
};

// ---------------------------------------------------------------------------
// BLOSUM62 data
//
// Row/column order matches zeacel's amino alphabet:
//   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
// ---------------------------------------------------------------------------

const blosum62_data: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2 }, // A
    .{  0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2 }, // C
    .{ -2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3 }, // D
    .{ -1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2 }, // E
    .{ -2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3 }, // F
    .{  0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3 }, // G
    .{ -2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2 }, // H
    .{ -1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1 }, // I
    .{ -1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2 }, // K
    .{ -1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1 }, // L
    .{ -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1 }, // M
    .{ -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2 }, // N
    .{ -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3 }, // P
    .{ -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1 }, // Q
    .{ -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2 }, // R
    .{  1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2 }, // S
    .{  0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2 }, // T
    .{  0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1 }, // V
    .{ -3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2 }, // W
    .{ -2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7 }, // Y
};

/// BLOSUM62 substitution matrix (Henikoff & Henikoff, 1992).
/// Widely used for protein sequence alignment.
pub const blosum62 = ScoreMatrix{
    .data = &blosum62_data,
    .name = "BLOSUM62",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// Identity matrix: +1 for match, -1 for mismatch.
// Useful for testing and simple DNA-like scoring.
// ---------------------------------------------------------------------------

const identity_data: [20][20]i8 = blk: {
    var m: [20][20]i8 = undefined;
    var i: usize = 0;
    while (i < 20) : (i += 1) {
        var j: usize = 0;
        while (j < 20) : (j += 1) {
            m[i][j] = if (i == j) 1 else -1;
        }
    }
    break :blk m;
};

/// Simple identity substitution matrix: +1 for match, -1 for mismatch.
pub const identity = ScoreMatrix{
    .data = &identity_data,
    .name = "IDENTITY",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "blosum62: diagonal values" {
    // A-A (0,0) = 4
    try std.testing.expectEqual(@as(i8, 4), blosum62.score(0, 0));
    // C-C (1,1) = 9
    try std.testing.expectEqual(@as(i8, 9), blosum62.score(1, 1));
    // W-W (18,18) = 11 (highest diagonal)
    try std.testing.expectEqual(@as(i8, 11), blosum62.score(18, 18));
}

test "blosum62: off-diagonal values" {
    // A-C (0,1) = 0
    try std.testing.expectEqual(@as(i8, 0), blosum62.score(0, 1));
}

test "blosum62: symmetry" {
    var i: u8 = 0;
    while (i < 20) : (i += 1) {
        var j: u8 = 0;
        while (j < 20) : (j += 1) {
            try std.testing.expectEqual(blosum62.score(i, j), blosum62.score(j, i));
        }
    }
}

test "blosum62: scoreByChar valid" {
    try std.testing.expectEqual(@as(i8, 4), try blosum62.scoreByChar('A', 'A'));
    try std.testing.expectEqual(@as(i8, 11), try blosum62.scoreByChar('W', 'W'));
    try std.testing.expectEqual(@as(i8, 0), try blosum62.scoreByChar('A', 'C'));
}

test "blosum62: scoreByChar case-insensitive" {
    try std.testing.expectEqual(@as(i8, 4), try blosum62.scoreByChar('a', 'a'));
    try std.testing.expectEqual(@as(i8, 11), try blosum62.scoreByChar('w', 'W'));
}

test "blosum62: scoreByChar invalid character" {
    try std.testing.expectError(error.InvalidCharacter, blosum62.scoreByChar('!', 'A'));
    try std.testing.expectError(error.InvalidCharacter, blosum62.scoreByChar('A', '!'));
}

test "blosum62: scoreByChar non-canonical (gap) returns error" {
    // '-' encodes to code 20, which is >= abc.k (20), so it's an error
    try std.testing.expectError(error.InvalidCharacter, blosum62.scoreByChar('-', 'A'));
}

test "identity: match and mismatch" {
    try std.testing.expectEqual(@as(i8, 1), identity.score(0, 0));
    try std.testing.expectEqual(@as(i8, -1), identity.score(0, 1));
    try std.testing.expectEqual(@as(i8, 1), identity.score(19, 19));
    try std.testing.expectEqual(@as(i8, -1), identity.score(0, 19));
}

test "identity: symmetry" {
    var i: u8 = 0;
    while (i < 20) : (i += 1) {
        var j: u8 = 0;
        while (j < 20) : (j += 1) {
            try std.testing.expectEqual(identity.score(i, j), identity.score(j, i));
        }
    }
}

test "blosum62: known off-diagonal values" {
    // D-E (2,3) = 2
    try std.testing.expectEqual(@as(i8, 2), blosum62.score(2, 3));
    // F-Y (4,19) = 3
    try std.testing.expectEqual(@as(i8, 3), blosum62.score(4, 19));
    // A-R (0,14) = -1
    try std.testing.expectEqual(@as(i8, -1), blosum62.score(0, 14));
    // H-H (6,6) = 8
    try std.testing.expectEqual(@as(i8, 8), blosum62.score(6, 6));
    // P-P (12,12) = 7
    try std.testing.expectEqual(@as(i8, 7), blosum62.score(12, 12));
}

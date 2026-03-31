// Substitution score matrices for biological sequences.
//
// Canonical residue order (amino acid digital codes 0-19):
//   A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
//   M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19

const std = @import("std");
const math = std.math;
const alphabet_mod = @import("alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;
const composition = @import("composition.zig");

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

    /// Check if the matrix is symmetric (S[i][j] == S[j][i] for all i,j).
    pub fn isSymmetric(self: ScoreMatrix) bool {
        for (0..20) |i| {
            for (i + 1..20) |j| {
                if (self.data[i][j] != self.data[j][i]) return false;
            }
        }
        return true;
    }

    /// Return the maximum score value in the matrix.
    pub fn maxScore(self: ScoreMatrix) i8 {
        var m: i8 = self.data[0][0];
        for (self.data) |row| {
            for (row) |val| {
                if (val > m) m = val;
            }
        }
        return m;
    }

    /// Return the minimum score value in the matrix.
    pub fn minScore(self: ScoreMatrix) i8 {
        var m: i8 = self.data[0][0];
        for (self.data) |row| {
            for (row) |val| {
                if (val < m) m = val;
            }
        }
        return m;
    }

    /// Compute the expected score E = sum_ij f_i * S_ij * f_j given background frequencies.
    pub fn expectedScore(self: ScoreMatrix, bg: *const [20]f64) f64 {
        var e: f64 = 0;
        for (0..20) |i| {
            for (0..20) |j| {
                e += bg[i] * @as(f64, @floatFromInt(self.data[i][j])) * bg[j];
            }
        }
        return e;
    }

    /// Compute relative entropy H = sum_ij p_ij * log2(p_ij / (f_i * f_j))
    /// given background frequencies, where p_ij = f_i * f_j * exp(lambda * s_ij).
    /// Returns the relative entropy in bits, or null if lambda cannot be estimated.
    pub fn relativeEntropy(self: ScoreMatrix, bg: *const [20]f64, lambda: f64) f64 {
        // First compute the normalizing constant Z = sum_ij f_i * f_j * exp(lambda * s_ij)
        var z: f64 = 0;
        for (0..20) |i| {
            for (0..20) |j| {
                z += bg[i] * bg[j] * @exp(lambda * @as(f64, @floatFromInt(self.data[i][j])));
            }
        }

        // H = sum_ij p_ij * log2(p_ij / (f_i * f_j))
        //   = sum_ij p_ij * (lambda * s_ij * log2(e) - log2(Z))
        var h: f64 = 0;
        for (0..20) |i| {
            for (0..20) |j| {
                const s: f64 = @floatFromInt(self.data[i][j]);
                const p_ij = bg[i] * bg[j] * @exp(lambda * s) / z;
                if (p_ij > 0) {
                    h += p_ij * (lambda * s * math.log2e - @log2(z));
                }
            }
        }
        return h;
    }
};

/// Look up a named score matrix. Returns null if name is not recognized.
pub fn getByName(name: []const u8) ?ScoreMatrix {
    if (std.ascii.eqlIgnoreCase(name, "BLOSUM62") or std.ascii.eqlIgnoreCase(name, "BL62")) return blosum62;
    if (std.ascii.eqlIgnoreCase(name, "BLOSUM45")) return blosum45;
    if (std.ascii.eqlIgnoreCase(name, "BLOSUM80")) return blosum80;
    if (std.ascii.eqlIgnoreCase(name, "PAM30")) return pam30;
    if (std.ascii.eqlIgnoreCase(name, "PAM70")) return pam70;
    if (std.ascii.eqlIgnoreCase(name, "PAM120")) return pam120;
    if (std.ascii.eqlIgnoreCase(name, "PAM250")) return pam250;
    if (std.ascii.eqlIgnoreCase(name, "IDENTITY")) return identity;
    return null;
}

// ---------------------------------------------------------------------------
// BLOSUM62 data
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

pub const blosum62 = ScoreMatrix{
    .data = &blosum62_data,
    .name = "BLOSUM62",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// BLOSUM45 data
// ---------------------------------------------------------------------------

const blosum45_data: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  5, -1, -2, -1, -2,  0, -2, -1, -1, -1, -1, -1, -1, -1, -2,  1,  0,  0, -2, -2 }, // A
    .{ -1, 12, -3, -3, -2, -3, -3, -3, -3, -2, -2, -2, -4, -3, -3, -1, -1, -1, -5, -3 }, // C
    .{ -2, -3,  7,  2, -4, -1,  0, -4,  0, -3, -3,  2, -1,  0, -1,  0, -1, -3, -4, -2 }, // D
    .{ -1, -3,  2,  6, -3, -2,  0, -3,  1, -2, -2,  0,  0,  2,  0,  0, -1, -3, -3, -2 }, // E
    .{ -2, -2, -4, -3,  8, -3, -2,  0, -3,  1,  0, -2, -3, -4, -2, -2, -1,  0,  1,  3 }, // F
    .{  0, -3, -1, -2, -3,  7, -2, -4, -2, -3, -2,  0, -2, -2, -2,  0, -2, -3, -2, -3 }, // G
    .{ -2, -3,  0,  0, -2, -2, 10, -3, -1, -2,  0,  1, -2,  1,  0, -1, -2, -3, -3,  2 }, // H
    .{ -1, -3, -4, -3,  0, -4, -3,  5, -3,  2,  2, -2, -2, -2, -3, -2, -1,  3, -2,  0 }, // I
    .{ -1, -3,  0,  1, -3, -2, -1, -3,  5, -3, -1,  0, -1,  1,  3, -1, -1, -2, -2, -1 }, // K
    .{ -1, -2, -3, -2,  1, -3, -2,  2, -3,  5,  2, -3, -3, -2, -2, -3, -1,  1, -2,  0 }, // L
    .{ -1, -2, -3, -2,  0, -2,  0,  2, -1,  2,  6, -2, -2,  0, -1, -2, -1,  1, -2,  0 }, // M
    .{ -1, -2,  2,  0, -2,  0,  1, -2,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2 }, // N
    .{ -1, -4, -1,  0, -3, -2, -2, -2, -1, -3, -2, -2,  9, -1, -2, -1, -1, -3, -3, -3 }, // P
    .{ -1, -3,  0,  2, -4, -2,  1, -2,  1, -2,  0,  0, -1,  6,  1,  0, -1, -3, -2, -1 }, // Q
    .{ -2, -3, -1,  0, -2, -2,  0, -3,  3, -2, -1,  0, -2,  1,  7, -1, -1, -2, -2, -1 }, // R
    .{  1, -1,  0,  0, -2,  0, -1, -2, -1, -3, -2,  1, -1,  0, -1,  4,  2, -1, -4, -2 }, // S
    .{  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  2,  5,  0, -3, -1 }, // T
    .{  0, -1, -3, -3,  0, -3, -3,  3, -2,  1,  1, -3, -3, -3, -2, -1,  0,  5, -3, -1 }, // V
    .{ -2, -5, -4, -3,  1, -2, -3, -2, -2, -2, -2, -4, -3, -2, -2, -4, -3, -3, 15,  3 }, // W
    .{ -2, -3, -2, -2,  3, -3,  2,  0, -1,  0,  0, -2, -3, -1, -1, -2, -1, -1,  3,  8 }, // Y
};

pub const blosum45 = ScoreMatrix{
    .data = &blosum45_data,
    .name = "BLOSUM45",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// BLOSUM80 data
// ---------------------------------------------------------------------------

const blosum80_data: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  7, -1, -3, -2, -4,  0, -3, -3, -1, -3, -2, -3, -1, -2, -3,  2,  0, -1, -5, -4 }, // A
    .{ -1, 13, -7, -7, -4, -6, -7, -2, -6, -3, -3, -5, -6, -5, -6, -2, -2, -2, -5, -5 }, // C
    .{ -3, -7, 10,  2, -6, -3, -2, -5, -2, -7, -6,  2, -3, -1, -3, -1, -2, -6, -8, -6 }, // D
    .{ -2, -7,  2,  8, -6, -4,  0, -6,  1, -6, -4, -1, -2,  3, -1, -1, -2, -4, -6, -5 }, // E
    .{ -4, -4, -6, -6, 10, -6, -2, -1, -5,  2,  0, -6, -6, -5, -5, -4, -4, -2,  0,  4 }, // F
    .{  0, -6, -3, -4, -6,  9, -4, -7, -3, -7, -5, -1, -5, -4, -4, -1, -3, -6, -6, -6 }, // G
    .{ -3, -7, -2,  0, -2, -4, 12, -6, -1, -5, -4,  1, -4,  1,  0, -2, -3, -5, -4,  3 }, // H
    .{ -3, -2, -5, -6, -1, -7, -6,  7, -5,  2,  2, -6, -5, -5, -5, -4, -2,  4, -5, -3 }, // I
    .{ -1, -6, -2,  1, -5, -3, -1, -5,  8, -4, -3,  0, -2,  2,  3, -1, -1, -4, -6, -4 }, // K
    .{ -3, -3, -7, -6,  2, -7, -5,  2, -4,  6,  3, -6, -5, -4, -4, -4, -3,  1, -4, -2 }, // L
    .{ -2, -3, -6, -4,  0, -5, -4,  2, -3,  3,  9, -4, -4, -1, -3, -3, -1,  1, -3, -3 }, // M
    .{ -3, -5,  2, -1, -6, -1,  1, -6,  0, -6, -4,  9, -4,  0, -1,  1,  0, -5, -7, -4 }, // N
    .{ -1, -6, -3, -2, -6, -5, -4, -5, -2, -5, -4, -4, 12, -3, -3, -2, -3, -4, -7, -6 }, // P
    .{ -2, -5, -1,  3, -5, -4,  1, -5,  2, -4, -1,  0, -3,  9,  1, -1, -1, -4, -4, -3 }, // Q
    .{ -3, -6, -3, -1, -5, -4,  0, -5,  3, -4, -3, -1, -3,  1,  9, -2, -2, -4, -5, -4 }, // R
    .{  2, -2, -1, -1, -4, -1, -2, -4, -1, -4, -3,  1, -2, -1, -2,  7,  2, -3, -6, -3 }, // S
    .{  0, -2, -2, -2, -4, -3, -3, -2, -1, -3, -1,  0, -3, -1, -2,  2,  8, -1, -5, -3 }, // T
    .{ -1, -2, -6, -4, -2, -6, -5,  4, -4,  1,  1, -5, -4, -4, -4, -3, -1,  7, -5, -3 }, // V
    .{ -5, -5, -8, -6,  0, -6, -4, -5, -6, -4, -3, -7, -7, -4, -5, -6, -5, -5, 16,  3 }, // W
    .{ -4, -5, -6, -5,  4, -6,  3, -3, -4, -2, -3, -4, -6, -3, -4, -3, -3, -3,  3, 11 }, // Y
};

pub const blosum80 = ScoreMatrix{
    .data = &blosum80_data,
    .name = "BLOSUM80",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// PAM30 data
// ---------------------------------------------------------------------------

const pam30_data: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  6, -6, -3, -2, -8,  0, -7, -5, -7, -6, -5, -4, -2, -4, -7,  0, -1, -2,-13, -8 }, // A
    .{ -6, 10,-14,-14,-13, -9, -7, -6,-14,-15,-13,-11, -8,-14, -8, -3, -8, -6,-15, -4 }, // C
    .{ -3,-14,  8,  2,-15, -3, -4, -7, -4,-12,-11,  2, -8, -2,-10, -4, -5, -8,-15,-11 }, // D
    .{ -2,-14,  2,  8,-14, -4, -5, -5, -4, -9, -7, -2, -5,  1, -9, -4, -6, -6,-17, -8 }, // E
    .{ -8,-13,-15,-14,  9, -9, -6, -2,-14, -3, -4,-11,-10,-13, -9, -6, -9, -8, -4,  2 }, // F
    .{  0, -9, -3, -4, -9,  6, -9,-11, -7,-10, -8, -3, -6, -7, -9, -2, -6, -5,-15,-14 }, // G
    .{ -7, -7, -4, -5, -6, -9,  9, -9, -6, -6,-10, -3, -4,  1, -2, -6, -7, -6, -7, -3 }, // H
    .{ -5, -6, -7, -5, -2,-11, -9,  8, -6, -1, -1, -5, -8, -8, -5, -7, -2,  2,-14, -6 }, // I
    .{ -7,-14, -4, -4,-14, -7, -6, -6,  7, -8, -2, -1, -6, -3,  0, -4, -3, -9, -6, -9 }, // K
    .{ -6,-15,-12, -9, -3,-10, -6, -1, -8,  7,  1, -7, -7, -5, -8, -8, -7, -2, -6, -7 }, // L
    .{ -5,-13,-11, -7, -4, -8,-10, -1, -2,  1, 11, -9, -8, -4, -4, -5, -4, -1,-13,-11 }, // M
    .{ -4,-11,  2, -2,-11, -3, -3, -5, -1, -7, -9,  8, -6, -3, -6,  0, -2, -8, -8, -4 }, // N
    .{ -2, -8, -8, -5,-10, -6, -4, -8, -6, -7, -8, -6,  8, -3, -4, -2, -4, -6,-14,-13 }, // P
    .{ -4,-14, -2,  1,-13, -7,  1, -8, -3, -5, -4, -3, -3,  8, -2, -5, -5, -7,-13,-12 }, // Q
    .{ -7, -8,-10, -9, -9, -9, -2, -5,  0, -8, -4, -6, -4, -2,  8, -3, -6, -8, -2,-10 }, // R
    .{  0, -3, -4, -4, -6, -2, -6, -7, -4, -8, -5,  0, -2, -5, -3,  6,  0, -6, -5, -7 }, // S
    .{ -1, -8, -5, -6, -9, -6, -7, -2, -3, -7, -4, -2, -4, -5, -6,  0,  7, -3,-13, -6 }, // T
    .{ -2, -6, -8, -6, -8, -5, -6,  2, -9, -2, -1, -8, -6, -7, -8, -6, -3,  7,-15, -7 }, // V
    .{-13,-15,-15,-17, -4,-15, -7,-14, -6, -6,-13, -8,-14,-13, -2, -5,-13,-15, 13, -5 }, // W
    .{ -8, -4,-11, -8,  2,-14, -3, -6, -9, -7,-11, -4,-13,-12,-10, -7, -6, -7, -5, 10 }, // Y
};

pub const pam30 = ScoreMatrix{
    .data = &pam30_data,
    .name = "PAM30",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// PAM70 data
// ---------------------------------------------------------------------------

const pam70_data: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  5, -4, -1, -1, -6,  0, -4, -2, -4, -4, -3, -2, -1, -2, -4,  1,  1, -1, -9, -5 }, // A
    .{ -4,  9, -9, -9, -8, -6, -5, -4, -9,-10, -9, -7, -5, -9, -5, -1, -5, -4,-11, -2 }, // C
    .{ -1, -9,  7,  3,-10, -1, -1, -5, -2, -8, -7,  3, -4, -1, -5, -1, -2, -5,-10, -7 }, // D
    .{ -1, -9,  3,  7, -9, -3, -2, -3, -2, -6, -4, -1, -3,  2, -5, -2, -3, -4,-11, -6 }, // E
    .{ -6, -8,-10, -9,  8, -7, -4, -1, -9, -1, -2, -6, -7, -9, -7, -4, -6, -5, -2,  4 }, // F
    .{  0, -6, -1, -3, -7,  6, -6, -6, -5, -7, -6, -1, -3, -4, -6, -1, -3, -3,-10, -9 }, // G
    .{ -4, -5, -1, -2, -4, -6,  8, -6, -3, -4, -6, -1, -2,  1, -1, -3, -4, -4, -5, -1 }, // H
    .{ -2, -4, -5, -3, -1, -6, -6,  7, -4, -1,  0, -3, -5, -5, -3, -4, -1,  3, -9, -5 }, // I
    .{ -4, -9, -2, -2, -9, -5, -3, -4,  6, -5, -1, -1, -4, -1,  1, -2, -1, -6, -7, -7 }, // K
    .{ -4,-10, -8, -6, -1, -7, -4, -1, -5,  6,  1, -5, -5, -3, -6, -6, -4, -1, -4, -4 }, // L
    .{ -3, -9, -7, -4, -2, -6, -6,  0, -1,  1,  8, -5, -5, -2, -2, -3, -2,  0, -8, -7 }, // M
    .{ -2, -7,  3, -1, -6, -1, -1, -3, -1, -5, -5,  6, -3, -1, -3,  1, -1, -5, -6, -3 }, // N
    .{ -1, -5, -4, -3, -7, -3, -2, -5, -4, -5, -5, -3,  7, -2, -2, -1, -2, -3, -9, -9 }, // P
    .{ -2, -9, -1,  2, -9, -4,  1, -5, -1, -3, -2, -1, -2,  7, -1, -3, -3, -4, -8, -8 }, // Q
    .{ -4, -5, -5, -5, -7, -6, -1, -3,  1, -6, -2, -3, -2, -1,  8, -2, -4, -5, -1, -7 }, // R
    .{  1, -1, -1, -2, -4, -1, -3, -4, -2, -6, -3,  1, -1, -3, -2,  5,  2, -3, -3, -5 }, // S
    .{  1, -5, -2, -3, -6, -3, -4, -1, -1, -4, -2, -1, -2, -3, -4,  2,  6, -1, -8, -4 }, // T
    .{ -1, -4, -5, -4, -5, -3, -4,  3, -6, -1,  0, -5, -3, -4, -5, -3, -1,  6,-10, -5 }, // V
    .{ -9,-11,-10,-11, -2,-10, -5, -9, -7, -4, -8, -6, -9, -8, -1, -3, -8,-10, 13, -3 }, // W
    .{ -5, -2, -7, -6,  4, -9, -1, -5, -7, -4, -7, -3, -9, -8, -7, -5, -4, -5, -3,  9 }, // Y
};

pub const pam70 = ScoreMatrix{
    .data = &pam70_data,
    .name = "PAM70",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// PAM120 and PAM250 (abbreviated — diagonal + key off-diagonals)
// ---------------------------------------------------------------------------

const pam120_data: [20][20]i8 = .{
    .{  3, -3,  0,  0, -4,  1, -3, -1, -2, -3, -2, -1,  1, -1, -3,  1,  1,  0, -7, -4 }, // A
    .{ -3,  9, -7, -7, -6, -5, -4, -3, -7, -7, -6, -5, -4, -7, -4,  0, -3, -3, -8, -1 }, // C
    .{  0, -7,  5,  3, -7, -1,  0, -3, -1, -5, -4,  2, -3,  1, -3,  0, -1, -3, -8, -5 }, // D
    .{  0, -7,  3,  5, -7, -2, -1, -3, -1, -4, -3,  1, -2,  2, -3, -1, -2, -3, -8, -5 }, // E
    .{ -4, -6, -7, -7,  8, -5, -2,  0, -7,  0, -1, -4, -5, -6, -5, -3, -4, -3, -1,  4 }, // F
    .{  1, -5, -1, -2, -5,  5, -4, -4, -3, -5, -4,  0, -2, -3, -4,  1, -1, -2, -8, -6 }, // G
    .{ -3, -4,  0, -1, -2, -4,  7, -4, -2, -3, -4,  2, -1,  3,  1, -2, -3, -3, -3, -1 }, // H
    .{ -1, -3, -3, -3,  0, -4, -4,  6, -3,  1,  1, -2, -3, -3, -2, -2,  0,  3, -6, -2 }, // I
    .{ -2, -7, -1, -1, -7, -3, -2, -3,  5, -4,  0,  1, -2,  0,  2, -1, -1, -4, -5, -5 }, // K
    .{ -3, -7, -5, -4,  0, -5, -3,  1, -4,  5,  3, -4, -3, -2, -4, -4, -3,  1, -3, -2 }, // L
    .{ -2, -6, -4, -3, -1, -4, -4,  1,  0,  3,  8, -3, -3, -1, -1, -2, -1,  1, -6, -4 }, // M
    .{ -1, -5,  2,  1, -4,  0,  2, -2,  1, -4, -3,  4, -2,  0, -1,  1,  0, -3, -4, -2 }, // N
    .{  1, -4, -3, -2, -5, -2, -1, -3, -2, -3, -3, -2,  6, -1, -1,  1, -1, -2, -7, -6 }, // P
    .{ -1, -7,  1,  2, -6, -3,  3, -3,  0, -2, -1,  0, -1,  6,  1, -2, -2, -3, -6, -5 }, // Q
    .{ -3, -4, -3, -3, -5, -4,  1, -2,  2, -4, -1, -1, -1,  1,  6, -1, -2, -3,  1, -5 }, // R
    .{  1,  0,  0, -1, -3,  1, -2, -2, -1, -4, -2,  1,  1, -2, -1,  3,  2, -2, -2, -3 }, // S
    .{  1, -3, -1, -2, -4, -1, -3,  0, -1, -3, -1,  0, -1, -2, -2,  2,  4,  0, -6, -3 }, // T
    .{  0, -3, -3, -3, -3, -2, -3,  3, -4,  1,  1, -3, -2, -3, -3, -2,  0,  5, -8, -3 }, // V
    .{ -7, -8, -8, -8, -1, -8, -3, -6, -5, -3, -6, -4, -7, -6,  1, -2, -6, -8, 12, -2 }, // W
    .{ -4, -1, -5, -5,  4, -6, -1, -2, -5, -2, -4, -2, -6, -5, -5, -3, -3, -3, -2,  8 }, // Y
};

pub const pam120 = ScoreMatrix{
    .data = &pam120_data,
    .name = "PAM120",
    .abc = &alphabet_mod.amino,
};

const pam250_data: [20][20]i8 = .{
    .{  2, -2,  0,  0, -3,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3 }, // A
    .{ -2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0 }, // C
    .{  0, -5,  4,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4 }, // D
    .{  0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4 }, // E
    .{ -3, -4, -6, -5,  9, -5, -2,  1, -5,  2,  0, -3, -5, -5, -4, -3, -3, -1,  0,  7 }, // F
    .{  1, -3,  1,  0, -5,  5, -2, -3, -2, -4, -3,  0,  0, -1, -3,  1,  0, -1, -7, -5 }, // G
    .{ -1, -3,  1,  1, -2, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0 }, // H
    .{ -1, -2, -2, -2,  1, -3, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1 }, // I
    .{ -1, -5,  0,  0, -5, -2,  0, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4 }, // K
    .{ -2, -6, -4, -3,  2, -4, -2,  2, -3,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1 }, // L
    .{ -1, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6, -2, -2, -1,  0, -2, -1,  2, -4, -2 }, // M
    .{  0, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  2,  0,  1,  0,  1,  0, -2, -4, -2 }, // N
    .{  1, -3, -1, -1, -5,  0,  0, -2, -1, -3, -2,  0,  6,  0,  0,  1,  0, -1, -6, -5 }, // P
    .{  0, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4,  1, -1, -1, -2, -5, -4 }, // Q
    .{ -2, -4, -1, -1, -4, -3,  2, -2,  3, -3,  0,  0,  0,  1,  6,  0, -1, -2,  2, -4 }, // R
    .{  1,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1, -1,  0,  2,  1, -1, -2, -3 }, // S
    .{  1, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,  0, -5, -3 }, // T
    .{  0, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4, -6, -2 }, // V
    .{ -6, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17,  0 }, // W
    .{ -3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, 10 }, // Y
};

pub const pam250 = ScoreMatrix{
    .data = &pam250_data,
    .name = "PAM250",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// Identity matrix: +1 for match, -1 for mismatch.
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

pub const identity = ScoreMatrix{
    .data = &identity_data,
    .name = "IDENTITY",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "blosum62: diagonal values" {
    try std.testing.expectEqual(@as(i8, 4), blosum62.score(0, 0));
    try std.testing.expectEqual(@as(i8, 9), blosum62.score(1, 1));
    try std.testing.expectEqual(@as(i8, 11), blosum62.score(18, 18));
}

test "blosum62: off-diagonal values" {
    try std.testing.expectEqual(@as(i8, 0), blosum62.score(0, 1));
}

test "blosum62: symmetry" {
    try std.testing.expect(blosum62.isSymmetric());
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
    try std.testing.expectError(error.InvalidCharacter, blosum62.scoreByChar('-', 'A'));
}

test "identity: match and mismatch" {
    try std.testing.expectEqual(@as(i8, 1), identity.score(0, 0));
    try std.testing.expectEqual(@as(i8, -1), identity.score(0, 1));
    try std.testing.expectEqual(@as(i8, 1), identity.score(19, 19));
    try std.testing.expectEqual(@as(i8, -1), identity.score(0, 19));
}

test "identity: symmetry" {
    try std.testing.expect(identity.isSymmetric());
}

test "blosum62: known off-diagonal values" {
    try std.testing.expectEqual(@as(i8, 2), blosum62.score(2, 3));
    try std.testing.expectEqual(@as(i8, 3), blosum62.score(4, 19));
    try std.testing.expectEqual(@as(i8, -1), blosum62.score(0, 14));
    try std.testing.expectEqual(@as(i8, 8), blosum62.score(6, 6));
    try std.testing.expectEqual(@as(i8, 7), blosum62.score(12, 12));
}

test "maxScore and minScore" {
    try std.testing.expectEqual(@as(i8, 11), blosum62.maxScore());
    try std.testing.expectEqual(@as(i8, -4), blosum62.minScore());
}

test "expectedScore: BLOSUM62 with BL62 background is negative" {
    const e = blosum62.expectedScore(&composition.bl62);
    try std.testing.expect(e < 0);
}

test "getByName: known matrices" {
    try std.testing.expect(getByName("BLOSUM62") != null);
    try std.testing.expect(getByName("blosum62") != null);
    try std.testing.expect(getByName("BL62") != null);
    try std.testing.expect(getByName("BLOSUM45") != null);
    try std.testing.expect(getByName("BLOSUM80") != null);
    try std.testing.expect(getByName("PAM30") != null);
    try std.testing.expect(getByName("PAM70") != null);
    try std.testing.expect(getByName("PAM120") != null);
    try std.testing.expect(getByName("PAM250") != null);
    try std.testing.expect(getByName("IDENTITY") != null);
    try std.testing.expect(getByName("UNKNOWN") == null);
}

test "blosum45: symmetry" {
    try std.testing.expect(blosum45.isSymmetric());
}

test "blosum80: symmetry" {
    try std.testing.expect(blosum80.isSymmetric());
}

test "pam120: symmetry" {
    try std.testing.expect(pam120.isSymmetric());
}

test "pam250: symmetry" {
    try std.testing.expect(pam250.isSymmetric());
}

test "relativeEntropy: BLOSUM62 is positive" {
    const h = blosum62.relativeEntropy(&composition.bl62, 0.3466);
    try std.testing.expect(h > 0);
}

// Sequence distance calculations for aligned digital sequences.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("msa.zig").Msa;
const Alphabet = @import("alphabet.zig").Alphabet;

/// Percent identity between two aligned digital sequences.
/// Gaps in either sequence are excluded from comparison.
pub fn percentIdentity(abc: *const Alphabet, seq1: []const u8, seq2: []const u8) f64 {
    var matches: u32 = 0;
    var compared: u32 = 0;
    for (seq1, seq2) |a, b| {
        if (!abc.isGap(a) and !abc.isGap(b)) {
            compared += 1;
            if (a == b) matches += 1;
        }
    }
    if (compared == 0) return 0.0;
    return @as(f64, @floatFromInt(matches)) / @as(f64, @floatFromInt(compared));
}

/// Generalized Jukes-Cantor distance correction:
///   d = -((K-1)/K) * ln(1 - D*K/(K-1))
/// where D = 1 - percent_identity and K = alphabet_size.
/// For DNA, K=4 gives the classic -3/4 * ln(1 - 4/3 * D).
/// For protein, K=20. Returns inf if the argument to ln is <= 0.
pub fn jukesCantor(pid: f64, alphabet_size: u32) f64 {
    const k: f64 = @floatFromInt(alphabet_size);
    const d = 1.0 - pid;
    const scale = (k - 1.0) / k;
    const arg = 1.0 - d / scale;
    if (arg <= 0.0) return std.math.inf(f64);
    return -scale * @log(arg);
}

/// Kimura protein distance correction: d = -ln(1 - D - 0.2*D^2)
/// where D = 1 - percent_identity. Returns inf if argument <= 0.
pub fn kimura(pid: f64) f64 {
    const d = 1.0 - pid;
    const arg = 1.0 - d - 0.2 * d * d;
    if (arg <= 0.0) return std.math.inf(f64);
    return -@log(arg);
}

/// Correction model for pairwise distance calculation.
pub const Correction = enum { none, jukes_cantor, kimura };

/// Compute pairwise distance matrix from an MSA.
/// Returns a flat array of n*n f64 values (row-major). Caller owns memory.
pub fn pairwiseDistanceMatrix(allocator: Allocator, m: Msa, correction: Correction) ![]f64 {
    const n = m.nseq();
    const matrix = try allocator.alloc(f64, n * n);
    @memset(matrix, 0.0);
    for (0..n) |i| {
        for (i..n) |j| {
            const pid = percentIdentity(m.abc, m.seqs[i], m.seqs[j]);
            const dist = switch (correction) {
                .none => 1.0 - pid,
                .jukes_cantor => jukesCantor(pid, m.abc.k),
                .kimura => kimura(pid),
            };
            matrix[i * n + j] = dist;
            matrix[j * n + i] = dist;
        }
    }
    return matrix;
}

// --- Tests ---

test "percentIdentity: identical sequences returns 1.0" {
    const abc = &@import("alphabet.zig").dna;
    // A=0, C=1, G=2, T=3 — no gaps
    const seq = &[_]u8{ 0, 1, 2, 3 };
    const pid = percentIdentity(abc, seq, seq);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), pid, 1e-10);
}

test "percentIdentity: completely different sequences returns 0.0" {
    const abc = &@import("alphabet.zig").dna;
    // A=0 vs T=3 at every position
    const seq1 = &[_]u8{ 0, 0, 0, 0 };
    const seq2 = &[_]u8{ 3, 3, 3, 3 };
    const pid = percentIdentity(abc, seq1, seq2);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), pid, 1e-10);
}

test "percentIdentity: gaps are excluded from comparison" {
    const abc = &@import("alphabet.zig").dna;
    // gap code for DNA is abc.k = 4
    // seq1: A A - A  (positions 0,1,3 compared; position 2 is gap in seq2)
    // seq2: A A A -  (position 3 is gap in seq2, position 2 is gap in seq1)
    // Both have gaps at different places: positions 2 and 3 each have one gap.
    // Positions 0 and 1 are comparable: A==A and A==A -> 2 matches / 2 compared = 1.0
    const gap: u8 = 4;
    const seq1 = &[_]u8{ 0, 0, gap, 0 };
    const seq2 = &[_]u8{ 0, 0, 0, gap };
    const pid = percentIdentity(abc, seq1, seq2);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), pid, 1e-10);
}

test "percentIdentity: all gaps returns 0.0" {
    const abc = &@import("alphabet.zig").dna;
    const gap: u8 = 4;
    const seq1 = &[_]u8{ gap, gap };
    const seq2 = &[_]u8{ gap, gap };
    const pid = percentIdentity(abc, seq1, seq2);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), pid, 1e-10);
}

test "jukesCantor: DNA pid=1.0 gives distance 0.0" {
    const d = jukesCantor(1.0, 4);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), d, 1e-10);
}

test "jukesCantor: DNA pid=0.7 gives known value" {
    // D = 0.3; d = -0.75 * ln(1 - 4/3 * 0.3) = -0.75 * ln(0.6)
    const expected = -0.75 * @log(1.0 - (4.0 / 3.0) * 0.3);
    const d = jukesCantor(0.7, 4);
    try std.testing.expectApproxEqAbs(expected, d, 1e-10);
}

test "jukesCantor: DNA pid<=0.25 gives inf" {
    try std.testing.expectEqual(std.math.inf(f64), jukesCantor(0.25, 4));
    try std.testing.expectEqual(std.math.inf(f64), jukesCantor(0.0, 4));
}

test "jukesCantor: protein K=20 gives correct value" {
    // d = -(19/20) * ln(1 - D*20/19), D = 0.5
    const expected = -(19.0 / 20.0) * @log(1.0 - 0.5 * 20.0 / 19.0);
    const d = jukesCantor(0.5, 20);
    try std.testing.expectApproxEqAbs(expected, d, 1e-10);
}

test "kimura: pid=1.0 gives distance 0.0" {
    const d = kimura(1.0);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), d, 1e-10);
}

test "kimura: pid below threshold gives inf" {
    // arg = 1 - D - 0.2*D^2 <= 0
    // At D=1.0 (pid=0.0): arg = 1 - 1 - 0.2 = -0.2 <= 0 -> inf
    try std.testing.expectEqual(std.math.inf(f64), kimura(0.0));
}

test "pairwiseDistanceMatrix: 3 seqs, symmetric, diagonal zero" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "AAAA", "AATT", "TTTT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    const mat = try pairwiseDistanceMatrix(allocator, msa, .none);
    defer allocator.free(mat);

    const n = msa.nseq(); // 3

    // Diagonal must be 0.0
    for (0..n) |i| {
        try std.testing.expectApproxEqAbs(@as(f64, 0.0), mat[i * n + i], 1e-10);
    }

    // Matrix must be symmetric
    for (0..n) |i| {
        for (0..n) |j| {
            try std.testing.expectApproxEqAbs(mat[i * n + j], mat[j * n + i], 1e-10);
        }
    }

    // Verify known values:
    // s1 vs s2: "AAAA" vs "AATT" -> 2 matches / 4 compared = 0.5 pid -> dist = 0.5
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), mat[0 * n + 1], 1e-10);
    // s1 vs s3: "AAAA" vs "TTTT" -> 0 matches / 4 compared = 0.0 pid -> dist = 1.0
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), mat[0 * n + 2], 1e-10);
}

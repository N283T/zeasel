// Sequence distance calculations for aligned digital sequences.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("msa.zig").Msa;
const Alphabet = @import("alphabet.zig").Alphabet;

/// Percent identity between two aligned digital sequences.
/// Uses Easel's definition: denominator is MIN(len1, len2) where
/// len1/len2 are the number of non-gap residues in each sequence.
/// This matches esl_dst_XPairId() used by HMMER.
pub fn percentIdentity(abc: *const Alphabet, seq1: []const u8, seq2: []const u8) f64 {
    var matches: u32 = 0;
    var len1: u32 = 0;
    var len2: u32 = 0;
    for (seq1, seq2) |a, b| {
        const a_res = abc.isResidue(a);
        const b_res = abc.isResidue(b);
        if (a_res) len1 += 1;
        if (b_res) len2 += 1;
        if (a_res and b_res and a == b) matches += 1;
    }
    const denom = @min(len1, len2);
    if (denom == 0) return 0.0;
    return @as(f64, @floatFromInt(matches)) / @as(f64, @floatFromInt(denom));
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

/// Result of Jukes-Cantor distance estimation with variance.
pub const JukesCantorResult = struct {
    distance: f64,
    variance: f64,
};

/// Generalized Jukes-Cantor distance correction with variance estimate.
/// Returns both the corrected distance and its variance.
/// Variance formula: variance = exp(2*K*d/(K-1)) * D * (1-D) / N
/// where D = 1 - pid, d = JC distance, K = alphabet_size, N = n_compared.
/// See Easel esl_dst_jukescantor() for reference.
pub fn jukesCantorWithVariance(pid: f64, alphabet_size: u32, n_compared: u32) JukesCantorResult {
    if (n_compared == 0) return .{ .distance = std.math.inf(f64), .variance = std.math.inf(f64) };
    const k: f64 = @floatFromInt(alphabet_size);
    const big_d = 1.0 - pid;
    const scale = (k - 1.0) / k;
    const arg = 1.0 - big_d / scale;
    if (arg <= 0.0) return .{ .distance = std.math.inf(f64), .variance = std.math.inf(f64) };
    const distance = -scale * @log(arg);
    const n: f64 = @floatFromInt(n_compared);
    const variance = @exp(2.0 * k * distance / (k - 1.0)) * big_d * (1.0 - big_d) / n;
    return .{ .distance = distance, .variance = variance };
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

/// Compute average pairwise fractional identity for an MSA.
/// For small MSAs where N*(N-1)/2 <= max_pairs, computes all pairs exhaustively.
/// For large MSAs, subsamples up to max_pairs random pairs using a fixed seed
/// for reproducibility. Returns 1.0 for single-sequence MSAs by convention.
/// See Easel esl_dst_XAverageId() for reference.
pub fn averageIdentity(allocator: Allocator, m: Msa, max_pairs: usize) !f64 {
    _ = allocator;
    const n = m.nseq();
    if (n <= 1) return 1.0;

    const total_pairs = n * (n - 1) / 2;
    if (total_pairs <= max_pairs) {
        // Exhaustive calculation over all pairs.
        var sum: f64 = 0.0;
        for (0..n) |i| {
            for (i + 1..n) |j| {
                sum += percentIdentity(m.abc, m.seqs[i], m.seqs[j]);
            }
        }
        return sum / @as(f64, @floatFromInt(total_pairs));
    } else {
        // Stochastic subsampling with fixed seed for reproducibility.
        const Random = @import("util/random.zig").Random;
        var rng = Random.init(42);
        var sum: f64 = 0.0;
        for (0..max_pairs) |_| {
            const i = rng.uniformInt(@intCast(n));
            var j = rng.uniformInt(@intCast(n - 1));
            if (j >= i) j += 1;
            sum += percentIdentity(m.abc, m.seqs[i], m.seqs[j]);
        }
        return sum / @as(f64, @floatFromInt(max_pairs));
    }
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

test "percentIdentity: uses MIN(len1,len2) denominator" {
    const abc = &@import("alphabet.zig").dna;
    // seq1: A A - A  (len1=3)
    // seq2: A A A -  (len2=3)
    // matches at cols 0,1 = 2 (col 2: seq1 has gap; col 3: seq2 has gap)
    // pid = 2 / MIN(3,3) = 2/3
    const gap: u8 = 4;
    const seq1 = &[_]u8{ 0, 0, gap, 0 };
    const seq2 = &[_]u8{ 0, 0, 0, gap };
    const pid = percentIdentity(abc, seq1, seq2);
    try std.testing.expectApproxEqAbs(2.0 / 3.0, pid, 1e-10);
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

test "jukesCantorWithVariance: pid=1.0 gives zero distance and zero variance" {
    const result = jukesCantorWithVariance(1.0, 4, 100);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.distance, 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.variance, 1e-10);
}

test "jukesCantorWithVariance: distance matches jukesCantor" {
    const pid = 0.7;
    const result = jukesCantorWithVariance(pid, 4, 200);
    const expected_dist = jukesCantor(pid, 4);
    try std.testing.expectApproxEqAbs(expected_dist, result.distance, 1e-10);
    const expected_var = @exp(2.0 * 4.0 * expected_dist / 3.0) * 0.3 * 0.7 / 200.0;
    try std.testing.expectApproxEqAbs(expected_var, result.variance, 1e-10);
}

test "jukesCantorWithVariance: saturated returns inf for both" {
    const result = jukesCantorWithVariance(0.0, 4, 100);
    try std.testing.expectEqual(std.math.inf(f64), result.distance);
    try std.testing.expectEqual(std.math.inf(f64), result.variance);
}

test "jukesCantorWithVariance: n_compared=0 returns inf" {
    const result = jukesCantorWithVariance(0.8, 4, 0);
    try std.testing.expectEqual(std.math.inf(f64), result.distance);
    try std.testing.expectEqual(std.math.inf(f64), result.variance);
}

test "averageIdentity: identical sequences returns 1.0" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "AAAA", "AAAA" };
    var msa_data = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa_data.deinit();
    const avg = try averageIdentity(allocator, msa_data, 10000);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), avg, 1e-10);
}

test "averageIdentity: single sequence returns 1.0" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"AAAA"};
    var msa_data = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa_data.deinit();
    const avg = try averageIdentity(allocator, msa_data, 10000);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), avg, 1e-10);
}

test "averageIdentity: 3 seqs exhaustive matches manual calculation" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "AAAA", "AATT", "TTTT" };
    var msa_data = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa_data.deinit();
    // s1 vs s2: pid=0.5, s1 vs s3: pid=0.0, s2 vs s3: pid=0.5
    // average = (0.5 + 0.0 + 0.5) / 3 = 1/3
    const avg = try averageIdentity(allocator, msa_data, 10000);
    try std.testing.expectApproxEqAbs(1.0 / 3.0, avg, 1e-10);
}

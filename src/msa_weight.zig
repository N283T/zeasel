// MSA sequence weighting algorithms.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("msa.zig").Msa;

/// Compute position-based (PB/Henikoff & Henikoff 1994) weights for an MSA.
///
/// For each alignment column:
///   - Count r: the number of distinct canonical residue types
///     (gaps and degenerate codes are excluded).
///   - For each residue type, count n: how many sequences carry that residue.
///   - Each sequence at that column contributes 1 / (r * n) to its weight.
///
/// After accumulation, each weight is divided by the sequence's residue
/// count (per-sequence-length normalization, matching Easel's
/// esl_msaweight_PB_adv). Then weights are normalized to sum to nseq.
///
/// Returns an allocated slice of length nseq. Caller owns the memory.
pub fn positionBased(allocator: Allocator, m: Msa) ![]f64 {
    const n = m.nseq();
    const k = m.abc.k; // number of canonical residue codes
    const weights = try allocator.alloc(f64, n);
    @memset(weights, 0.0);

    // amino acid kp is 29 — the largest kp we support.
    const MAX_KP = 30;

    for (0..m.alen) |col| {
        // Count occurrences of each digital code in this column.
        // Only canonical residues (code < K) are counted; degenerate
        // codes (>= K), gaps, missing, and nonresidues are ignored.
        var type_counts: [MAX_KP]u32 = .{0} ** MAX_KP;
        for (0..n) |seq| {
            const code = m.seqs[seq][col];
            if (code < k) {
                type_counts[code] += 1;
            }
        }

        // Count how many distinct canonical types are present.
        var r: u32 = 0;
        for (type_counts[0..k]) |c| {
            if (c > 0) r += 1;
        }
        if (r == 0) continue; // all-gap/degenerate column

        // Accumulate weight contribution for each canonical-residue sequence.
        for (0..n) |seq| {
            const code = m.seqs[seq][col];
            if (code < k) {
                const n_i: f64 = @floatFromInt(type_counts[code]);
                const r_f: f64 = @floatFromInt(r);
                weights[seq] += 1.0 / (r_f * n_i);
            }
        }
    }

    // Per-sequence-length normalization: divide each weight by the number
    // of canonical residues in that sequence (matching Easel's behavior).
    for (0..n) |seq| {
        var rlen: u32 = 0;
        for (0..m.alen) |col| {
            if (m.seqs[seq][col] < k) rlen += 1;
        }
        if (rlen > 0) weights[seq] /= @floatFromInt(rlen);
    }

    // Normalize so that weights sum to nseq.
    var total: f64 = 0.0;
    for (weights) |w| total += w;
    if (total > 0.0) {
        const scale = @as(f64, @floatFromInt(n)) / total;
        for (weights) |*w| w.* *= scale;
    }

    return weights;
}

/// Compute BLOSUM-style clustering weights.
///
/// Sequences are clustered by single-linkage at the given percent-identity
/// threshold (e.g. 0.62 for 62%). Every sequence in a cluster is assigned
/// weight 1 / cluster_size. Weights are normalized to sum to nseq.
///
/// Returns an allocated slice of length nseq. Caller owns the memory.
pub fn blosum(allocator: Allocator, m: Msa, id_threshold: f64) ![]f64 {
    const n = m.nseq();

    // Union-find parent array: parent[i] == i means i is a root.
    var parent = try allocator.alloc(usize, n);
    defer allocator.free(parent);
    for (0..n) |i| parent[i] = i;

    // Single-linkage clustering: merge pairs that meet the threshold.
    for (0..n) |i| {
        for (i + 1..n) |j| {
            const pid = pairwiseIdentity(m, i, j);
            if (pid >= id_threshold) {
                unionSets(parent, i, j);
            }
        }
    }

    // Count how many sequences belong to each root's cluster.
    var cluster_size = try allocator.alloc(u32, n);
    defer allocator.free(cluster_size);
    @memset(cluster_size, 0);
    for (0..n) |i| {
        cluster_size[find(parent, i)] += 1;
    }

    // Each sequence gets weight 1 / cluster_size.
    const weights = try allocator.alloc(f64, n);
    for (0..n) |i| {
        weights[i] = 1.0 / @as(f64, @floatFromInt(cluster_size[find(parent, i)]));
    }

    // Normalize to sum to nseq.
    var total: f64 = 0.0;
    for (weights) |w| total += w;
    if (total > 0.0) {
        const scale = @as(f64, @floatFromInt(n)) / total;
        for (weights) |*w| w.* *= scale;
    }

    return weights;
}

// TODO: GSC (Gerstein-Sonnhammer-Chothia) weights require a guide tree.
// Depends on issue #13 (phylogenetic trees). Implement after that is available.

/// Compute percent identity between two aligned sequences.
/// Positions where either sequence has a gap are excluded from the comparison.
/// Returns 0.0 if there are no comparable (ungapped) positions.
/// Pairwise identity using Easel's definition: MIN(len1, len2) denominator.
pub fn pairwiseIdentity(m: Msa, i: usize, j: usize) f64 {
    var matches: u32 = 0;
    var len1: u32 = 0;
    var len2: u32 = 0;
    for (0..m.alen) |col| {
        const a = m.seqs[i][col];
        const b = m.seqs[j][col];
        const a_res = !m.abc.isGap(a);
        const b_res = !m.abc.isGap(b);
        if (a_res) len1 += 1;
        if (b_res) len2 += 1;
        if (a_res and b_res and a == b) matches += 1;
    }
    const denom = @min(len1, len2);
    if (denom == 0) return 0.0;
    return @as(f64, @floatFromInt(matches)) / @as(f64, @floatFromInt(denom));
}

// --- Union-Find helpers ---

fn find(parent: []usize, i: usize) usize {
    var x = i;
    while (parent[x] != x) x = parent[x];
    return x;
}

fn unionSets(parent: []usize, i: usize, j: usize) void {
    const ri = find(parent, i);
    const rj = find(parent, j);
    if (ri != rj) parent[ri] = rj;
}

// --- Tests ---

test "pairwiseIdentity: identical sequences return 1.0" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "ACGT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), pairwiseIdentity(m, 0, 1), 1e-9);
}

test "pairwiseIdentity: completely different sequences return 0.0" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "AAAA", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), pairwiseIdentity(m, 0, 1), 1e-9);
}

test "pairwiseIdentity: gaps are excluded from comparison" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    // Compared positions: 1 (A vs A match), 4 (T vs T match) = 2 matches / 2 compared
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "A-G-T", "A-C-T" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    // Ungapped columns: 0 (A=A match), 2 (G vs C mismatch), 4 (T=T match)
    // identity = 2/3
    try std.testing.expectApproxEqAbs(2.0 / 3.0, pairwiseIdentity(m, 0, 1), 1e-9);
}

test "positionBased: 3 identical sequences get equal weights" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    // All sequences identical -> each column has r=1, n=3 -> contribution = 1/(1*3)
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "ACGT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try positionBased(allocator, m);
    defer allocator.free(w);
    // After normalization: each weight = 1.0 (sum = 3 = nseq)
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[1], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[2], 1e-9);
}

test "positionBased: 2 completely different sequences each get weight 1.0" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    // Each column has r=2 residue types, n=1 for each -> contribution = 1/(2*1) = 0.5 per col
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "AAAA", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try positionBased(allocator, m);
    defer allocator.free(w);
    // Both sequences accumulate the same raw weight -> normalize to 1.0 each
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[1], 1e-9);
}

test "positionBased: known diversity example" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    // 3 seqs, length 4:
    //   s1: ACGT
    //   s2: ACGT  (same as s1)
    //   s3: TTTT  (different)
    //
    // Column 0 (A,A,T): r=2, n_A=2, n_T=1
    //   s1 += 1/(2*2) = 0.25,  s2 += 0.25,  s3 += 1/(2*1) = 0.5
    // Column 1 (C,C,T): r=2, n_C=2, n_T=1
    //   s1 += 0.25,             s2 += 0.25,  s3 += 0.5
    // Column 2 (G,G,T): r=2, n_G=2, n_T=1
    //   s1 += 0.25,             s2 += 0.25,  s3 += 0.5
    // Column 3 (T,T,T): r=1, n_T=3
    //   s1 += 1/(1*3) ≈ 0.333, s2 += 0.333, s3 += 0.333
    //
    // Raw: s1 = 0.25*3 + 0.333 = 1.083,  s2 = 1.083,  s3 = 0.5*3 + 0.333 = 1.833
    // Sum = 1.083 + 1.083 + 1.833 = 4.0
    // scale = 3/4 = 0.75
    // w[0] = 1.083 * 0.75 = 0.8125,  w[1] = 0.8125,  w[2] = 1.833 * 0.75 = 1.375
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try positionBased(allocator, m);
    defer allocator.free(w);

    // Recompute expected values precisely using fractions.
    // raw_s1 = raw_s2 = 3*(1/4) + 1/3 = 3/4 + 1/3 = 9/12 + 4/12 = 13/12
    // raw_s3 = 3*(1/2) + 1/3 = 3/2 + 1/3 = 9/6 + 2/6 = 11/6
    // sum = 13/12 + 13/12 + 11/6 = 13/12 + 13/12 + 22/12 = 48/12 = 4
    // scale = 3/4
    // w[0] = w[1] = (13/12) * (3/4) = 39/48 = 13/16 = 0.8125
    // w[2] = (11/6) * (3/4) = 33/24 = 11/8 = 1.375
    try std.testing.expectApproxEqAbs(13.0 / 16.0, w[0], 1e-9);
    try std.testing.expectApproxEqAbs(13.0 / 16.0, w[1], 1e-9);
    try std.testing.expectApproxEqAbs(11.0 / 8.0, w[2], 1e-9);
}

test "blosum: 3 identical sequences all in one cluster" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "ACGT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try blosum(allocator, m, 0.62);
    defer allocator.free(w);
    // All in one cluster of size 3: raw weight = 1/3 each.
    // Normalize: scale = 3 / (3 * 1/3) = 3/1 = 3 -> w[i] = (1/3)*3 = 1.0
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[1], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[2], 1e-9);
}

test "blosum: 2 identical + 1 different at threshold 0.62" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    // s1 and s2 are identical (pid=1.0 >= 0.62) -> one cluster of size 2.
    // s3 differs from both (pid=0.0 < 0.62) -> singleton cluster.
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try blosum(allocator, m, 0.62);
    defer allocator.free(w);
    // raw: s1=0.5, s2=0.5, s3=1.0  -> sum = 2.0
    // scale = 3/2
    // w[0] = w[1] = 0.5 * 1.5 = 0.75
    // w[2] = 1.0 * 1.5 = 1.5
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), w[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), w[1], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), w[2], 1e-9);
}

// MSA clustering and manipulation operations.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("msa.zig").Msa;
const Random = @import("util/random.zig").Random;
const pairwiseIdentity = @import("msa_weight.zig").pairwiseIdentity;

// --- Union-Find helpers (local copies — not exported from msa_weight) ---

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

/// Single-linkage clustering of MSA sequences by percent identity.
/// Returns cluster assignments: result[i] = cluster ID for sequence i.
/// Cluster IDs are 0-indexed and contiguous. Caller owns the returned slice.
pub fn clusterByIdentity(allocator: Allocator, m: Msa, id_threshold: f64) ![]usize {
    const n = m.nseq();

    // Union-find parent array.
    var parent = try allocator.alloc(usize, n);
    defer allocator.free(parent);
    for (0..n) |i| parent[i] = i;

    // Single-linkage: union any pair whose identity meets the threshold.
    for (0..n) |i| {
        for (i + 1..n) |j| {
            const pid = pairwiseIdentity(m, i, j);
            if (pid >= id_threshold) {
                unionSets(parent, i, j);
            }
        }
    }

    // Map each root to a sequential cluster ID.
    const cluster_id = try allocator.alloc(usize, n);
    errdefer allocator.free(cluster_id);

    // root_to_id: maps root index -> assigned cluster ID.
    // We reuse a temporary array of size n, initialized to max (sentinel = unassigned).
    const root_to_id = try allocator.alloc(usize, n);
    defer allocator.free(root_to_id);
    @memset(root_to_id, std.math.maxInt(usize));

    var next_id: usize = 0;
    for (0..n) |i| {
        const root = find(parent, i);
        if (root_to_id[root] == std.math.maxInt(usize)) {
            root_to_id[root] = next_id;
            next_id += 1;
        }
        cluster_id[i] = root_to_id[root];
    }

    return cluster_id;
}

/// Count residue composition across all sequences in an MSA.
/// Returns array of length abc.kp with raw counts for each digital code.
/// Caller owns the returned slice.
pub fn composition(allocator: Allocator, m: Msa) ![]f64 {
    const kp: usize = @intCast(m.abc.kp);
    const counts = try allocator.alloc(f64, kp);
    @memset(counts, 0.0);

    for (m.seqs) |seq| {
        for (seq) |code| {
            if (code < kp) {
                counts[code] += 1.0;
            }
        }
    }

    return counts;
}

/// Remove columns where fraction of gaps exceeds max_gap_fraction.
/// Returns a new Msa containing only the retained columns. Caller owns it.
pub fn removeGappyColumns(allocator: Allocator, m: Msa, max_gap_fraction: f64) !Msa {
    const n = m.nseq();
    const mask = try allocator.alloc(bool, m.alen);
    defer allocator.free(mask);

    for (0..m.alen) |col| {
        var gap_count: usize = 0;
        for (0..n) |seq| {
            if (m.abc.isGap(m.seqs[seq][col])) gap_count += 1;
        }
        const gap_frac: f64 = @as(f64, @floatFromInt(gap_count)) / @as(f64, @floatFromInt(n));
        mask[col] = gap_frac <= max_gap_fraction;
    }

    return m.selectColumns(mask);
}

/// Mark sequences as fragments if they have fewer than threshold fraction of
/// non-gap residues relative to the alignment length.
/// Returns a boolean slice: true = fragment. Caller owns it.
pub fn markFragments(allocator: Allocator, m: Msa, threshold: f64) ![]bool {
    const n = m.nseq();
    const flags = try allocator.alloc(bool, n);

    for (0..n) |i| {
        var non_gap: usize = 0;
        for (m.seqs[i]) |code| {
            if (!m.abc.isGap(code)) non_gap += 1;
        }
        const frac: f64 = @as(f64, @floatFromInt(non_gap)) / @as(f64, @floatFromInt(m.alen));
        flags[i] = frac < threshold;
    }

    return flags;
}

/// Shuffle columns of an MSA in place using Fisher-Yates (for null model testing).
pub fn shuffleColumns(rng: *Random, m: *Msa) void {
    if (m.alen <= 1) return;

    var i: usize = m.alen - 1;
    while (i > 0) : (i -= 1) {
        // j in [0, i]
        const j = rng.uniformInt(@intCast(i + 1));
        // Swap column i with column j across all rows.
        for (m.seqs) |seq| {
            const tmp = seq[i];
            seq[i] = seq[j];
            seq[j] = tmp;
        }
    }
}

// --- Tests ---

test "clusterByIdentity: 2 identical + 1 different -> 2 clusters" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    // s1 and s2 are identical (pid=1.0), s3 is completely different (pid=0.0).
    const seqs = [_][]const u8{ "ACGT", "ACGT", "TTTT" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    const clusters = try clusterByIdentity(allocator, m, 0.62);
    defer allocator.free(clusters);

    try std.testing.expectEqual(@as(usize, 3), clusters.len);
    // s1 and s2 must be in the same cluster.
    try std.testing.expectEqual(clusters[0], clusters[1]);
    // s3 must be in a different cluster.
    try std.testing.expect(clusters[0] != clusters[2]);
    // Exactly 2 distinct cluster IDs.
    const id0 = clusters[0];
    const id2 = clusters[2];
    const distinct = (id0 != id2);
    try std.testing.expect(distinct);
}

test "clusterByIdentity: all identical -> 1 cluster" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "ACGT" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    const clusters = try clusterByIdentity(allocator, m, 0.62);
    defer allocator.free(clusters);

    try std.testing.expectEqual(clusters[0], clusters[1]);
    try std.testing.expectEqual(clusters[1], clusters[2]);
}

test "clusterByIdentity: all different -> 3 clusters" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "AAAA", "CCCC", "GGGG" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    const clusters = try clusterByIdentity(allocator, m, 0.62);
    defer allocator.free(clusters);

    // All pairwise identities are 0.0 < 0.62, so each gets its own cluster.
    try std.testing.expect(clusters[0] != clusters[1]);
    try std.testing.expect(clusters[1] != clusters[2]);
    try std.testing.expect(clusters[0] != clusters[2]);
    // IDs should be 0, 1, 2 in order of first encounter.
    try std.testing.expectEqual(@as(usize, 0), clusters[0]);
    try std.testing.expectEqual(@as(usize, 1), clusters[1]);
    try std.testing.expectEqual(@as(usize, 2), clusters[2]);
}

test "composition: known MSA -> expected counts" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // Two sequences: ACGT and ACGT -> 2 A, 2 C, 2 G, 2 T, 0 gaps.
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "ACGT" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    const counts = try composition(allocator, m);
    defer allocator.free(counts);

    // DNA codes: A=0, C=1, G=2, T=3, gap=4
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), counts[0], 1e-9); // A
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), counts[1], 1e-9); // C
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), counts[2], 1e-9); // G
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), counts[3], 1e-9); // T
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), counts[4], 1e-9); // gap
}

test "composition: counts gaps correctly" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // A-GT: gap code is 4.
    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"A-GT"};

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    const counts = try composition(allocator, m);
    defer allocator.free(counts);

    try std.testing.expectApproxEqAbs(@as(f64, 1.0), counts[0], 1e-9); // A
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), counts[4], 1e-9); // gap
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), counts[2], 1e-9); // G
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), counts[3], 1e-9); // T
}

test "removeGappyColumns: gap-heavy column is removed" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // 4 seqs, column 2 is all-gap.
    const names = [_][]const u8{ "s1", "s2", "s3", "s4" };
    const seqs = [_][]const u8{ "AC-T", "AC-T", "AC-T", "AC-T" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    // max_gap_fraction=0.5 -> column with 100% gaps should be removed.
    var trimmed = try removeGappyColumns(allocator, m, 0.5);
    defer trimmed.deinit();

    try std.testing.expectEqual(@as(usize, 3), trimmed.alen);

    const text0 = try trimmed.abc.textize(allocator, trimmed.seqs[0]);
    defer allocator.free(text0);
    try std.testing.expectEqualStrings("ACT", text0);
}

test "removeGappyColumns: column below threshold is kept" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // 4 seqs, column 2 has 1 gap out of 4 = 25% gaps.
    const names = [_][]const u8{ "s1", "s2", "s3", "s4" };
    const seqs = [_][]const u8{ "ACT", "ACT", "ACT", "A-T" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    // max_gap_fraction=0.5 -> 25% < 50%, column should be kept.
    var trimmed = try removeGappyColumns(allocator, m, 0.5);
    defer trimmed.deinit();

    try std.testing.expectEqual(@as(usize, 3), trimmed.alen);
}

test "markFragments: short sequence marked as fragment" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // s1 has 1 non-gap residue out of 4 columns = 25% coverage.
    // s2 has 4 non-gap residues = 100% coverage.
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "---A", "ACGT" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    // threshold=0.5: s1 (0.25) is a fragment, s2 (1.0) is not.
    const flags = try markFragments(allocator, m, 0.5);
    defer allocator.free(flags);

    try std.testing.expect(flags[0]); // s1 is a fragment
    try std.testing.expect(!flags[1]); // s2 is not
}

test "markFragments: no fragments when all sequences are full length" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    const flags = try markFragments(allocator, m, 0.5);
    defer allocator.free(flags);

    try std.testing.expect(!flags[0]);
    try std.testing.expect(!flags[1]);
}

test "shuffleColumns: composition preserved after shuffle" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    // Record composition before shuffle.
    const before = try composition(allocator, m);
    defer allocator.free(before);

    var rng = Random.init(42);
    shuffleColumns(&rng, &m);

    // Composition should be identical after shuffle (only order changes).
    const after = try composition(allocator, m);
    defer allocator.free(after);

    for (before, after) |b, a| {
        try std.testing.expectApproxEqAbs(b, a, 1e-9);
    }

    // Alignment length should be unchanged.
    try std.testing.expectEqual(@as(usize, 4), m.alen);
}

test "shuffleColumns: single column MSA is unchanged" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"A"};

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    var rng = Random.init(1);
    shuffleColumns(&rng, &m);

    try std.testing.expectEqual(@as(usize, 1), m.alen);
    // A = code 0
    try std.testing.expectEqual(@as(u8, 0), m.seqs[0][0]);
}

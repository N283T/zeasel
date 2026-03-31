// MSA clustering and manipulation operations.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("alphabet.zig").Alphabet;
const Msa = @import("msa.zig").Msa;
const Random = @import("util/random.zig").Random;
const pairwiseIdentity = @import("msa_weight.zig").pairwiseIdentity;
const wuss = @import("wuss.zig");

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

/// Bootstrap resampling: resample columns WITH replacement to create a
/// bootstrap alignment sample. The result has the same number of columns
/// (alen) but columns are randomly selected with replacement.
/// Returns a new Msa; caller owns it.
pub fn bootstrap(allocator: Allocator, m: Msa, rng: *Random) !Msa {
    const n = m.nseq();
    const alen = m.alen;
    if (alen == 0) return error.InvalidInput;

    // Generate column indices by sampling with replacement.
    const col_indices = try allocator.alloc(usize, alen);
    defer allocator.free(col_indices);
    for (0..alen) |i| {
        col_indices[i] = rng.uniformInt(@intCast(alen));
    }

    // Build new names.
    var new_names = try allocator.alloc([]const u8, n);
    var names_done: usize = 0;
    errdefer {
        for (0..names_done) |i| allocator.free(new_names[i]);
        allocator.free(new_names);
    }

    // Build new sequences.
    var new_seqs = try allocator.alloc([]u8, n);
    var seqs_done: usize = 0;
    errdefer {
        for (0..seqs_done) |i| allocator.free(new_seqs[i]);
        allocator.free(new_seqs);
    }

    for (0..n) |i| {
        new_names[i] = try allocator.dupe(u8, m.names[i]);
        names_done += 1;

        new_seqs[i] = try allocator.alloc(u8, alen);
        seqs_done += 1;
        for (0..alen) |c| {
            new_seqs[i][c] = m.seqs[i][col_indices[c]];
        }
    }

    // Copy per-column annotations using the same resampling.
    var new_ss: ?[]const u8 = null;
    errdefer if (new_ss) |s| allocator.free(s);
    if (m.consensus_ss) |ss| {
        const buf = try allocator.alloc(u8, alen);
        for (0..alen) |c| buf[c] = ss[col_indices[c]];
        new_ss = buf;
    }

    var new_rf: ?[]const u8 = null;
    errdefer if (new_rf) |s| allocator.free(s);
    if (m.reference) |rf| {
        const buf = try allocator.alloc(u8, alen);
        for (0..alen) |c| buf[c] = rf[col_indices[c]];
        new_rf = buf;
    }

    var new_weights: ?[]f64 = null;
    errdefer if (new_weights) |w| allocator.free(w);
    if (m.weights) |w| {
        new_weights = try allocator.dupe(f64, w);
    }

    return Msa{
        .names = new_names,
        .seqs = new_seqs,
        .alen = alen,
        .abc = m.abc,
        .allocator = allocator,
        .name = if (m.name) |v| try allocator.dupe(u8, v) else null,
        .accession = if (m.accession) |v| try allocator.dupe(u8, v) else null,
        .description = if (m.description) |v| try allocator.dupe(u8, v) else null,
        .author = if (m.author) |v| try allocator.dupe(u8, v) else null,
        .weights = new_weights,
        .consensus_ss = new_ss,
        .reference = new_rf,
    };
}

/// Vertical shuffle: shuffle residues within each column independently,
/// preserving per-column composition but destroying sequence identity.
/// Returns a new Msa; caller owns it.
pub fn shuffleVertical(allocator: Allocator, m: Msa, rng: *Random) !Msa {
    const n = m.nseq();
    const alen = m.alen;

    // Build new names.
    var new_names = try allocator.alloc([]const u8, n);
    var names_done: usize = 0;
    errdefer {
        for (0..names_done) |i| allocator.free(new_names[i]);
        allocator.free(new_names);
    }

    // Deep-copy sequences (we will shuffle in place per-column).
    var new_seqs = try allocator.alloc([]u8, n);
    var seqs_done: usize = 0;
    errdefer {
        for (0..seqs_done) |i| allocator.free(new_seqs[i]);
        allocator.free(new_seqs);
    }

    for (0..n) |i| {
        new_names[i] = try allocator.dupe(u8, m.names[i]);
        names_done += 1;
        new_seqs[i] = try allocator.dupe(u8, m.seqs[i]);
        seqs_done += 1;
    }

    // Fisher-Yates shuffle within each column.
    for (0..alen) |col| {
        if (n <= 1) continue;
        var i: usize = n - 1;
        while (i > 0) : (i -= 1) {
            const j = rng.uniformInt(@intCast(i + 1));
            const tmp = new_seqs[i][col];
            new_seqs[i][col] = new_seqs[j][col];
            new_seqs[j][col] = tmp;
        }
    }

    // Copy per-column annotations unchanged.
    var new_ss: ?[]const u8 = null;
    errdefer if (new_ss) |s| allocator.free(s);
    if (m.consensus_ss) |ss| {
        new_ss = try allocator.dupe(u8, ss);
    }

    var new_rf: ?[]const u8 = null;
    errdefer if (new_rf) |s| allocator.free(s);
    if (m.reference) |rf| {
        new_rf = try allocator.dupe(u8, rf);
    }

    var new_weights: ?[]f64 = null;
    errdefer if (new_weights) |w| allocator.free(w);
    if (m.weights) |w| {
        new_weights = try allocator.dupe(f64, w);
    }

    return Msa{
        .names = new_names,
        .seqs = new_seqs,
        .alen = alen,
        .abc = m.abc,
        .allocator = allocator,
        .name = if (m.name) |v| try allocator.dupe(u8, v) else null,
        .accession = if (m.accession) |v| try allocator.dupe(u8, v) else null,
        .description = if (m.description) |v| try allocator.dupe(u8, v) else null,
        .author = if (m.author) |v| try allocator.dupe(u8, v) else null,
        .weights = new_weights,
        .consensus_ss = new_ss,
        .reference = new_rf,
    };
}

/// Flush-left insert regions: normalize insert columns so gaps follow
/// residues (flush-left). Consensus columns (identified by RF annotation
/// where the RF character is not '.' or '~') are left unchanged.
/// If no RF annotation is present, all columns are treated as consensus
/// and the MSA is returned unchanged (as a copy).
/// Returns a new Msa; caller owns it.
pub fn flushLeftInserts(allocator: Allocator, m: Msa) !Msa {
    const n = m.nseq();
    const alen = m.alen;

    // Determine which columns are insert columns from RF annotation.
    const rf = m.reference;

    // Deep-copy names.
    var new_names = try allocator.alloc([]const u8, n);
    var names_done: usize = 0;
    errdefer {
        for (0..names_done) |i| allocator.free(new_names[i]);
        allocator.free(new_names);
    }

    // Deep-copy sequences.
    var new_seqs = try allocator.alloc([]u8, n);
    var seqs_done: usize = 0;
    errdefer {
        for (0..seqs_done) |i| allocator.free(new_seqs[i]);
        allocator.free(new_seqs);
    }

    for (0..n) |i| {
        new_names[i] = try allocator.dupe(u8, m.names[i]);
        names_done += 1;
        new_seqs[i] = try allocator.dupe(u8, m.seqs[i]);
        seqs_done += 1;
    }

    // If RF annotation exists, flush-left each insert block.
    if (rf) |rf_ann| {
        var col: usize = 0;
        while (col < alen) {
            // Skip consensus columns.
            if (!isInsertColumn(rf_ann[col])) {
                col += 1;
                continue;
            }
            // Found start of an insert block.
            const block_start = col;
            while (col < alen and isInsertColumn(rf_ann[col])) {
                col += 1;
            }
            const block_end = col;
            const block_len = block_end - block_start;

            // For each sequence, flush residues left within this block.
            for (0..n) |seq_idx| {
                flushLeftBlock(new_seqs[seq_idx], block_start, block_len, m.abc);
            }
        }
    }

    // Copy per-column annotations.
    var new_ss: ?[]const u8 = null;
    errdefer if (new_ss) |s| allocator.free(s);
    if (m.consensus_ss) |ss| {
        new_ss = try allocator.dupe(u8, ss);
    }

    var new_rf: ?[]const u8 = null;
    errdefer if (new_rf) |s| allocator.free(s);
    if (rf) |rf_ann| {
        new_rf = try allocator.dupe(u8, rf_ann);
    }

    var new_weights: ?[]f64 = null;
    errdefer if (new_weights) |w| allocator.free(w);
    if (m.weights) |w| {
        new_weights = try allocator.dupe(f64, w);
    }

    return Msa{
        .names = new_names,
        .seqs = new_seqs,
        .alen = alen,
        .abc = m.abc,
        .allocator = allocator,
        .name = if (m.name) |v| try allocator.dupe(u8, v) else null,
        .accession = if (m.accession) |v| try allocator.dupe(u8, v) else null,
        .description = if (m.description) |v| try allocator.dupe(u8, v) else null,
        .author = if (m.author) |v| try allocator.dupe(u8, v) else null,
        .weights = new_weights,
        .consensus_ss = new_ss,
        .reference = new_rf,
    };
}

/// Returns true if the RF character marks an insert column.
fn isInsertColumn(rf_char: u8) bool {
    return rf_char == '.' or rf_char == '~';
}

/// Flush residues to the left within a block of an aligned sequence,
/// moving all gap characters to the right end of the block.
fn flushLeftBlock(seq: []u8, start: usize, len: usize, abc: *const Alphabet) void {
    var write_pos: usize = 0;
    for (0..len) |i| {
        if (!abc.isGap(seq[start + i])) {
            seq[start + write_pos] = seq[start + i];
            write_pos += 1;
        }
    }
    const gap = abc.gapCode();
    while (write_pos < len) : (write_pos += 1) {
        seq[start + write_pos] = gap;
    }
}

/// Remove broken base pairs from a WUSS secondary structure annotation.
/// Given a WUSS string and a column-keep mask, find all base pairs using
/// `wuss.parseToPairs()`, then replace characters of broken pairs (where
/// one partner is removed) with '.'.
/// Returns a new string; caller owns it.
pub fn removeBrokenBasepairs(allocator: Allocator, ss: []const u8, keep: []const bool) ![]u8 {
    std.debug.assert(ss.len == keep.len);

    const pairs = try wuss.parseToPairs(allocator, ss);
    defer allocator.free(pairs);

    const result = try allocator.alloc(u8, ss.len);
    @memcpy(result, ss);

    for (0..ss.len) |i| {
        if (pairs[i] >= 0) {
            const partner: usize = @intCast(pairs[i]);
            // If either partner is removed, break both sides.
            if (!keep[i] or !keep[partner]) {
                result[i] = '.';
                result[partner] = '.';
            }
        }
    }

    return result;
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

test "bootstrap: result has same dimensions as original" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AACC" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    var rng = Random.init(99);
    var bs = try bootstrap(allocator, m, &rng);
    defer bs.deinit();

    try std.testing.expectEqual(m.nseq(), bs.nseq());
    try std.testing.expectEqual(m.alen, bs.alen);
}

test "bootstrap: each column in result comes from original" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    var rng = Random.init(7);
    var bs = try bootstrap(allocator, m, &rng);
    defer bs.deinit();

    // Each column in bs must match some column in m.
    for (0..bs.alen) |col| {
        var found = false;
        for (0..m.alen) |orig_col| {
            var matches = true;
            for (0..m.nseq()) |seq| {
                if (bs.seqs[seq][col] != m.seqs[seq][orig_col]) {
                    matches = false;
                    break;
                }
            }
            if (matches) {
                found = true;
                break;
            }
        }
        try std.testing.expect(found);
    }
}

test "shuffleVertical: preserves per-column composition" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2", "s3", "s4" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AACC", "GGTT" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    // Record per-column composition before shuffle.
    var col_comp_before: [4][6]usize = undefined;
    for (0..m.alen) |col| {
        @memset(&col_comp_before[col], 0);
        for (0..m.nseq()) |seq| {
            col_comp_before[col][m.seqs[seq][col]] += 1;
        }
    }

    var rng = Random.init(42);
    var shuffled = try shuffleVertical(allocator, m, &rng);
    defer shuffled.deinit();

    // Per-column composition must be identical.
    for (0..shuffled.alen) |col| {
        var col_comp_after: [6]usize = undefined;
        @memset(&col_comp_after, 0);
        for (0..shuffled.nseq()) |seq| {
            col_comp_after[shuffled.seqs[seq][col]] += 1;
        }
        for (0..6) |code| {
            try std.testing.expectEqual(col_comp_before[col][code], col_comp_after[code]);
        }
    }
}

test "shuffleVertical: result has same dimensions" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    var rng = Random.init(1);
    var shuffled = try shuffleVertical(allocator, m, &rng);
    defer shuffled.deinit();

    try std.testing.expectEqual(m.nseq(), shuffled.nseq());
    try std.testing.expectEqual(m.alen, shuffled.alen);
}

test "flushLeftInserts: gaps move right within insert blocks" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // 6 columns: RF = "x..x.x"
    // Columns 0, 3, 5 are consensus (x), columns 1-2 and 4 are insert (.).
    // s1: A--GCT  -> insert block [1,2]: "--" stays "--"; block [4]: "C" stays "C"
    // s2: A-CGAT  -> insert block [1,2]: "-C" -> "C-"; block [4]: "A" stays "A"
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "A--GCT", "A-CGAT" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    // Set RF annotation.
    const rf_copy = try allocator.dupe(u8, "x..x.x");
    m.reference = rf_copy;

    var flushed = try flushLeftInserts(allocator, m);
    defer flushed.deinit();

    // s1 insert block [1,2] was "--" -> stays "--" (all gaps)
    const text0 = try flushed.abc.textize(allocator, flushed.seqs[0]);
    defer allocator.free(text0);
    try std.testing.expectEqualStrings("A--GCT", text0);

    // s2 insert block [1,2] was "-C" -> "C-" (residue flushed left)
    const text1 = try flushed.abc.textize(allocator, flushed.seqs[1]);
    defer allocator.free(text1);
    try std.testing.expectEqualStrings("AC-GAT", text1);
}

test "flushLeftInserts: no RF annotation returns unchanged copy" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "A-GT", "ACGT" };

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    var flushed = try flushLeftInserts(allocator, m);
    defer flushed.deinit();

    // Without RF, no reordering should happen.
    const text0 = try flushed.abc.textize(allocator, flushed.seqs[0]);
    defer allocator.free(text0);
    try std.testing.expectEqualStrings("A-GT", text0);

    const text1 = try flushed.abc.textize(allocator, flushed.seqs[1]);
    defer allocator.free(text1);
    try std.testing.expectEqualStrings("ACGT", text1);
}

test "flushLeftInserts: multiple residues in insert block" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // RF = "x....x" -> columns 1-4 are insert.
    // s1: A--CGT -> insert block [1,2,3,4]: "--CG" -> "CG--"
    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"A--CGT"};

    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    const rf_copy = try allocator.dupe(u8, "x....x");
    m.reference = rf_copy;

    var flushed = try flushLeftInserts(allocator, m);
    defer flushed.deinit();

    const text = try flushed.abc.textize(allocator, flushed.seqs[0]);
    defer allocator.free(text);
    try std.testing.expectEqualStrings("ACG--T", text);
}

test "removeBrokenBasepairs: broken pair replaced with dots" {
    const allocator = std.testing.allocator;
    // SS: <<..>> means positions 0-1 paired with 5-4.
    // Keep mask: remove column 0 -> pair (0,5) is broken.
    const ss = "<<..>>";
    const keep = [_]bool{ false, true, true, true, true, true };

    const result = try removeBrokenBasepairs(allocator, ss, &keep);
    defer allocator.free(result);

    // Position 0 and 5 should be replaced with '.'.
    // Position 1 and 4 remain paired.
    try std.testing.expectEqualStrings(".<..>.", result);
}

test "removeBrokenBasepairs: all pairs intact" {
    const allocator = std.testing.allocator;
    const ss = "<<..>>";
    const keep = [_]bool{ true, true, true, true, true, true };

    const result = try removeBrokenBasepairs(allocator, ss, &keep);
    defer allocator.free(result);

    try std.testing.expectEqualStrings("<<..>>", result);
}

test "removeBrokenBasepairs: no pairs" {
    const allocator = std.testing.allocator;
    const ss = "......";
    const keep = [_]bool{ true, false, true, false, true, true };

    const result = try removeBrokenBasepairs(allocator, ss, &keep);
    defer allocator.free(result);

    try std.testing.expectEqualStrings("......", result);
}

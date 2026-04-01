// MSA sequence weighting algorithms.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("msa.zig").Msa;
const distance = @import("distance.zig");
const tree_mod = @import("tree.zig");

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
/// When `symfrac` is non-null, only consensus columns are used:
///   - If `m.reference` exists, consensus columns are those marked 'x' in it.
///   - Otherwise, `m.reasonableRF(symfrac)` is called to compute them.
///   - Non-consensus columns (marked '.') are skipped entirely, including
///     for the per-sequence-length (rlen) normalization.
///
/// Returns an allocated slice of length nseq. Caller owns the memory.
pub fn positionBased(allocator: Allocator, m: Msa, symfrac: ?f64) ![]f64 {
    const n = m.nseq();
    const k = m.abc.k; // number of canonical residue codes
    const weights = try allocator.alloc(f64, n);
    @memset(weights, 0.0);

    // Determine RF mask: null means use all columns.
    // When symfrac is provided and m.reference exists, use it directly.
    // Otherwise compute via reasonableRF. The computed slice is owned by
    // m.allocator and must be freed after use.
    var computed_rf: ?[]u8 = null;
    defer if (computed_rf) |rf| m.allocator.free(rf);

    const rf: ?[]const u8 = blk: {
        if (symfrac) |sf| {
            if (m.reference) |ref| {
                break :blk ref;
            } else {
                computed_rf = try m.reasonableRF(sf);
                break :blk computed_rf.?;
            }
        }
        break :blk null;
    };

    // amino acid kp is 29 — the largest kp we support.
    const MAX_KP = 30;

    for (0..m.alen) |col| {
        // Skip non-consensus columns when RF filtering is active.
        if (rf) |mask| {
            if (mask[col] == '.') continue;
        }

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
    // of canonical residues in that sequence in consensus columns only
    // (matching Easel's behavior; respects the same RF mask as above).
    for (0..n) |seq| {
        var rlen: u32 = 0;
        for (0..m.alen) |col| {
            if (rf) |mask| {
                if (mask[col] == '.') continue;
            }
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

/// Filter sequences by percent identity using a greedy algorithm.
/// Returns a boolean mask where `keep[i] = true` means sequence i is retained.
/// Sequences are considered in input order; the first sequence is always kept,
/// and subsequent sequences are removed if they exceed `max_id` identity with
/// any already-kept sequence.
///
/// When `prefer_conscover` is true, the sequence with fewer non-gap residues
/// in consensus columns (as marked by `msa.reference`) is removed instead of
/// always removing the later sequence. If no reference annotation is available,
/// falls back to removing the later sequence (same as prefer_conscover=false).
///
/// Caller owns the returned slice.
pub fn idFilter(allocator: Allocator, m: Msa, max_id: f64, prefer_conscover: bool) ![]bool {
    const n = m.nseq();
    const keep = try allocator.alloc(bool, n);
    @memset(keep, true);

    // Precompute consensus coverage per sequence if needed.
    const conscover: ?[]u32 = if (prefer_conscover and m.reference != null) blk: {
        const cc = try allocator.alloc(u32, n);
        const ref = m.reference.?;
        for (0..n) |s| {
            var count: u32 = 0;
            for (0..m.alen) |col| {
                if (ref[col] == '.' or ref[col] == '~') continue; // insert column
                if (m.abc.isResidue(m.seqs[s][col])) count += 1;
            }
            cc[s] = count;
        }
        break :blk cc;
    } else null;
    defer if (conscover) |cc| allocator.free(cc);

    for (0..n) |i| {
        if (!keep[i]) continue;
        for (i + 1..n) |j| {
            if (!keep[j]) continue;
            if (pairwiseIdentity(m, i, j) >= max_id) {
                // Decide which to remove.
                if (conscover) |cc| {
                    // Remove the one with lower consensus coverage.
                    if (cc[j] < cc[i]) {
                        keep[j] = false;
                    } else if (cc[i] < cc[j]) {
                        // Swap: remove i, but since i is the outer loop seq
                        // and we need to stop comparing from i, mark i as removed
                        // and break out of the inner loop.
                        keep[i] = false;
                        break;
                    } else {
                        // Tie: remove later sequence (default behavior).
                        keep[j] = false;
                    }
                } else {
                    keep[j] = false;
                }
            }
        }
    }
    return keep;
}

/// Compute GSC (Gerstein-Sonnhammer-Chothia 1994) sequence weights.
///
/// Algorithm:
///   1. Compute a pairwise fractional-difference distance matrix (1 - pid).
///   2. Build a UPGMA guide tree from the distances.
///   3. Distribute branch lengths from root to leaves:
///      - Postorder: compute total branch length under each internal node.
///      - Preorder: at each internal node, split the weight from above in
///        proportion to left/right subtree branch lengths. When both subtree
///        lengths are zero, split in proportion to clade sizes instead.
///      - Each leaf's weight = its accumulated share + its own branch length.
///   4. Normalize weights to sum to nseq.
///
/// Complexity: O(N^2 * L) for the distance matrix, O(N^3) for UPGMA,
/// O(N) for the weight distribution. Memory: O(N^2) for the distance matrix.
///
/// Reference: Gerstein, Sonnhammer & Chothia, "A method to weight protein
/// sequences to correct for unequal representation", JMB 236:1067-1078, 1994.
///
/// Returns an allocated slice of length nseq. Caller owns the memory.
pub fn gsc(allocator: Allocator, m: Msa) ![]f64 {
    const n = m.nseq();
    if (n == 1) {
        const weights = try allocator.alloc(f64, 1);
        weights[0] = 1.0;
        return weights;
    }

    // Step 1: Compute pairwise fractional-difference distance matrix.
    // GSC uses raw fractional difference (1 - pid), not JC-corrected distances.
    const dist_matrix = try distance.pairwiseDistanceMatrix(allocator, m, .none);
    defer allocator.free(dist_matrix);

    // Step 2: Build UPGMA guide tree.
    var t = try tree_mod.upgma(allocator, dist_matrix, n, m.names);
    defer t.deinit();

    const total_nodes = t.n_nodes;
    const n_leaves = t.n_leaves;

    // Step 3: Compute clade sizes for each node.
    const clade_size = try allocator.alloc(u32, total_nodes);
    defer allocator.free(clade_size);
    for (0..n_leaves) |i| clade_size[i] = 1;

    // Compute clade sizes for internal nodes via postorder traversal.
    // Internal nodes are indexed n_leaves..total_nodes-1. In the UPGMA tree,
    // internal nodes are created in order, so node i always has children < i.
    // Thus iterating from n_leaves upward is a valid postorder.
    for (n_leaves..total_nodes) |i| {
        const lc: usize = @intCast(t.left[i]);
        const rc: usize = @intCast(t.right[i]);
        clade_size[i] = clade_size[lc] + clade_size[rc];
    }

    // x[i] stores different things in postorder vs preorder passes.
    const x = try allocator.alloc(f64, total_nodes);
    defer allocator.free(x);
    @memset(x, 0.0);

    // Postorder: compute total branch length under each internal node.
    for (n_leaves..total_nodes) |i| {
        const lc: usize = @intCast(t.left[i]);
        const rc: usize = @intCast(t.right[i]);
        x[i] = t.branch_length[lc] + t.branch_length[rc];
        if (lc >= n_leaves) x[i] += x[lc];
        if (rc >= n_leaves) x[i] += x[rc];
    }

    // Preorder: distribute weight from root to leaves.
    // x[i] is now repurposed to mean "weight accumulated above node i".
    const root = total_nodes - 1;
    x[root] = 0.0; // no branch above the root

    const weights = try allocator.alloc(f64, n);

    // Process internal nodes from root downward. Since internal node indices
    // increase from n_leaves to total_nodes-1, and parent index > child index,
    // iterating in reverse is a valid preorder.
    var i_signed: i64 = @intCast(root);
    while (i_signed >= @as(i64, @intCast(n_leaves))) : (i_signed -= 1) {
        const i: usize = @intCast(i_signed);
        const lc: usize = @intCast(t.left[i]);
        const rc: usize = @intCast(t.right[i]);

        // Total branch weight in left and right subtrees.
        var lw: f64 = t.branch_length[lc];
        if (lc >= n_leaves) lw += x[lc]; // x[lc] still holds subtree length from postorder
        var rw: f64 = t.branch_length[rc];
        if (rc >= n_leaves) rw += x[rc];

        var lx: f64 = undefined;
        var rx: f64 = undefined;

        if (lw + rw == 0.0) {
            // Special case: all branch lengths zero in subtree.
            // Split in proportion to clade sizes.
            const cs_i: f64 = @floatFromInt(clade_size[i]);
            if (lc >= n_leaves) {
                lx = x[i] * @as(f64, @floatFromInt(clade_size[lc])) / cs_i;
            } else {
                lx = x[i] / cs_i;
            }
            if (rc >= n_leaves) {
                rx = x[i] * @as(f64, @floatFromInt(clade_size[rc])) / cs_i;
            } else {
                rx = x[i] / cs_i;
            }
        } else {
            // Normal case: split in proportion to branch weight.
            lx = x[i] * lw / (lw + rw);
            rx = x[i] * rw / (lw + rw);
        }

        // Assign to children: if leaf, record final weight; if internal, store for next iteration.
        if (lc < n_leaves) {
            weights[lc] = lx + t.branch_length[lc];
        } else {
            x[lc] = lx + t.branch_length[lc];
        }

        if (rc < n_leaves) {
            weights[rc] = rx + t.branch_length[rc];
        } else {
            x[rc] = rx + t.branch_length[rc];
        }
    }

    // Step 4: Normalize weights to sum to nseq.
    // When all branch lengths are zero (e.g., all sequences identical),
    // all raw weights are zero; fall back to equal weights.
    var total: f64 = 0.0;
    for (weights) |w| total += w;
    if (total > 0.0) {
        const scale = @as(f64, @floatFromInt(n)) / total;
        for (weights) |*w| w.* *= scale;
    } else {
        for (weights) |*w| w.* = 1.0;
    }

    return weights;
}

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
        const a_res = m.abc.isResidue(a);
        const b_res = m.abc.isResidue(b);
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
    const w = try positionBased(allocator, m, null);
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
    const w = try positionBased(allocator, m, null);
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
    const w = try positionBased(allocator, m, null);
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

test "positionBased: symfrac filters insert columns" {
    // MSA with 3 sequences, 4 columns.
    //   s1: AACG
    //   s2: AACG
    //   s3: A-CG  (column 1 is a gap / insert for s3)
    //
    // RF: x.xx  (column 1 is insert; columns 0, 2, 3 are consensus)
    //
    // Without symfrac (null): all 4 columns used.
    //   col 0: A,A,A -> r=1, n=3 -> each gets 1/3
    //   col 1: A,A,- -> r=1(A), n=2 -> s1+=1/2, s2+=1/2, s3 skipped
    //   col 2: C,C,C -> r=1, n=3 -> each gets 1/3
    //   col 3: G,G,G -> r=1, n=3 -> each gets 1/3
    //   raw: s1=s2=1/3+1/2+1/3+1/3=3/2, s3=1/3+1/3+1/3=1
    //   rlen: s1=4 (A,A,C,G), s2=4, s3=3 (A,C,G)
    //   after rlen: s1=s2=(3/2)/4=3/8, s3=1/3
    //   sum=3/8+3/8+1/3=9/24+9/24+8/24=26/24=13/12
    //   scale=3/(13/12)=36/13
    //   w[0]=w[1]=(3/8)*(36/13)=27/26, w[2]=(1/3)*(36/13)=12/13
    //
    // With symfrac=0.6 (RF x.xx -> only cols 0, 2, 3 used):
    //   col 0: A,A,A -> r=1, n=3 -> each gets 1/3
    //   col 2: C,C,C -> r=1, n=3 -> each gets 1/3
    //   col 3: G,G,G -> r=1, n=3 -> each gets 1/3
    //   raw: all = 1
    //   rlen (consensus cols only): s1=3, s2=3, s3=3
    //   after rlen: all = 1/3 -> equal -> normalize to 1.0 each
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "AACG", "AACG", "A-CG" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    // Attach RF annotation: column 1 is insert, the rest are consensus.
    m.reference = try m.allocator.dupe(u8, "x.xx");

    // Without symfrac: the insert column boosts s1/s2 raw weight without a
    // proportional rlen increase relative to s3, so s1/s2 end up heavier.
    const w_all = try positionBased(allocator, m, null);
    defer allocator.free(w_all);

    // With symfrac=0.6 and the explicit RF above: col 1 is skipped.
    const w_sf = try positionBased(allocator, m, 0.6);
    defer allocator.free(w_sf);

    // symfrac path: all three sequences see identical consensus columns -> equal weights.
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w_sf[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w_sf[1], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w_sf[2], 1e-9);

    // No-symfrac path: s1/s2 get higher weight than s3 (insert col boosts them).
    try std.testing.expect(w_all[0] > w_all[2]);

    // Exact no-symfrac values: w[0]=w[1]=27/26, w[2]=12/13.
    try std.testing.expectApproxEqAbs(27.0 / 26.0, w_all[0], 1e-9);
    try std.testing.expectApproxEqAbs(27.0 / 26.0, w_all[1], 1e-9);
    try std.testing.expectApproxEqAbs(12.0 / 13.0, w_all[2], 1e-9);
}

test "positionBased: symfrac computes RF when reference is null" {
    // Same MSA but without a pre-set RF: reasonableRF should compute it.
    // Columns 0, 2, 3: occupancy = 3/3 = 1.0 >= symfrac -> 'x'
    // Column 1: 2/3 occupancy < symfrac=0.7 -> '.'
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "AACG", "AACG", "A-CG" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();

    // No reference set — reasonableRF will be called.
    const w_sf = try positionBased(allocator, m, 0.7);
    defer allocator.free(w_sf);

    // With only cols 0, 2, 3: all sequences see the same residues -> 1.0 each.
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w_sf[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w_sf[1], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w_sf[2], 1e-9);
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

test "gsc: single sequence returns weight 1.0" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGT"};
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try gsc(allocator, m);
    defer allocator.free(w);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[0], 1e-9);
}

test "gsc: 2 identical sequences get equal weights" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "ACGT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try gsc(allocator, m);
    defer allocator.free(w);
    // Identical sequences should get equal weights, each = 1.0 (sum = nseq = 2)
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[1], 1e-9);
}

test "gsc: 2 completely different sequences get equal weights" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "AAAA", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try gsc(allocator, m);
    defer allocator.free(w);
    // Two sequences always get equal weight after normalization
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[1], 1e-9);
}

test "gsc: weights sum to nseq" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2", "s3", "s4" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "TTTT", "AATT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try gsc(allocator, m);
    defer allocator.free(w);
    var total: f64 = 0.0;
    for (w) |wt| total += wt;
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), total, 1e-9);
}

test "gsc: all weights non-negative" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2", "s3", "s4" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "TTTT", "AATT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try gsc(allocator, m);
    defer allocator.free(w);
    for (w) |wt| {
        try std.testing.expect(wt >= 0.0);
    }
}

test "gsc: unique sequence gets higher weight than duplicated ones" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    // s1 and s2 are identical (close on the tree), s3 is unique.
    // GSC should upweight s3 relative to s1 and s2.
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try gsc(allocator, m);
    defer allocator.free(w);
    // s1 and s2 should have equal weights (identical sequences)
    try std.testing.expectApproxEqAbs(w[0], w[1], 1e-9);
    // s3 should have higher weight than s1 or s2
    try std.testing.expect(w[2] > w[0]);
}

test "gsc: 3 identical sequences get equal weights" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "ACGT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const w = try gsc(allocator, m);
    defer allocator.free(w);
    // All identical -> all equal weight, each = 1.0
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[1], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), w[2], 1e-9);
}

test "idFilter: basic filtering removes similar sequences" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    // s1 and s2 are identical (pid=1.0), s3 is different.
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const keep = try idFilter(allocator, m, 0.9, false);
    defer allocator.free(keep);
    // s1 kept, s2 removed (identical to s1), s3 kept (different)
    try std.testing.expect(keep[0]);
    try std.testing.expect(!keep[1]);
    try std.testing.expect(keep[2]);
}

test "idFilter: all unique sequences kept" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "AAAA", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    const keep = try idFilter(allocator, m, 0.5, false);
    defer allocator.free(keep);
    try std.testing.expect(keep[0]);
    try std.testing.expect(keep[1]);
}

test "idFilter: conscover prefers sequence with more consensus residues" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    // s1: AC-T  (3 residues in consensus cols 0,2,3: A,-,T -> 2 consensus residues)
    // s2: ACGT  (4 residues in consensus cols 0,2,3: A,G,T -> 3 consensus residues)
    // s3: TTTT  (different enough to keep)
    // RF: x.xx  (col 1 is insert)
    // s1 and s2 are similar (pid high). With conscover, s1 should be removed
    // because s2 has more consensus coverage.
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "AC-T", "ACGT", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    m.reference = try m.allocator.dupe(u8, "x.xx");

    // s1 vs s2: A,C,gap,T vs A,C,G,T
    // len1=3 (A,C,T), len2=4 (A,C,G,T), matches=3 (cols 0,1,3), pid=3/3=1.0
    // With conscover: s1 has 2 consensus residues (A,T), s2 has 3 (A,G,T)
    // -> remove s1 (lower coverage)
    const keep = try idFilter(allocator, m, 0.9, true);
    defer allocator.free(keep);
    try std.testing.expect(!keep[0]); // s1 removed (lower conscover)
    try std.testing.expect(keep[1]); // s2 kept (higher conscover)
    try std.testing.expect(keep[2]); // s3 kept (different)
}

test "idFilter: conscover falls back to default without reference" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;
    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "ACGT", "TTTT" };
    var m = try Msa.fromText(allocator, abc, &names, &seqs);
    defer m.deinit();
    // No reference set -> conscover has no effect, later seq removed.
    const keep = try idFilter(allocator, m, 0.9, true);
    defer allocator.free(keep);
    try std.testing.expect(keep[0]);
    try std.testing.expect(!keep[1]);
    try std.testing.expect(keep[2]);
}

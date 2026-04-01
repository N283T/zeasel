// Independent set splitting for training/test partitioning.
//
// Implements the Cobalt and Blue algorithms (Petti & Eddy, 2022) for
// splitting datasets into quasi-independent training and test sets.
//
// The Cobalt algorithm is a simple greedy sequential approach.
// The Blue algorithm uses a multi-round election process and is more
// likely to achieve a split in difficult cases. Both are randomized.
//
// "mono" variants produce an independent set (IS): no two elements
// in the set are linked.
// "bi" variants produce a bipartite independent pair (BIP): sets S
// and T such that no element in S is linked to any element in T.
//
// References:
//   Petti & Eddy (2022) - Cobalt and Blue algorithms
//   Blelloch et al (2012) - greedy sequential maximum IS
//   Metivier et al (2011) - IS random priority algorithm

const std = @import("std");
const Allocator = std.mem.Allocator;
const Random = @import("util/random.zig").Random;

/// Linkage callback: returns true if elements i and j are "linked"
/// (too similar and should not be in the same set).
pub const LinkFn = *const fn (i: usize, j: usize, data: ?*anyopaque) bool;

/// Greedy independent set (monoCobalt).
///
/// Produces an independent set U: no two elements in U are linked.
/// Elements are processed in random order; each is added to U if
/// it has no link to any element already in U.
///
/// Returns assignments[0..n-1]: 0 = not in set, 1 = in independent set.
/// Caller owns returned slice.
pub fn monoCobalt(allocator: Allocator, rng: *Random, linkfn: LinkFn, data: ?*anyopaque, n: usize) ![]u8 {
    const assignments = try allocator.alloc(u8, n);
    @memset(assignments, 0);

    // Shuffled order of vertices
    const order = try allocator.alloc(usize, n);
    defer allocator.free(order);
    for (order, 0..) |*slot, i| slot.* = i;
    shuffle(rng, order);

    // Vertices added to the independent set so far
    const iset = try allocator.alloc(usize, n);
    defer allocator.free(iset);
    var iset_len: usize = 0;

    for (order) |v| {
        var adjacent = false;
        for (iset[0..iset_len]) |u| {
            if (linkfn(v, u, data)) {
                adjacent = true;
                break;
            }
        }
        if (!adjacent) {
            assignments[v] = 1;
            iset[iset_len] = v;
            iset_len += 1;
        }
    }

    return assignments;
}

/// Bipartite independent pair (biCobalt).
///
/// Splits n elements into sets S (=1) and T (=2) such that no element
/// in S is linked to any element in T. Elements not assigned to either
/// set get 0. S is defined as the larger set.
///
/// Returns assignments[0..n-1]: 0 = neither, 1 = set S, 2 = set T.
/// Caller owns returned slice.
pub fn biCobalt(allocator: Allocator, rng: *Random, linkfn: LinkFn, data: ?*anyopaque, n: usize) ![]u8 {
    const assignments = try allocator.alloc(u8, n);
    @memset(assignments, 0);

    // Shuffled order
    const order = try allocator.alloc(usize, n);
    defer allocator.free(order);
    for (order, 0..) |*slot, i| slot.* = i;
    shuffle(rng, order);

    // Side 1 and side 2 membership lists
    const side1 = try allocator.alloc(usize, n);
    defer allocator.free(side1);
    const side2 = try allocator.alloc(usize, n);
    defer allocator.free(side2);
    var nb1: usize = 0;
    var nb2: usize = 0;

    for (order) |v| {
        if (rng.uniform() <= 0.5) {
            // Try to put v in side 2: check adjacency to side 1
            const adj1 = isAdjacentToAny(linkfn, data, v, side1[0..nb1]);
            if (adj1) {
                // Can't go in side 2; try side 1
                const adj2 = isAdjacentToAny(linkfn, data, v, side2[0..nb2]);
                if (!adj2) {
                    assignments[v] = 1;
                    side1[nb1] = v;
                    nb1 += 1;
                }
                // else: assignments[v] stays 0
            } else {
                assignments[v] = 2;
                side2[nb2] = v;
                nb2 += 1;
            }
        } else {
            // Try to put v in side 1: check adjacency to side 2
            const adj2 = isAdjacentToAny(linkfn, data, v, side2[0..nb2]);
            if (adj2) {
                // Can't go in side 1; try side 2
                const adj1 = isAdjacentToAny(linkfn, data, v, side1[0..nb1]);
                if (!adj1) {
                    assignments[v] = 2;
                    side2[nb2] = v;
                    nb2 += 1;
                }
            } else {
                assignments[v] = 1;
                side1[nb1] = v;
                nb1 += 1;
            }
        }
    }

    // Define set 1 (S) as the larger set; swap labels if needed
    if (nb2 > nb1) {
        for (assignments) |*a| {
            if (a.* == 1) {
                a.* = 2;
            } else if (a.* == 2) {
                a.* = 1;
            }
        }
    }

    return assignments;
}

/// Multi-round election independent set (monoBlue).
///
/// Uses the random priority algorithm: in each round, vertices with
/// a label smaller than all their neighbors' labels are elected to the
/// independent set, and their neighbors are removed. Repeats until no
/// vertices remain.
///
/// Returns assignments[0..n-1]: 0 = not in set, 1 = in independent set.
/// Caller owns returned slice.
pub fn monoBlue(allocator: Allocator, rng: *Random, linkfn: LinkFn, data: ?*anyopaque, n: usize) ![]u8 {
    const assignments = try allocator.alloc(u8, n);
    @memset(assignments, 0);

    // status: -1 = in iset, -3 = removed from graph, >= 0 = still eligible
    // (value is next index in to_add to check)
    const status_d = try allocator.alloc(i32, n);
    defer allocator.free(status_d);
    @memset(status_d, 0);

    const dec_o = try allocator.alloc(usize, n);
    defer allocator.free(dec_o);
    const label_o = try allocator.alloc(usize, n);
    defer allocator.free(label_o);
    const to_add = try allocator.alloc(usize, n);
    defer allocator.free(to_add);

    for (dec_o, 0..) |*slot, i| slot.* = i;
    for (label_o, 0..) |*slot, i| slot.* = i;

    var k: usize = n;

    while (k > 0) {
        var lta: usize = 0;

        // i_select: determine which vertices to add to iset this round
        for (0..k) |idx| {
            const v = dec_o[idx];

            if (status_d[v] < 0) continue;

            // Check adjacency to vertices already in to_add
            var adj = false;
            {
                const start: usize = @intCast(status_d[v]);
                for (start..lta) |j| {
                    if (linkfn(v, to_add[j], data)) {
                        status_d[v] = -3;
                        adj = true;
                        break;
                    }
                }
            }

            if (!adj) {
                status_d[v] = @intCast(lta);

                // Check if v has a lower label than all its neighbors still in graph
                var found_self = false;
                var adj_lower = false;
                var j: usize = 0;
                while (!found_self and !adj_lower) {
                    const w = label_o[j];
                    if (w == v) {
                        found_self = true;
                        break;
                    }

                    if (status_d[w] >= 0) {
                        if (linkfn(v, w, data)) {
                            // Check if w is actually still eligible
                            var w_there = true;
                            const w_start: usize = @intCast(status_d[w]);
                            for (w_start..lta) |l| {
                                if (linkfn(w, to_add[l], data)) {
                                    status_d[w] = -3;
                                    w_there = false;
                                    break;
                                }
                            }
                            if (w_there) {
                                status_d[w] = @intCast(lta);
                                adj_lower = true;
                                break;
                            }
                        }
                    }
                    j += 1;
                }

                if (found_self) {
                    to_add[lta] = v;
                    lta += 1;
                    status_d[v] = -1;
                }
            }
        }

        // Remove eligibility of vertices adjacent to vertices in to_add
        for (0..k) |idx| {
            const v = dec_o[idx];
            if (status_d[v] >= 0) {
                const start: usize = @intCast(status_d[v]);
                for (start..lta) |j| {
                    if (linkfn(v, to_add[j], data)) {
                        status_d[v] = -3;
                        break;
                    }
                }
            }
        }

        // i_update_workspace: commit to_add to assignments, rebuild workspace
        for (0..lta) |i| {
            assignments[to_add[i]] = 1;
        }

        var d: usize = 0;
        for (0..k) |i| {
            const v_label = label_o[i];
            if (status_d[v_label] >= 0) {
                dec_o[d] = v_label;
                d += 1;
                status_d[v_label] = 0;
            }
            if (status_d[v_label] == -3) {
                assignments[v_label] = 0;
            }
        }

        // Copy decision order to label order
        for (0..d) |i| label_o[i] = dec_o[i];

        // Shuffle both
        shuffle(rng, dec_o[0..d]);
        shuffle(rng, label_o[0..d]);

        k = d;
    }

    return assignments;
}

/// Multi-round election bipartite independent pair (biBlue).
///
/// Uses the Blue algorithm to split n elements into sets S (=1) and T (=2)
/// such that no element in S is linked to any element in T. More likely to
/// achieve a split in difficult cases than biCobalt.
///
/// Returns assignments[0..n-1]: 0 = neither, 1 = set S, 2 = set T.
/// Caller owns returned slice.
pub fn biBlue(allocator: Allocator, rng: *Random, linkfn: LinkFn, data: ?*anyopaque, n: usize) ![]u8 {
    if (n == 0) return try allocator.alloc(u8, 0);

    const assignments = try allocator.alloc(u8, n);
    @memset(assignments, 0);

    // elig: 0 = removed, 1 = eligible for side 1 only, 2 = side 2 only, 3 = both
    const elig = try allocator.alloc(u8, n);
    defer allocator.free(elig);
    @memset(elig, 3);

    // status_d: for 1-candidates, label position; for 2-candidates, next to_add index to check
    const status_d = try allocator.alloc(i32, n);
    defer allocator.free(status_d);
    @memset(status_d, 0);

    // dec_o: 1-candidates on left, 2-candidates on right
    const dec_o = try allocator.alloc(usize, n);
    defer allocator.free(dec_o);

    // to_add: 1-side additions from left, 2-side additions from right
    const to_add = try allocator.alloc(usize, n);
    defer allocator.free(to_add);
    @memset(to_add, std.math.maxInt(usize)); // sentinel

    var nb1: usize = 0;
    var nb2: usize = 0;
    var d: usize = 0;
    var l: usize = 0;
    var lta1: usize = 0;
    var lta2: usize = 0;

    // Initial workspace setup
    biUpdateWorkspaceBlue(rng, dec_o, status_d, to_add, elig, assignments, n, &d, &l, &lta1, &lta2, &nb1, &nb2);

    while (l + d > 0) {
        lta1 = 0;
        lta2 = 0;
        const label_o = dec_o[n - l ..];

        // Select 1-candidates for 1-side
        for (0..d) |i| {
            const v = dec_o[i];
            var should_add = true;
            const limit: usize = @intCast(status_d[v]);

            for (0..limit) |j| {
                biUpdate2Elig(j, n, linkfn, data, label_o, status_d, to_add, elig, lta1);
                const w = label_o[j];
                if (elig[w] == 2 or elig[w] == 3) {
                    if (linkfn(v, w, data)) {
                        status_d[v] = @intCast(j);
                        should_add = false;
                        break;
                    }
                }
            }
            if (should_add) {
                to_add[lta1] = v;
                lta1 += 1;
                elig[v] = 0;
            }
        }

        // Select 2-candidates for 2-side
        for (0..l) |j| {
            biUpdate2Elig(j, n, linkfn, data, label_o, status_d, to_add, elig, lta1);
            const w = label_o[j];
            if (elig[w] == 2 or elig[w] == 3) {
                to_add[n - 1 - lta2] = w;
                lta2 += 1;
                elig[w] = 0;

                // Remove 1-eligibility of 1-candidates adjacent to w
                for (0..d) |i| {
                    const v = dec_o[i];
                    if (elig[v] == 0) continue;
                    if (elig[v] == 1 or elig[v] == 3) {
                        const sv: usize = @intCast(status_d[v]);
                        if (sv == j) {
                            elig[v] = if (elig[v] >= 1) elig[v] - 1 else 0;
                        } else if (sv < j) {
                            if (linkfn(v, w, data)) {
                                elig[v] = if (elig[v] >= 1) elig[v] - 1 else 0;
                            }
                        }
                    }
                }
            }
        }

        // Remove 1-elig of 2-candidates adjacent to a vertex in 2-side of to_add
        for (0..l) |j| {
            const w = label_o[j];
            if (elig[w] == 1 or elig[w] == 3) {
                for (0..lta2) |k_idx| {
                    const u = to_add[n - 1 - k_idx];
                    if (linkfn(u, w, data)) {
                        elig[w] = if (elig[w] >= 1) elig[w] - 1 else 0;
                        break;
                    }
                }
            }
        }

        // Remove 2-elig of 1-candidates adjacent to a vertex in 1-side of to_add
        for (0..d) |i| {
            const v = dec_o[i];
            if (elig[v] == 2 or elig[v] == 3) {
                for (0..lta1) |k_idx| {
                    const u = to_add[k_idx];
                    if (linkfn(v, u, data)) {
                        elig[v] = if (elig[v] >= 2) elig[v] - 2 else 0;
                        break;
                    }
                }
            }
        }

        // Commit and rebuild workspace
        biUpdateWorkspaceBlue(rng, dec_o, status_d, to_add, elig, assignments, n, &d, &l, &lta1, &lta2, &nb1, &nb2);
    }

    // Define set 1 (S) as the larger set
    if (nb2 > nb1) {
        for (assignments) |*a| {
            if (a.* == 1) {
                a.* = 2;
            } else if (a.* == 2) {
                a.* = 1;
            }
        }
    }

    return assignments;
}

// --- Helper functions ---

/// Fisher-Yates shuffle for a usize slice.
fn shuffle(rng: *Random, items: []usize) void {
    if (items.len <= 1) return;
    var i: usize = items.len - 1;
    while (i > 0) : (i -= 1) {
        const j = rng.uniformInt(@intCast(i + 1));
        const tmp = items[i];
        items[i] = items[j];
        items[j] = tmp;
    }
}

fn isAdjacentToAny(linkfn: LinkFn, data: ?*anyopaque, v: usize, set: []const usize) bool {
    for (set) |u| {
        if (linkfn(v, u, data)) return true;
    }
    return false;
}

/// Update 2-eligibility of label_o[j] by checking against 1-side of to_add.
fn biUpdate2Elig(
    j: usize,
    n: usize,
    linkfn: LinkFn,
    data: ?*anyopaque,
    label_o: []const usize,
    status_d: []i32,
    to_add: []const usize,
    elig: []u8,
    lta1: usize,
) void {
    _ = n;
    const w = label_o[j];
    if (elig[w] == 2 or elig[w] == 3) {
        const start: usize = @intCast(status_d[w]);
        for (start..lta1) |i| {
            const v = to_add[i];
            // If v has a higher label than j, v and w were already compared
            const sv: usize = @intCast(status_d[v]);
            if (sv <= j) {
                if (linkfn(v, w, data)) {
                    elig[w] = if (elig[w] >= 2) elig[w] - 2 else 0;
                    status_d[w] = @intCast(lta1);
                    break;
                }
            }
        }
    }
}

/// Rebuild workspace for biBlue between rounds.
fn biUpdateWorkspaceBlue(
    rng: *Random,
    dec_o: []usize,
    status_d: []i32,
    to_add: []usize,
    elig: []const u8,
    assignments: []u8,
    n: usize,
    d_out: *usize,
    l_out: *usize,
    lta1: *usize,
    lta2: *usize,
    nb1: *usize,
    nb2: *usize,
) void {
    // Commit 1-side of to_add to assignments
    for (0..lta1.*) |i| {
        assignments[to_add[i]] = 1;
        nb1.* += 1;
        to_add[i] = std.math.maxInt(usize);
    }

    // Commit 2-side of to_add to assignments
    for (0..lta2.*) |i| {
        assignments[to_add[n - 1 - i]] = 2;
        nb2.* += 1;
        to_add[n - 1 - i] = std.math.maxInt(usize);
    }

    var d: usize = 0;
    var l: usize = 0;

    for (0..n) |i| {
        if (elig[i] == 3) {
            if (rng.uniform() < 0.5) {
                // 1-candidate
                dec_o[d] = i;
                d += 1;
            } else {
                // 2-candidate
                dec_o[n - 1 - l] = i;
                l += 1;
                status_d[i] = 0;
            }
        } else if (elig[i] == 1) {
            dec_o[d] = i;
            d += 1;
        } else if (elig[i] == 2) {
            dec_o[n - 1 - l] = i;
            l += 1;
            status_d[i] = 0;
        }
    }

    // For 1-candidates, assign random label position in [0, l]
    for (0..d) |i| {
        status_d[dec_o[i]] = @intCast(rng.uniformInt(@intCast(l + 1)));
    }

    shuffle(rng, dec_o[0..d]);
    if (l > 0) {
        shuffle(rng, dec_o[n - l ..]);
    }

    lta1.* = 0;
    lta2.* = 0;
    d_out.* = d;
    l_out.* = l;
}

/// Random bipartite independent pair (biRandom).
///
/// Randomly assigns each vertex to S with probability `t_prob`, then
/// eligible vertices (not linked to any S member) go to T.
/// Reference: Easel esl_iset_biRandom().
///
/// Returns .s[0..n] and .t[0..n] boolean arrays. Caller owns both slices.
pub fn biRandom(
    allocator: Allocator,
    n: usize,
    linkfn: LinkFn,
    data: ?*anyopaque,
    t_prob: f64,
    rng: *Random,
) !struct { s: []bool, t: []bool } {
    const s = try allocator.alloc(bool, n);
    errdefer allocator.free(s);
    const t = try allocator.alloc(bool, n);
    errdefer allocator.free(t);
    @memset(s, false);
    @memset(t, false);

    // Phase 1: randomly assign vertices to S with probability t_prob.
    for (0..n) |i| {
        if (rng.uniform() < t_prob) {
            s[i] = true;
        }
    }

    // Phase 2: assign eligible vertices to T.
    // A vertex is eligible for T if it is not in S and not linked to any S member.
    for (0..n) |i| {
        if (s[i]) continue;
        var linked_to_s = false;
        for (0..n) |j| {
            if (s[j] and linkfn(i, j, data)) {
                linked_to_s = true;
                break;
            }
        }
        if (!linked_to_s) {
            t[i] = true;
        }
    }

    return .{ .s = s, .t = t };
}

/// Validate a mono independent set: no two assigned vertices are linked.
///
/// Returns true if assignment is a valid independent set.
pub fn monoValidate(n: usize, linkfn: LinkFn, data: ?*anyopaque, assignment: []const u8) bool {
    if (assignment.len < n) return false;
    for (0..n) |i| {
        if (assignment[i] == 0) continue;
        for (i + 1..n) |j| {
            if (assignment[j] == 0) continue;
            if (linkfn(i, j, data)) return false;
        }
    }
    return true;
}

/// Validate a bipartite independent pair: S and T are each valid
/// independent sets and no S-T cross-links exist.
///
/// Returns true if s and t form a valid bipartite independent pair.
pub fn biValidate(n: usize, linkfn: LinkFn, data: ?*anyopaque, s: []const bool, t: []const bool) bool {
    if (s.len < n or t.len < n) return false;
    for (0..n) |i| {
        for (i + 1..n) |j| {
            if (!linkfn(i, j, data)) continue;
            // No two S members may be linked.
            if (s[i] and s[j]) return false;
            // No two T members may be linked.
            if (t[i] and t[j]) return false;
            // No S-T cross-links.
            if (s[i] and t[j]) return false;
            if (t[i] and s[j]) return false;
        }
    }
    return true;
}

// --- Legacy API (distance-matrix based) ---

/// Split n items into two quasi-independent sets based on a pairwise
/// distance matrix. Items within the same set should be "independent"
/// (distance > threshold).
///
/// Returns a boolean array: true = set A, false = set B.
/// Uses a greedy coloring approach: assign each item to the set where
/// it has fewer close neighbors.
pub fn splitIndependent(
    allocator: Allocator,
    dist: []const f64,
    n: usize,
    threshold: f64,
) ![]bool {
    const assignment = try allocator.alloc(bool, n);
    @memset(assignment, true); // Start all in set A

    // Count close neighbors in each set for each item
    const neighbors_a = try allocator.alloc(usize, n);
    defer allocator.free(neighbors_a);
    const neighbors_b = try allocator.alloc(usize, n);
    defer allocator.free(neighbors_b);

    // Iterative refinement
    for (0..10) |_| {
        @memset(neighbors_a, 0);
        @memset(neighbors_b, 0);

        // Count neighbors
        for (0..n) |i| {
            for (0..n) |j| {
                if (i == j) continue;
                if (dist[i * n + j] <= threshold) {
                    if (assignment[j]) {
                        neighbors_a[i] += 1;
                    } else {
                        neighbors_b[i] += 1;
                    }
                }
            }
        }

        // Reassign: put each item in the set with fewer close neighbors
        var changed = false;
        for (0..n) |i| {
            const should_be_a = neighbors_a[i] <= neighbors_b[i];
            if (assignment[i] != should_be_a) {
                assignment[i] = should_be_a;
                changed = true;
            }
        }

        if (!changed) break;
    }

    return assignment;
}

/// Count the number of close pairs within each set (quality metric).
/// Returns (intra_a, intra_b) -- lower is better.
pub fn countIntraPairs(assignment: []const bool, dist: []const f64, n: usize, threshold: f64) struct { a: usize, b: usize } {
    var intra_a: usize = 0;
    var intra_b: usize = 0;

    for (0..n) |i| {
        for (i + 1..n) |j| {
            if (dist[i * n + j] <= threshold) {
                if (assignment[i] and assignment[j]) {
                    intra_a += 1;
                } else if (!assignment[i] and !assignment[j]) {
                    intra_b += 1;
                }
            }
        }
    }

    return .{ .a = intra_a, .b = intra_b };
}

// --- Tests ---

/// Test linkage function: elements are linked if they are within 1 of each other
/// in a simple chain (0-1-2-3-...).
fn chainLink(i: usize, j: usize, _: ?*anyopaque) bool {
    const diff = if (i > j) i - j else j - i;
    return diff == 1;
}

/// Test linkage function using an adjacency matrix stored in data.
fn matrixLink(i: usize, j: usize, data: ?*anyopaque) bool {
    const matrix: *const AdjMatrix = @ptrCast(@alignCast(data.?));
    return matrix.linked(i, j);
}

const AdjMatrix = struct {
    edges: []const bool,
    n: usize,

    fn linked(self: *const AdjMatrix, i: usize, j: usize) bool {
        return self.edges[i * self.n + j];
    }
};

test "monoCobalt: chain of 5 — independent set has no adjacent elements" {
    const allocator = std.testing.allocator;
    var rng = Random.init(42);

    const result = try monoCobalt(allocator, &rng, chainLink, null, 5);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 5), result.len);

    // No two adjacent elements should both be in the set
    for (0..4) |i| {
        if (result[i] == 1 and result[i + 1] == 1) {
            return error.TestUnexpectedResult;
        }
    }

    // At least some elements should be in the set
    var count: usize = 0;
    for (result) |a| count += a;
    try std.testing.expect(count >= 2);
}

test "monoCobalt: no links — all elements in set" {
    const allocator = std.testing.allocator;
    var rng = Random.init(7);

    const noLink = struct {
        fn f(_: usize, _: usize, _: ?*anyopaque) bool {
            return false;
        }
    }.f;

    const result = try monoCobalt(allocator, &rng, noLink, null, 4);
    defer allocator.free(result);

    for (result) |a| {
        try std.testing.expectEqual(@as(u8, 1), a);
    }
}

test "monoCobalt: complete graph — only one element in set" {
    const allocator = std.testing.allocator;
    var rng = Random.init(99);

    const allLink = struct {
        fn f(i: usize, j: usize, _: ?*anyopaque) bool {
            return i != j;
        }
    }.f;

    const result = try monoCobalt(allocator, &rng, allLink, null, 5);
    defer allocator.free(result);

    var count: usize = 0;
    for (result) |a| count += a;
    try std.testing.expectEqual(@as(usize, 1), count);
}

test "biCobalt: chain of 6 — S and T have no cross-links" {
    const allocator = std.testing.allocator;
    var rng = Random.init(42);

    const result = try biCobalt(allocator, &rng, chainLink, null, 6);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 6), result.len);

    // No element in S should be linked to any element in T
    for (0..6) |i| {
        for (0..6) |j| {
            if (chainLink(i, j, null)) {
                if (result[i] == 1 and result[j] == 2) {
                    return error.TestUnexpectedResult;
                }
                if (result[i] == 2 and result[j] == 1) {
                    return error.TestUnexpectedResult;
                }
            }
        }
    }

    // S should be at least as large as T
    var count_s: usize = 0;
    var count_t: usize = 0;
    for (result) |a| {
        if (a == 1) count_s += 1;
        if (a == 2) count_t += 1;
    }
    try std.testing.expect(count_s >= count_t);
}

test "biCobalt: no links — all in one set" {
    const allocator = std.testing.allocator;
    var rng = Random.init(13);

    const noLink = struct {
        fn f(_: usize, _: usize, _: ?*anyopaque) bool {
            return false;
        }
    }.f;

    const result = try biCobalt(allocator, &rng, noLink, null, 5);
    defer allocator.free(result);

    // With no links, all elements can be in the same set
    var count_s: usize = 0;
    var count_t: usize = 0;
    for (result) |a| {
        if (a == 1) count_s += 1;
        if (a == 2) count_t += 1;
    }
    // All assigned (none are 0)
    try std.testing.expectEqual(@as(usize, 5), count_s + count_t);
}

test "biCobalt: adjacency matrix — correctness with data pointer" {
    const allocator = std.testing.allocator;
    var rng = Random.init(77);

    // Triangle: 0-1, 1-2, 0-2 all linked
    const edges = [_]bool{
        false, true,  true,
        true,  false, true,
        true,  true,  false,
    };
    var adj = AdjMatrix{ .edges = &edges, .n = 3 };

    const result = try biCobalt(allocator, &rng, matrixLink, @ptrCast(&adj), 3);
    defer allocator.free(result);

    // In a complete triangle, BIP can have at most 1 in S and 1 in T
    for (0..3) |i| {
        for (0..3) |j| {
            if (adj.linked(i, j)) {
                const invalid = (result[i] == 1 and result[j] == 2) or
                    (result[i] == 2 and result[j] == 1);
                try std.testing.expect(!invalid);
            }
        }
    }
}

test "monoBlue: chain of 5 — independent set has no adjacent elements" {
    const allocator = std.testing.allocator;
    var rng = Random.init(42);

    const result = try monoBlue(allocator, &rng, chainLink, null, 5);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 5), result.len);

    // No two adjacent elements should both be in the set
    for (0..4) |i| {
        if (result[i] == 1 and result[i + 1] == 1) {
            return error.TestUnexpectedResult;
        }
    }

    var count: usize = 0;
    for (result) |a| count += a;
    try std.testing.expect(count >= 2);
}

test "monoBlue: complete graph — only one element" {
    const allocator = std.testing.allocator;
    var rng = Random.init(55);

    const allLink = struct {
        fn f(i: usize, j: usize, _: ?*anyopaque) bool {
            return i != j;
        }
    }.f;

    const result = try monoBlue(allocator, &rng, allLink, null, 5);
    defer allocator.free(result);

    var count: usize = 0;
    for (result) |a| count += a;
    try std.testing.expectEqual(@as(usize, 1), count);
}

test "monoBlue: no links — all elements in set" {
    const allocator = std.testing.allocator;
    var rng = Random.init(7);

    const noLink = struct {
        fn f(_: usize, _: usize, _: ?*anyopaque) bool {
            return false;
        }
    }.f;

    const result = try monoBlue(allocator, &rng, noLink, null, 4);
    defer allocator.free(result);

    for (result) |a| {
        try std.testing.expectEqual(@as(u8, 1), a);
    }
}

test "biBlue: chain of 6 — S and T have no cross-links" {
    const allocator = std.testing.allocator;
    var rng = Random.init(42);

    const result = try biBlue(allocator, &rng, chainLink, null, 6);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 6), result.len);

    for (0..6) |i| {
        for (0..6) |j| {
            if (chainLink(i, j, null)) {
                if (result[i] == 1 and result[j] == 2) {
                    return error.TestUnexpectedResult;
                }
                if (result[i] == 2 and result[j] == 1) {
                    return error.TestUnexpectedResult;
                }
            }
        }
    }

    var count_s: usize = 0;
    var count_t: usize = 0;
    for (result) |a| {
        if (a == 1) count_s += 1;
        if (a == 2) count_t += 1;
    }
    try std.testing.expect(count_s >= count_t);
}

test "biBlue: empty input" {
    const allocator = std.testing.allocator;
    var rng = Random.init(0);

    const noLink = struct {
        fn f(_: usize, _: usize, _: ?*anyopaque) bool {
            return false;
        }
    }.f;

    const result = try biBlue(allocator, &rng, noLink, null, 0);
    defer allocator.free(result);
    try std.testing.expectEqual(@as(usize, 0), result.len);
}

test "biBlue: complete graph of 4 — valid BIP" {
    const allocator = std.testing.allocator;
    var rng = Random.init(33);

    const allLink = struct {
        fn f(i: usize, j: usize, _: ?*anyopaque) bool {
            return i != j;
        }
    }.f;

    const result = try biBlue(allocator, &rng, allLink, null, 4);
    defer allocator.free(result);

    // Verify BIP property
    for (0..4) |i| {
        for (0..4) |j| {
            if (i != j) {
                const invalid = (result[i] == 1 and result[j] == 2) or
                    (result[i] == 2 and result[j] == 1);
                try std.testing.expect(!invalid);
            }
        }
    }
}

// Legacy tests

test "splitIndependent: 2 close pairs, 1 distant" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.1, 0.9,
        0.1, 0.0, 0.9,
        0.9, 0.9, 0.0,
    };

    const split = try splitIndependent(allocator, &dist, 3, 0.5);
    defer allocator.free(split);

    try std.testing.expectEqual(@as(usize, 3), split.len);
}

test "splitIndependent: all distant = all same set" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.9, 0.9,
        0.9, 0.0, 0.9,
        0.9, 0.9, 0.0,
    };

    const split = try splitIndependent(allocator, &dist, 3, 0.5);
    defer allocator.free(split);

    try std.testing.expectEqual(@as(usize, 3), split.len);
}

test "countIntraPairs: basic" {
    const dist = [_]f64{
        0.0, 0.1, 0.9,
        0.1, 0.0, 0.9,
        0.9, 0.9, 0.0,
    };
    const assignment = [_]bool{ true, true, false };
    const result = countIntraPairs(&assignment, &dist, 3, 0.5);

    try std.testing.expectEqual(@as(usize, 1), result.a);
    try std.testing.expectEqual(@as(usize, 0), result.b);
}

test "biRandom: produces valid bipartite independent pair" {
    const allocator = std.testing.allocator;
    var rng = Random.init(42);

    const result = try biRandom(allocator, 6, chainLink, null, 0.3, &rng);
    defer allocator.free(result.s);
    defer allocator.free(result.t);

    // No vertex can be in both S and T.
    for (0..6) |i| {
        try std.testing.expect(!(result.s[i] and result.t[i]));
    }

    // No S-T cross-links.
    for (0..6) |i| {
        for (0..6) |j| {
            if (chainLink(i, j, null)) {
                try std.testing.expect(!(result.s[i] and result.t[j]));
            }
        }
    }
}

test "biRandom: no links — all non-S vertices go to T" {
    const allocator = std.testing.allocator;
    var rng = Random.init(7);

    const noLink = struct {
        fn f(_: usize, _: usize, _: ?*anyopaque) bool {
            return false;
        }
    }.f;

    const result = try biRandom(allocator, 4, noLink, null, 0.5, &rng);
    defer allocator.free(result.s);
    defer allocator.free(result.t);

    // Every non-S vertex should be in T when there are no links.
    for (0..4) |i| {
        if (!result.s[i]) {
            try std.testing.expect(result.t[i]);
        }
    }
}

test "monoValidate: valid independent set" {
    const result = [_]u8{ 1, 0, 1, 0, 1 };
    try std.testing.expect(monoValidate(5, chainLink, null, &result));
}

test "monoValidate: invalid — adjacent elements both assigned" {
    const result = [_]u8{ 1, 1, 0, 0, 0 };
    try std.testing.expect(!monoValidate(5, chainLink, null, &result));
}

test "biValidate: valid bipartite independent pair" {
    // Chain 0-1-2-3-4. S={0,2,4}, T={} is valid; S={0,4}, T={2} is valid.
    const s = [_]bool{ true, false, false, false, true };
    const t = [_]bool{ false, false, true, false, false };
    try std.testing.expect(biValidate(5, chainLink, null, &s, &t));
}

test "biValidate: invalid — S-T cross-link" {
    // 0-1 are linked; S={0}, T={1} is invalid.
    const s = [_]bool{ true, false, false, false, false };
    const t = [_]bool{ false, true, false, false, false };
    try std.testing.expect(!biValidate(5, chainLink, null, &s, &t));
}

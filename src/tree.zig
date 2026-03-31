// Phylogenetic tree construction and serialisation.
// Currently implements UPGMA clustering from a flat pairwise distance matrix
// and Newick format output.

const std = @import("std");
const Allocator = std.mem.Allocator;

/// A rooted binary tree stored in arrays indexed by node ID.
/// Node IDs 0..(n_leaves-1) are leaves; n_leaves..(n_nodes-1) are internal.
pub const Tree = struct {
    /// parent[i] = index of parent node, or -1 for the root.
    parent: []i32,
    /// left[i] = left child index for internal nodes, -1 for leaves.
    left: []i32,
    /// right[i] = right child index for internal nodes, -1 for leaves.
    right: []i32,
    /// Branch length from node i to its parent.
    branch_length: []f64,
    /// Leaf name for leaf nodes; null for internal nodes.
    names: []?[]const u8,
    /// Total number of nodes (leaves + internal).
    n_nodes: usize,
    /// Number of leaf nodes.
    n_leaves: usize,
    allocator: Allocator,

    pub fn deinit(self: *Tree) void {
        self.allocator.free(self.parent);
        self.allocator.free(self.left);
        self.allocator.free(self.right);
        self.allocator.free(self.branch_length);
        for (self.names) |n| if (n) |name| self.allocator.free(name);
        self.allocator.free(self.names);
    }
};

/// Build a UPGMA tree from a distance matrix.
///
/// `dist` is a flat n*n row-major distance matrix.
/// `leaf_names[0..n-1]` are the names of the n leaves.
///
/// The UPGMA algorithm iteratively merges the two closest clusters.
/// Distances between the new cluster and all remaining clusters are
/// computed as the weighted average of the component cluster distances.
pub fn upgma(allocator: Allocator, dist: []const f64, n: usize, leaf_names: []const []const u8) !Tree {
    const total_nodes = 2 * n - 1; // n leaves + n-1 internal nodes

    var parent = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(parent);
    @memset(parent, -1);

    var left_arr = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(left_arr);
    @memset(left_arr, -1);

    var right_arr = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(right_arr);
    @memset(right_arr, -1);

    var bl = try allocator.alloc(f64, total_nodes);
    errdefer allocator.free(bl);
    @memset(bl, 0.0);

    var names_arr = try allocator.alloc(?[]const u8, total_nodes);
    errdefer {
        for (names_arr) |nm| if (nm) |name| allocator.free(name);
        allocator.free(names_arr);
    }
    @memset(names_arr, null);

    // Duplicate leaf names into our allocation.
    for (0..n) |i| {
        names_arr[i] = try allocator.dupe(u8, leaf_names[i]);
    }

    // Working distance matrix, indexed over total_nodes (only 0..next_node-1 used).
    var d = try allocator.alloc(f64, total_nodes * total_nodes);
    defer allocator.free(d);
    @memset(d, std.math.inf(f64));
    for (0..n) |i| {
        for (0..n) |j| {
            d[i * total_nodes + j] = dist[i * n + j];
        }
    }

    // Cluster sizes (number of leaves in each cluster).
    var size = try allocator.alloc(f64, total_nodes);
    defer allocator.free(size);
    for (0..total_nodes) |i| size[i] = 1.0;

    // Track which cluster indices are still active.
    var active = try allocator.alloc(bool, total_nodes);
    defer allocator.free(active);
    @memset(active, false);
    for (0..n) |i| active[i] = true;

    // Height of each node (distance from node to any of its leaves).
    // Initialised to zero; updated when a node is merged into a parent.
    var height = try allocator.alloc(f64, total_nodes);
    defer allocator.free(height);
    @memset(height, 0.0);

    var next_node: usize = n;

    for (0..n - 1) |_| {
        // Find the pair of active clusters with minimum distance.
        var min_d: f64 = std.math.inf(f64);
        var mi: usize = 0;
        var mj: usize = 0;
        for (0..next_node) |i| {
            if (!active[i]) continue;
            for (i + 1..next_node) |j| {
                if (!active[j]) continue;
                const dij = d[i * total_nodes + j];
                if (dij < min_d) {
                    min_d = dij;
                    mi = i;
                    mj = j;
                }
            }
        }

        // Create a new internal node joining mi and mj.
        const new_node = next_node;
        left_arr[new_node] = @intCast(mi);
        right_arr[new_node] = @intCast(mj);
        parent[mi] = @intCast(new_node);
        parent[mj] = @intCast(new_node);

        // UPGMA height: half the distance between the two clusters.
        const node_height = min_d / 2.0;
        height[new_node] = node_height;
        bl[mi] = node_height - height[mi];
        bl[mj] = node_height - height[mj];

        // Update distances: weighted average by cluster size.
        for (0..next_node) |k| {
            if (!active[k] or k == mi or k == mj) continue;
            const di = d[mi * total_nodes + k];
            const dj = d[mj * total_nodes + k];
            const new_d = (size[mi] * di + size[mj] * dj) / (size[mi] + size[mj]);
            d[new_node * total_nodes + k] = new_d;
            d[k * total_nodes + new_node] = new_d;
        }

        size[new_node] = size[mi] + size[mj];
        active[mi] = false;
        active[mj] = false;
        active[new_node] = true;
        next_node += 1;
    }

    return Tree{
        .parent = parent,
        .left = left_arr,
        .right = right_arr,
        .branch_length = bl,
        .names = names_arr,
        .n_nodes = total_nodes,
        .n_leaves = n,
        .allocator = allocator,
    };
}

/// Write the tree in Newick format to `dest`.
/// The root is the last node (index n_nodes - 1).
pub fn writeNewick(tree: Tree, dest: std.io.AnyWriter) !void {
    try writeNewickNode(tree, dest, @intCast(tree.n_nodes - 1));
    try dest.writeAll(";\n");
}

fn writeNewickNode(tree: Tree, dest: std.io.AnyWriter, node: i32) !void {
    if (node < 0) return;
    const idx: usize = @intCast(node);
    if (tree.left[idx] >= 0) {
        // Internal node: recurse into children.
        try dest.writeByte('(');
        try writeNewickNode(tree, dest, tree.left[idx]);
        try dest.writeByte(',');
        try writeNewickNode(tree, dest, tree.right[idx]);
        try dest.writeByte(')');
    } else if (tree.names[idx]) |name| {
        try dest.writeAll(name);
    }
    // Print branch length for every node except the root (which has no parent).
    if (tree.parent[idx] >= 0) {
        try dest.print(":{d:.6}", .{tree.branch_length[idx]});
    }
}

/// Parse a Newick format string into a Tree.
/// Supports labeled leaves, branch lengths, and nested clades.
/// Example: "((A:0.1,B:0.2):0.3,C:0.4);"
pub fn readNewick(allocator: Allocator, input: []const u8) !Tree {
    // First pass: count leaves and internal nodes
    var n_leaves: usize = 0;
    var n_internal: usize = 0;
    var i: usize = 0;
    while (i < input.len) : (i += 1) {
        switch (input[i]) {
            '(' => n_internal += 1,
            ',' => {},
            ')' => {},
            ';', '\n', '\r', ' ', '\t' => {},
            ':' => {
                // Skip branch length
                i += 1;
                while (i < input.len and input[i] != ',' and input[i] != ')' and input[i] != ';' and input[i] != '(' and input[i] != ' ') : (i += 1) {}
                if (i < input.len) i -= 1; // will be incremented by loop
            },
            else => {
                // Start of a leaf name
                n_leaves += 1;
                while (i + 1 < input.len and input[i + 1] != ':' and input[i + 1] != ',' and input[i + 1] != ')' and input[i + 1] != ';' and input[i + 1] != '(' and input[i + 1] != ' ') : (i += 1) {}
            },
        }
    }

    if (n_leaves == 0) return error.InvalidInput;
    const total_nodes = n_leaves + n_internal;

    var parent = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(parent);
    @memset(parent, -1);

    var left_arr = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(left_arr);
    @memset(left_arr, -1);

    var right_arr = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(right_arr);
    @memset(right_arr, -1);

    var bl = try allocator.alloc(f64, total_nodes);
    errdefer allocator.free(bl);
    @memset(bl, 0.0);

    var names_arr = try allocator.alloc(?[]const u8, total_nodes);
    errdefer {
        for (names_arr) |nm| if (nm) |name| allocator.free(name);
        allocator.free(names_arr);
    }
    @memset(names_arr, null);

    // Stack for parsing
    var stack: std.ArrayList(usize) = .empty;
    defer stack.deinit(allocator);

    var next_leaf: usize = 0;
    var next_internal: usize = n_leaves;
    var current_node: ?usize = null;

    i = 0;
    while (i < input.len) : (i += 1) {
        switch (input[i]) {
            '(' => {
                // Push current context, start new internal node
                const new_node = next_internal;
                next_internal += 1;
                try stack.append(allocator, new_node);
                current_node = null;
            },
            ',' => {
                // Attach current_node as a child of the top-of-stack internal node
                if (current_node) |cn| {
                    if (stack.items.len > 0) {
                        const parent_node = stack.items[stack.items.len - 1];
                        if (left_arr[parent_node] < 0) {
                            left_arr[parent_node] = @intCast(cn);
                            parent[cn] = @intCast(parent_node);
                        } else if (right_arr[parent_node] < 0) {
                            right_arr[parent_node] = @intCast(cn);
                            parent[cn] = @intCast(parent_node);
                        } else {
                            // Multifurcation: more than 2 children
                            return error.MultifurcationNotSupported;
                        }
                    }
                }
                current_node = null;
            },
            ')' => {
                // Close the internal node
                if (stack.items.len == 0) return error.InvalidInput;
                const internal_node = stack.pop().?;

                // Attach current_node as child
                if (current_node) |cn| {
                    if (left_arr[internal_node] < 0) {
                        left_arr[internal_node] = @intCast(cn);
                    } else if (right_arr[internal_node] < 0) {
                        right_arr[internal_node] = @intCast(cn);
                    } else {
                        return error.MultifurcationNotSupported;
                    }
                    parent[cn] = @intCast(internal_node);
                }

                current_node = internal_node;
            },
            ':' => {
                // Parse branch length for current_node
                i += 1;
                const start = i;
                while (i < input.len and input[i] != ',' and input[i] != ')' and input[i] != ';' and input[i] != '(' and input[i] != ' ' and input[i] != '\n') : (i += 1) {}
                if (current_node) |cn| {
                    bl[cn] = std.fmt.parseFloat(f64, input[start..i]) catch 0.0;
                }
                if (i < input.len) i -= 1;
            },
            ';', '\n', '\r', ' ', '\t' => {},
            else => {
                // Leaf name
                const start = i;
                while (i + 1 < input.len and input[i + 1] != ':' and input[i + 1] != ',' and input[i + 1] != ')' and input[i + 1] != ';' and input[i + 1] != '(') : (i += 1) {}
                const leaf_node = next_leaf;
                next_leaf += 1;
                names_arr[leaf_node] = try allocator.dupe(u8, input[start .. i + 1]);

                current_node = leaf_node;
            },
        }
    }

    return Tree{
        .parent = parent,
        .left = left_arr,
        .right = right_arr,
        .branch_length = bl,
        .names = names_arr,
        .n_nodes = total_nodes,
        .n_leaves = n_leaves,
        .allocator = allocator,
    };
}

/// Build a tree using single-linkage clustering from a distance matrix.
/// At each step, merges the two closest clusters. Distance to merged cluster
/// is the minimum of distances to the two components.
pub fn singleLinkage(allocator: Allocator, dist: []const f64, n: usize, leaf_names: []const []const u8) !Tree {
    return genericLinkage(allocator, dist, n, leaf_names, .single);
}

/// Build a tree using complete-linkage clustering from a distance matrix.
/// Distance to merged cluster is the maximum of distances to the two components.
pub fn completeLinkage(allocator: Allocator, dist: []const f64, n: usize, leaf_names: []const []const u8) !Tree {
    return genericLinkage(allocator, dist, n, leaf_names, .complete);
}

/// Build a tree using WPGMA (weighted pair group method with arithmetic mean).
/// Like UPGMA but uses unweighted average (each subtree weighted equally regardless of size).
pub fn wpgma(allocator: Allocator, dist: []const f64, n: usize, leaf_names: []const []const u8) !Tree {
    return genericLinkage(allocator, dist, n, leaf_names, .wpgma);
}

const LinkageType = enum { single, complete, wpgma };

fn genericLinkage(allocator: Allocator, dist: []const f64, n: usize, leaf_names: []const []const u8, linkage: LinkageType) !Tree {
    const total_nodes = 2 * n - 1;

    var parent_arr = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(parent_arr);
    @memset(parent_arr, -1);

    var left_arr = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(left_arr);
    @memset(left_arr, -1);

    var right_arr = try allocator.alloc(i32, total_nodes);
    errdefer allocator.free(right_arr);
    @memset(right_arr, -1);

    var bl = try allocator.alloc(f64, total_nodes);
    errdefer allocator.free(bl);
    @memset(bl, 0.0);

    var names_arr = try allocator.alloc(?[]const u8, total_nodes);
    errdefer {
        for (names_arr) |nm| if (nm) |name_v| allocator.free(name_v);
        allocator.free(names_arr);
    }
    @memset(names_arr, null);
    for (0..n) |idx| {
        names_arr[idx] = try allocator.dupe(u8, leaf_names[idx]);
    }

    var d = try allocator.alloc(f64, total_nodes * total_nodes);
    defer allocator.free(d);
    @memset(d, std.math.inf(f64));
    for (0..n) |ii| {
        for (0..n) |jj| {
            d[ii * total_nodes + jj] = dist[ii * n + jj];
        }
    }

    var active = try allocator.alloc(bool, total_nodes);
    defer allocator.free(active);
    @memset(active, false);
    for (0..n) |idx| active[idx] = true;

    var height = try allocator.alloc(f64, total_nodes);
    defer allocator.free(height);
    @memset(height, 0.0);

    var next_node: usize = n;

    for (0..n - 1) |_| {
        var min_d: f64 = std.math.inf(f64);
        var mi: usize = 0;
        var mj: usize = 0;
        for (0..next_node) |ii| {
            if (!active[ii]) continue;
            for (ii + 1..next_node) |jj| {
                if (!active[jj]) continue;
                const dij = d[ii * total_nodes + jj];
                if (dij < min_d) {
                    min_d = dij;
                    mi = ii;
                    mj = jj;
                }
            }
        }

        const new_node = next_node;
        left_arr[new_node] = @intCast(mi);
        right_arr[new_node] = @intCast(mj);
        parent_arr[mi] = @intCast(new_node);
        parent_arr[mj] = @intCast(new_node);

        const node_height = min_d / 2.0;
        height[new_node] = node_height;
        bl[mi] = node_height - height[mi];
        bl[mj] = node_height - height[mj];

        for (0..next_node) |k| {
            if (!active[k] or k == mi or k == mj) continue;
            const di = d[mi * total_nodes + k];
            const dj = d[mj * total_nodes + k];
            const new_d = switch (linkage) {
                .single => @min(di, dj),
                .complete => @max(di, dj),
                .wpgma => (di + dj) / 2.0,
            };
            d[new_node * total_nodes + k] = new_d;
            d[k * total_nodes + new_node] = new_d;
        }

        active[mi] = false;
        active[mj] = false;
        active[new_node] = true;
        next_node += 1;
    }

    return Tree{
        .parent = parent_arr,
        .left = left_arr,
        .right = right_arr,
        .branch_length = bl,
        .names = names_arr,
        .n_nodes = total_nodes,
        .n_leaves = n,
        .allocator = allocator,
    };
}

// --- Tests ---

test "upgma: 3 sequences with known distances" {
    const allocator = std.testing.allocator;

    // Symmetric 3x3 distance matrix:
    //      A    B    C
    // A  [ 0   0.4  0.8 ]
    // B  [ 0.4  0   0.8 ]
    // C  [ 0.8  0.8  0  ]
    //
    // Minimum distance is A-B = 0.4, merged first at height 0.2.
    // Then (AB) vs C: distance = (0.8 + 0.8) / 2 = 0.8, merged at height 0.4.
    const dist = [_]f64{
        0.0, 0.4, 0.8,
        0.4, 0.0, 0.8,
        0.8, 0.8, 0.0,
    };
    const names = [_][]const u8{ "A", "B", "C" };

    var tree = try upgma(allocator, &dist, 3, &names);
    defer tree.deinit();

    // Total nodes = 2*3 - 1 = 5
    try std.testing.expectEqual(@as(usize, 5), tree.n_nodes);
    try std.testing.expectEqual(@as(usize, 3), tree.n_leaves);

    // The root is the last internal node (index 4).
    const root: usize = tree.n_nodes - 1;
    try std.testing.expectEqual(@as(i32, -1), tree.parent[@intCast(root)]);

    // The first internal node (index 3) merges A(0) and B(1).
    // Branch lengths from A and B to node 3 should both be 0.2.
    const ab_node: usize = 3;
    try std.testing.expectEqual(@as(i32, 0), tree.left[ab_node]);
    try std.testing.expectEqual(@as(i32, 1), tree.right[ab_node]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), tree.branch_length[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), tree.branch_length[1], 1e-10);

    // Branch length from node 3 (AB) to root should be 0.2 (0.4 - 0.2).
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), tree.branch_length[ab_node], 1e-10);

    // Branch length from C(2) to root should be 0.4.
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), tree.branch_length[2], 1e-10);
}

test "upgma: Newick output" {
    const allocator = std.testing.allocator;

    // Simple 2-leaf star tree: dist(A,B) = 0.6 -> each branch = 0.3
    const dist = [_]f64{
        0.0, 0.6,
        0.6, 0.0,
    };
    const names = [_][]const u8{ "A", "B" };

    var tree = try upgma(allocator, &dist, 2, &names);
    defer tree.deinit();

    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try writeNewick(tree, buf.writer(allocator).any());

    // Expected: "(A:0.300000,B:0.300000);\n"
    const newick = buf.items;
    try std.testing.expect(std.mem.startsWith(u8, newick, "("));
    try std.testing.expect(std.mem.endsWith(u8, newick, ";\n"));
    try std.testing.expect(std.mem.indexOf(u8, newick, "A:") != null);
    try std.testing.expect(std.mem.indexOf(u8, newick, "B:") != null);
}

test "upgma: single pair branch lengths sum to distance" {
    const allocator = std.testing.allocator;
    const d = 1.0;
    const dist = [_]f64{ 0.0, d, d, 0.0 };
    const names = [_][]const u8{ "X", "Y" };
    var tree = try upgma(allocator, &dist, 2, &names);
    defer tree.deinit();
    // Each leaf's branch length should be d/2
    try std.testing.expectApproxEqAbs(d / 2.0, tree.branch_length[0], 1e-10);
    try std.testing.expectApproxEqAbs(d / 2.0, tree.branch_length[1], 1e-10);
}

test "readNewick: simple 2-leaf tree" {
    const allocator = std.testing.allocator;
    var tree = try readNewick(allocator, "(A:0.1,B:0.2);");
    defer tree.deinit();

    try std.testing.expectEqual(@as(usize, 2), tree.n_leaves);
    try std.testing.expectEqual(@as(usize, 3), tree.n_nodes);
    try std.testing.expectEqualStrings("A", tree.names[0].?);
    try std.testing.expectEqualStrings("B", tree.names[1].?);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), tree.branch_length[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), tree.branch_length[1], 1e-10);
}

test "readNewick: nested tree" {
    const allocator = std.testing.allocator;
    var tree = try readNewick(allocator, "((A:0.1,B:0.2):0.3,C:0.4);");
    defer tree.deinit();

    try std.testing.expectEqual(@as(usize, 3), tree.n_leaves);
    try std.testing.expectEqual(@as(usize, 5), tree.n_nodes);
    try std.testing.expectEqualStrings("A", tree.names[0].?);
    try std.testing.expectEqualStrings("B", tree.names[1].?);
    try std.testing.expectEqualStrings("C", tree.names[2].?);
}

test "readNewick: round trip with writeNewick" {
    const allocator = std.testing.allocator;
    var tree = try readNewick(allocator, "(A:0.100000,B:0.200000);");
    defer tree.deinit();

    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try writeNewick(tree, buf.writer(allocator).any());
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "A:") != null);
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "B:") != null);
}

test "singleLinkage: 3 sequences" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.4, 0.8,
        0.4, 0.0, 0.8,
        0.8, 0.8, 0.0,
    };
    const names = [_][]const u8{ "A", "B", "C" };
    var tree = try singleLinkage(allocator, &dist, 3, &names);
    defer tree.deinit();

    try std.testing.expectEqual(@as(usize, 5), tree.n_nodes);
    try std.testing.expectEqual(@as(usize, 3), tree.n_leaves);
    // A-B merged first at min distance 0.4
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), tree.branch_length[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), tree.branch_length[1], 1e-10);
}

test "completeLinkage: 3 sequences" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.4, 0.8,
        0.4, 0.0, 0.8,
        0.8, 0.8, 0.0,
    };
    const names = [_][]const u8{ "A", "B", "C" };
    var tree = try completeLinkage(allocator, &dist, 3, &names);
    defer tree.deinit();

    try std.testing.expectEqual(@as(usize, 5), tree.n_nodes);
    // With complete linkage, distance from {A,B} to C = max(0.8, 0.8) = 0.8
    // Same as single linkage for this symmetric case
}

test "wpgma: 3 sequences" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.4, 0.8,
        0.4, 0.0, 0.6,
        0.8, 0.6, 0.0,
    };
    const names = [_][]const u8{ "A", "B", "C" };
    var tree = try wpgma(allocator, &dist, 3, &names);
    defer tree.deinit();

    try std.testing.expectEqual(@as(usize, 5), tree.n_nodes);
    // WPGMA: distance from {A,B} to C = (d(A,C) + d(B,C)) / 2 = (0.8 + 0.6) / 2 = 0.7
}

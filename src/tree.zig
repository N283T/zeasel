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

    /// Compare two trees for topological equivalence.
    /// Two trees are equal if they have the same number of leaves, the same
    /// leaf names, and the same branching pattern (same parent-child
    /// relationships when leaves are matched by name). Branch lengths are
    /// compared within the given tolerance.
    /// Uses the Goodman M(g) tree-mapping function via post-order traversal
    /// (see Easel esl_tree_Compare / Zmasek-Eddy 2001 SDI algorithm).
    /// Compare two trees for topological equivalence (heap-allocated version).
    /// Returns error.OutOfMemory if allocation fails, or true/false for match.
    pub fn compare(self: Tree, other: Tree, tolerance: f64) !bool {
        if (self.n_leaves != other.n_leaves) return false;
        if (self.n_nodes != other.n_nodes) return false;

        const n = self.n_leaves;
        const alloc = self.allocator;

        // Build taxon mapping: for each leaf in self, find the matching leaf
        // in other by name. Both trees must have leaf names.
        const taxa_map = try alloc.alloc(usize, n);
        defer alloc.free(taxa_map);
        for (0..n) |a| {
            const name_a = self.names[a] orelse return false;
            var found: ?usize = null;
            for (0..n) |b| {
                const name_b = other.names[b] orelse return false;
                if (std.mem.eql(u8, name_a, name_b)) {
                    found = b;
                    break;
                }
            }
            taxa_map[a] = found orelse return false;
        }

        // mg maps each node in self (leaf or internal) to the corresponding
        // node in other. For leaves, use taxa_map directly. For internal nodes,
        // compute via post-order traversal using the SDI algorithm.
        const mg = try alloc.alloc(i32, self.n_nodes);
        defer alloc.free(mg);

        // Initialize leaf mappings.
        for (0..n) |i| mg[i] = @intCast(taxa_map[i]);

        // Post-order traversal using a stack.
        const order = try alloc.alloc(usize, self.n_nodes);
        defer alloc.free(order);
        var order_len: usize = 0;
        {
            const stack = try alloc.alloc(usize, self.n_nodes);
            defer alloc.free(stack);
            const visited = try alloc.alloc(bool, self.n_nodes);
            defer alloc.free(visited);
            @memset(visited, false);
            // Find root.
            var root: usize = 0;
            for (0..self.n_nodes) |i| {
                if (self.parent[i] == -1) { root = i; break; }
            }
            var sp: usize = 0;
            stack[sp] = root;
            sp += 1;
            while (sp > 0) {
                const node = stack[sp - 1];
                const lc = self.left[node];
                const rc = self.right[node];
                const lc_done = (lc < 0 or visited[@intCast(lc)]);
                const rc_done = (rc < 0 or visited[@intCast(rc)]);
                if (lc_done and rc_done) {
                    order[order_len] = node;
                    order_len += 1;
                    visited[node] = true;
                    sp -= 1;
                } else {
                    if (!rc_done) { stack[sp] = @intCast(rc); sp += 1; }
                    if (!lc_done) { stack[sp] = @intCast(lc); sp += 1; }
                }
            }
        }

        // Process internal nodes in post-order.
        for (order[0..order_len]) |node_idx| {
            if (node_idx < n) continue; // skip leaves
            const lc: usize = @intCast(self.left[node_idx]);
            const rc: usize = @intCast(self.right[node_idx]);

            // Find parent of each child's mapped node in other.
            const a: i32 = if (mg[lc] >= 0) other.parent[@intCast(mg[lc])] else mg[lc];
            const b: i32 = if (mg[rc] >= 0) other.parent[@intCast(mg[rc])] else mg[rc];

            if (a != b) return false;
            mg[node_idx] = a; // mg[node] = corresponding internal node in other
        }

        // Compare branch lengths using the mapping.
        for (0..self.n_nodes) |i| {
            if (self.parent[i] < 0) continue; // root, no branch to compare
            const mapped = mg[i];
            if (mapped < 0) continue; // mapped to root in other
            const other_i: usize = @intCast(mapped);
            if (other.parent[other_i] < 0) continue; // root in other
            if (@abs(self.branch_length[i] - other.branch_length[other_i]) > tolerance) return false;
        }

        return true;
    }

    /// Validate internal consistency of the tree.
    /// Checks: parent/child links are reciprocal, all nodes reachable from
    /// root, no cycles, correct n_nodes/n_leaves counts, non-negative branch
    /// lengths.
    pub fn validate(self: Tree) bool {
        if (self.n_nodes == 0) return false;
        if (self.n_leaves == 0) return false;
        if (self.n_nodes != 2 * self.n_leaves - 1) return false;

        if (self.parent.len != self.n_nodes) return false;
        if (self.left.len != self.n_nodes) return false;
        if (self.right.len != self.n_nodes) return false;
        if (self.branch_length.len != self.n_nodes) return false;
        if (self.names.len != self.n_nodes) return false;

        // Find root: exactly one node should have parent == -1.
        var root_count: usize = 0;
        var root_idx: usize = 0;
        for (0..self.n_nodes) |i| {
            if (self.parent[i] == -1) {
                root_count += 1;
                root_idx = i;
            }
        }
        if (root_count != 1) return false;

        // Check parent-child reciprocity for all non-root nodes.
        for (0..self.n_nodes) |i| {
            if (self.parent[i] < 0) continue;
            const p: usize = @intCast(self.parent[i]);
            if (p >= self.n_nodes) return false;
            const is_left = (self.left[p] >= 0 and @as(usize, @intCast(self.left[p])) == i);
            const is_right = (self.right[p] >= 0 and @as(usize, @intCast(self.right[p])) == i);
            if (!is_left and !is_right) return false;
        }

        // Leaves have no children; internal nodes have two.
        for (0..self.n_leaves) |i| {
            if (self.left[i] != -1 or self.right[i] != -1) return false;
        }
        for (self.n_leaves..self.n_nodes) |i| {
            if (self.left[i] < 0 or self.right[i] < 0) return false;
            const lc: usize = @intCast(self.left[i]);
            const rc: usize = @intCast(self.right[i]);
            if (lc >= self.n_nodes or rc >= self.n_nodes) return false;
        }

        // BFS reachability from root and cycle detection.
        var visited: [512]bool = undefined;
        if (self.n_nodes > 512) return false;
        for (0..self.n_nodes) |i| visited[i] = false;

        var queue: [512]usize = undefined;
        var head: usize = 0;
        var tail: usize = 0;
        queue[tail] = root_idx;
        tail += 1;
        visited[root_idx] = true;

        while (head < tail) {
            const node = queue[head];
            head += 1;
            if (self.left[node] >= 0) {
                const child: usize = @intCast(self.left[node]);
                if (visited[child]) return false;
                visited[child] = true;
                queue[tail] = child;
                tail += 1;
            }
            if (self.right[node] >= 0) {
                const child: usize = @intCast(self.right[node]);
                if (visited[child]) return false;
                visited[child] = true;
                queue[tail] = child;
                tail += 1;
            }
        }

        for (0..self.n_nodes) |i| {
            if (!visited[i]) return false;
        }

        // Branch lengths must be non-negative.
        for (0..self.n_nodes) |i| {
            if (self.branch_length[i] < 0.0) return false;
        }

        return true;
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
                    bl[cn] = std.fmt.parseFloat(f64, input[start..i]) catch return error.InvalidFormat;
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

        // For linkage trees (single/complete), height = minD (raw linkage value).
        // For additive trees (WPGMA/UPGMA), height = minD / 2.
        // Reference: Easel esl_tree.c:1648-1649.
        const node_height = switch (linkage) {
            .single, .complete => min_d,
            .wpgma => min_d / 2.0,
        };
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

/// Generate a random binary tree by random sequential addition.
/// Each new leaf is attached at a random existing branch, splitting it
/// into two segments. Branch lengths are drawn from an exponential
/// distribution (uniform random, so the tree is not ultrametric).
///
/// `n_leaves` must be >= 2. Leaf names are "leaf0", "leaf1", etc.
pub fn simulate(allocator: Allocator, n_leaves: usize, rng: *Random) !Tree {
    if (n_leaves < 2) return error.InvalidInput;

    const total_nodes = 2 * n_leaves - 1;

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

    // Generate leaf names.
    for (0..n_leaves) |i| {
        names_arr[i] = try std.fmt.allocPrint(allocator, "leaf{d}", .{i});
    }

    // Start with a tree of 2 leaves joined by an internal node.
    const first_internal = n_leaves; // index of first internal node
    left_arr[first_internal] = 0;
    right_arr[first_internal] = 1;
    parent_arr[0] = @intCast(first_internal);
    parent_arr[1] = @intCast(first_internal);
    bl[0] = rng.uniform() + 0.01;
    bl[1] = rng.uniform() + 0.01;

    var next_internal = first_internal + 1;

    // Sequentially add remaining leaves.
    for (2..n_leaves) |leaf| {
        // Pick a random existing node (excluding the current root) to split.
        // Existing nodes are 0..(leaf-1) leaves + first_internal..(next_internal-1) internals.
        // We pick from all non-root nodes (those with parent != -1).
        // Actually simpler: pick a random edge. An edge is identified by its
        // child node (every node except root has exactly one incoming edge).
        // Count non-root nodes so far.
        const n_edges = leaf + (next_internal - first_internal) - 1;
        const pick = rng.uniformInt(@intCast(n_edges));

        // Enumerate non-root nodes to find the picked one.
        var target: usize = 0;
        var count: u32 = 0;
        for (0..next_internal + 1) |node| {
            if (node >= total_nodes) continue;
            if (parent_arr[node] == -1) continue;
            if (count == pick) {
                target = node;
                break;
            }
            count += 1;
        }

        // Insert a new internal node on the edge from target to its parent.
        const new_internal = next_internal;
        const old_parent = parent_arr[target];
        const old_parent_u: usize = @intCast(old_parent);

        // Replace target in old_parent's children with new_internal.
        if (left_arr[old_parent_u] == @as(i32, @intCast(target))) {
            left_arr[old_parent_u] = @intCast(new_internal);
        } else {
            right_arr[old_parent_u] = @intCast(new_internal);
        }

        parent_arr[new_internal] = old_parent;
        left_arr[new_internal] = @intCast(target);
        right_arr[new_internal] = @intCast(leaf);
        parent_arr[target] = @intCast(new_internal);
        parent_arr[leaf] = @intCast(new_internal);

        // Split the original branch length.
        const orig_bl = bl[target];
        const frac = rng.uniform();
        bl[new_internal] = orig_bl * frac;
        bl[target] = orig_bl * (1.0 - frac);
        bl[leaf] = rng.uniform() + 0.01;

        next_internal += 1;
    }

    return Tree{
        .parent = parent_arr,
        .left = left_arr,
        .right = right_arr,
        .branch_length = bl,
        .names = names_arr,
        .n_nodes = total_nodes,
        .n_leaves = n_leaves,
        .allocator = allocator,
    };
}

const Random = @import("util/random.zig").Random;

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

test "readNewick: invalid branch length returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(error.InvalidFormat, readNewick(allocator, "(A:abc,B:0.2);"));
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
    // A-B merged first at min distance 0.4.
    // For linkage trees, height = minD (not minD/2), so branch_length = 0.4 - 0 = 0.4.
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), tree.branch_length[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), tree.branch_length[1], 1e-10);
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

test "validate: valid UPGMA tree passes" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.4, 0.8,
        0.4, 0.0, 0.8,
        0.8, 0.8, 0.0,
    };
    const names = [_][]const u8{ "A", "B", "C" };
    var tree = try upgma(allocator, &dist, 3, &names);
    defer tree.deinit();
    try std.testing.expect(tree.validate());
}

test "validate: valid parsed Newick tree passes" {
    const allocator = std.testing.allocator;
    var tree = try readNewick(allocator, "((A:0.1,B:0.2):0.3,C:0.4);");
    defer tree.deinit();
    try std.testing.expect(tree.validate());
}

test "validate: detects broken parent link" {
    const allocator = std.testing.allocator;
    var tree = try readNewick(allocator, "((A:0.1,B:0.2):0.3,C:0.4);");
    defer tree.deinit();
    tree.parent[0] = -1;
    try std.testing.expect(!tree.validate());
}

test "compare: identical trees are equal" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.4, 0.8,
        0.4, 0.0, 0.8,
        0.8, 0.8, 0.0,
    };
    const names = [_][]const u8{ "A", "B", "C" };
    var t1 = try upgma(allocator, &dist, 3, &names);
    defer t1.deinit();
    var t2 = try upgma(allocator, &dist, 3, &names);
    defer t2.deinit();
    try std.testing.expect(try t1.compare(t2, 1e-6));
}

test "compare: different topologies are not equal" {
    const allocator = std.testing.allocator;
    var t1 = try readNewick(allocator, "((A:0.1,B:0.2):0.3,C:0.4);");
    defer t1.deinit();
    var t2 = try readNewick(allocator, "((A:0.1,C:0.2):0.3,B:0.4);");
    defer t2.deinit();
    try std.testing.expect(!try t1.compare(t2, 1e-6));
}

test "compare: same topology different branch lengths fails strict tolerance" {
    const allocator = std.testing.allocator;
    var t1 = try readNewick(allocator, "((A:0.1,B:0.2):0.3,C:0.4);");
    defer t1.deinit();
    var t2 = try readNewick(allocator, "((A:0.5,B:0.2):0.3,C:0.4);");
    defer t2.deinit();
    try std.testing.expect(!try t1.compare(t2, 1e-6));
    try std.testing.expect(try t1.compare(t2, 1.0));
}

test "compare: round-trip Newick preserves topology" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.6,
        0.6, 0.0,
    };
    const names = [_][]const u8{ "X", "Y" };
    var t1 = try upgma(allocator, &dist, 2, &names);
    defer t1.deinit();

    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try writeNewick(t1, buf.writer(allocator).any());

    var t2 = try readNewick(allocator, buf.items);
    defer t2.deinit();
    try std.testing.expect(try t1.compare(t2, 1e-6));
}

test "simulate: random tree structure is valid" {
    const allocator = std.testing.allocator;
    var rng = Random.init(42);
    var tree = try simulate(allocator, 10, &rng);
    defer tree.deinit();

    try std.testing.expectEqual(@as(usize, 10), tree.n_leaves);
    try std.testing.expectEqual(@as(usize, 19), tree.n_nodes);
    try std.testing.expect(tree.validate());
}

test "simulate: minimum 2 leaves" {
    const allocator = std.testing.allocator;
    var rng = Random.init(1);
    var tree = try simulate(allocator, 2, &rng);
    defer tree.deinit();

    try std.testing.expectEqual(@as(usize, 2), tree.n_leaves);
    try std.testing.expectEqual(@as(usize, 3), tree.n_nodes);
    try std.testing.expect(tree.validate());
}

test "simulate: error on fewer than 2 leaves" {
    var rng = Random.init(1);
    try std.testing.expectError(error.InvalidInput, simulate(std.testing.allocator, 1, &rng));
    try std.testing.expectError(error.InvalidInput, simulate(std.testing.allocator, 0, &rng));
}

// Red-black tree — a self-balancing binary search tree with f64 keys.
//
// Used by HMMER's sparse DP matrix management for efficient
// ordered key lookup, insertion, and deletion.

const std = @import("std");
const Allocator = std.mem.Allocator;

const Color = enum { red, black };

const Node = struct {
    key: f64,
    value: usize,
    color: Color,
    left: ?*Node,
    right: ?*Node,
    parent: ?*Node,
};

pub const RedBlackTree = struct {
    root: ?*Node,
    count: usize,
    allocator: Allocator,

    pub fn init(allocator: Allocator) RedBlackTree {
        return .{ .root = null, .count = 0, .allocator = allocator };
    }

    /// Insert a key-value pair. Replaces value if key exists.
    pub fn insert(self: *RedBlackTree, key: f64, value: usize) !void {
        const new_node = try self.allocator.create(Node);
        new_node.* = .{
            .key = key,
            .value = value,
            .color = .red,
            .left = null,
            .right = null,
            .parent = null,
        };

        if (self.root == null) {
            new_node.color = .black;
            self.root = new_node;
            self.count += 1;
            return;
        }

        // BST insert
        var current = self.root;
        var par: ?*Node = null;
        while (current) |c| {
            par = c;
            if (key < c.key) {
                current = c.left;
            } else if (key > c.key) {
                current = c.right;
            } else {
                // Key exists — update value, free new node
                c.value = value;
                self.allocator.destroy(new_node);
                return;
            }
        }

        new_node.parent = par;
        if (par) |p| {
            if (key < p.key) {
                p.left = new_node;
            } else {
                p.right = new_node;
            }
        }
        self.count += 1;

        // Fix red-black properties
        self.insertFixup(new_node);
    }

    /// Search for a key. Returns the value or null.
    pub fn search(self: RedBlackTree, key: f64) ?usize {
        var current = self.root;
        while (current) |c| {
            if (key < c.key) {
                current = c.left;
            } else if (key > c.key) {
                current = c.right;
            } else {
                return c.value;
            }
        }
        return null;
    }

    /// Return the minimum key in the tree, or null if empty.
    pub fn min(self: RedBlackTree) ?f64 {
        var current = self.root orelse return null;
        while (current.left) |l| current = l;
        return current.key;
    }

    /// Return the maximum key in the tree, or null if empty.
    pub fn max(self: RedBlackTree) ?f64 {
        var current = self.root orelse return null;
        while (current.right) |r| current = r;
        return current.key;
    }

    /// Free all nodes.
    pub fn deinit(self: *RedBlackTree) void {
        freeSubtree(self.allocator, self.root);
        self.root = null;
        self.count = 0;
    }

    fn freeSubtree(allocator: Allocator, node: ?*Node) void {
        const n = node orelse return;
        freeSubtree(allocator, n.left);
        freeSubtree(allocator, n.right);
        allocator.destroy(n);
    }

    fn insertFixup(self: *RedBlackTree, z_in: *Node) void {
        var z = z_in;
        while (z.parent) |p| {
            if (p.color != .red) break;
            const gp = p.parent orelse break;

            if (p == gp.left) {
                const uncle = gp.right;
                if (uncle != null and uncle.?.color == .red) {
                    p.color = .black;
                    uncle.?.color = .black;
                    gp.color = .red;
                    z = gp;
                } else {
                    if (z == p.right) {
                        z = p;
                        self.rotateLeft(z);
                    }
                    if (z.parent) |zp| {
                        zp.color = .black;
                        if (zp.parent) |zgp| {
                            zgp.color = .red;
                            self.rotateRight(zgp);
                        }
                    }
                }
            } else {
                const uncle = gp.left;
                if (uncle != null and uncle.?.color == .red) {
                    p.color = .black;
                    uncle.?.color = .black;
                    gp.color = .red;
                    z = gp;
                } else {
                    if (z == p.left) {
                        z = p;
                        self.rotateRight(z);
                    }
                    if (z.parent) |zp| {
                        zp.color = .black;
                        if (zp.parent) |zgp| {
                            zgp.color = .red;
                            self.rotateLeft(zgp);
                        }
                    }
                }
            }
        }
        if (self.root) |r| r.color = .black;
    }

    fn rotateLeft(self: *RedBlackTree, x: *Node) void {
        const y = x.right orelse return;
        x.right = y.left;
        if (y.left) |yl| yl.parent = x;
        y.parent = x.parent;
        if (x.parent == null) {
            self.root = y;
        } else if (x.parent) |xp| {
            if (x == xp.left) {
                xp.left = y;
            } else {
                xp.right = y;
            }
        }
        y.left = x;
        x.parent = y;
    }

    fn rotateRight(self: *RedBlackTree, y: *Node) void {
        const x = y.left orelse return;
        y.left = x.right;
        if (x.right) |xr| xr.parent = y;
        x.parent = y.parent;
        if (y.parent == null) {
            self.root = x;
        } else if (y.parent) |yp| {
            if (y == yp.left) {
                yp.left = x;
            } else {
                yp.right = x;
            }
        }
        x.right = y;
        y.parent = x;
    }
};

// --- Tests ---

test "insert and search" {
    const allocator = std.testing.allocator;
    var tree = RedBlackTree.init(allocator);
    defer tree.deinit();

    try tree.insert(5.0, 50);
    try tree.insert(3.0, 30);
    try tree.insert(7.0, 70);
    try tree.insert(1.0, 10);

    try std.testing.expectEqual(@as(?usize, 50), tree.search(5.0));
    try std.testing.expectEqual(@as(?usize, 30), tree.search(3.0));
    try std.testing.expectEqual(@as(?usize, 70), tree.search(7.0));
    try std.testing.expectEqual(@as(?usize, null), tree.search(99.0));
    try std.testing.expectEqual(@as(usize, 4), tree.count);
}

test "min and max" {
    const allocator = std.testing.allocator;
    var tree = RedBlackTree.init(allocator);
    defer tree.deinit();

    try tree.insert(5.0, 0);
    try tree.insert(2.0, 0);
    try tree.insert(8.0, 0);
    try tree.insert(1.0, 0);
    try tree.insert(9.0, 0);

    try std.testing.expectEqual(@as(?f64, 1.0), tree.min());
    try std.testing.expectEqual(@as(?f64, 9.0), tree.max());
}

test "insert duplicate updates value" {
    const allocator = std.testing.allocator;
    var tree = RedBlackTree.init(allocator);
    defer tree.deinit();

    try tree.insert(5.0, 10);
    try tree.insert(5.0, 20);

    try std.testing.expectEqual(@as(?usize, 20), tree.search(5.0));
    try std.testing.expectEqual(@as(usize, 1), tree.count);
}

test "empty tree" {
    const allocator = std.testing.allocator;
    var tree = RedBlackTree.init(allocator);
    defer tree.deinit();

    try std.testing.expectEqual(@as(?usize, null), tree.search(1.0));
    try std.testing.expectEqual(@as(?f64, null), tree.min());
    try std.testing.expectEqual(@as(?f64, null), tree.max());
}

test "many inserts maintain balance" {
    const allocator = std.testing.allocator;
    var tree = RedBlackTree.init(allocator);
    defer tree.deinit();

    // Insert 100 sequential keys — should not stack overflow (balanced)
    for (0..100) |i| {
        try tree.insert(@floatFromInt(i), i);
    }
    try std.testing.expectEqual(@as(usize, 100), tree.count);
    try std.testing.expectEqual(@as(?usize, 0), tree.search(0.0));
    try std.testing.expectEqual(@as(?usize, 99), tree.search(99.0));
}

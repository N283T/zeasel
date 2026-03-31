// Huffman coding for digitized alphabets.
//
// Builds a canonical Huffman code from symbol frequencies and provides
// encode/decode operations. Used for compressed sequence storage.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const HuffmanCode = struct {
    /// Code word for each symbol (bit pattern stored in u32).
    codes: []u32,
    /// Code length (bits) for each symbol.
    lengths: []u5,
    /// Number of symbols.
    n_symbols: usize,
    allocator: Allocator,

    /// Build a Huffman code from symbol frequencies.
    /// freq[i] = frequency of symbol i, for i in 0..n_symbols-1.
    pub fn build(allocator: Allocator, freq: []const f64) !HuffmanCode {
        const n = freq.len;
        if (n == 0) return error.InvalidInput;
        if (n == 1) {
            const codes = try allocator.alloc(u32, 1);
            codes[0] = 0;
            const lengths = try allocator.alloc(u5, 1);
            lengths[0] = 1;
            return .{ .codes = codes, .lengths = lengths, .n_symbols = 1, .allocator = allocator };
        }

        // Build code lengths using the package-merge algorithm approximation:
        // sort by frequency, assign lengths greedily.
        const indices = try allocator.alloc(usize, n);
        defer allocator.free(indices);
        for (0..n) |i| indices[i] = i;

        // Sort indices by frequency (ascending)
        std.mem.sort(usize, indices, freq, struct {
            fn lessThan(f: []const f64, a: usize, b: usize) bool {
                return f[a] < f[b];
            }
        }.lessThan);

        // Simple length assignment: use a binary tree simulation
        const lengths = try allocator.alloc(u5, n);
        errdefer allocator.free(lengths);
        @memset(lengths, 0);

        // Build tree bottom-up to determine lengths
        const tree_size = 2 * n - 1;
        const tree_freq = try allocator.alloc(f64, tree_size);
        defer allocator.free(tree_freq);
        const tree_left = try allocator.alloc(?usize, tree_size);
        defer allocator.free(tree_left);
        const tree_right = try allocator.alloc(?usize, tree_size);
        defer allocator.free(tree_right);

        // Initialize leaves
        for (0..n) |i| {
            tree_freq[i] = freq[indices[i]];
            tree_left[i] = null;
            tree_right[i] = null;
        }

        // Priority queue using sorted merge
        var active = try allocator.alloc(bool, tree_size);
        defer allocator.free(active);
        @memset(active, false);
        for (0..n) |i| active[i] = true;

        var next_internal: usize = n;
        for (0..n - 1) |_| {
            // Find two smallest active nodes
            var m1: usize = tree_size;
            var m2: usize = tree_size;
            var m1_freq: f64 = std.math.inf(f64);
            var m2_freq: f64 = std.math.inf(f64);

            for (0..next_internal) |i| {
                if (!active[i]) continue;
                if (tree_freq[i] < m1_freq) {
                    m2 = m1;
                    m2_freq = m1_freq;
                    m1 = i;
                    m1_freq = tree_freq[i];
                } else if (tree_freq[i] < m2_freq) {
                    m2 = i;
                    m2_freq = tree_freq[i];
                }
            }

            // Create internal node
            tree_freq[next_internal] = m1_freq + m2_freq;
            tree_left[next_internal] = m1;
            tree_right[next_internal] = m2;
            active[next_internal] = true;
            active[m1] = false;
            active[m2] = false;
            next_internal += 1;
        }

        // Traverse tree to assign lengths
        const depth = try allocator.alloc(u5, tree_size);
        defer allocator.free(depth);
        @memset(depth, 0);

        // BFS from root (last created node)
        var ii: usize = next_internal;
        while (ii > 0) {
            ii -= 1;
            if (tree_left[ii]) |l| depth[l] = depth[ii] +| 1;
            if (tree_right[ii]) |r| depth[r] = depth[ii] +| 1;
        }

        // Map leaf depths back to original symbols
        for (0..n) |i| {
            lengths[indices[i]] = if (depth[i] > 0) depth[i] else 1;
        }

        // Generate canonical codes from lengths
        const codes = try allocator.alloc(u32, n);
        errdefer allocator.free(codes);
        @memset(codes, 0);

        // Sort by length then index for canonical assignment
        const code_order = try allocator.alloc(usize, n);
        defer allocator.free(code_order);
        for (0..n) |i| code_order[i] = i;
        std.mem.sort(usize, code_order, lengths, struct {
            fn lessThan(l: []const u5, a: usize, b: usize) bool {
                if (l[a] != l[b]) return l[a] < l[b];
                return a < b;
            }
        }.lessThan);

        var code: u32 = 0;
        var prev_len: u5 = lengths[code_order[0]];
        codes[code_order[0]] = 0;
        for (1..n) |i| {
            code += 1;
            const curr_len = lengths[code_order[i]];
            if (curr_len > prev_len) {
                code <<= @intCast(curr_len - prev_len);
            }
            codes[code_order[i]] = code;
            prev_len = curr_len;
        }

        return .{ .codes = codes, .lengths = lengths, .n_symbols = n, .allocator = allocator };
    }

    /// Average code length (bits per symbol) given frequencies.
    pub fn avgLength(self: HuffmanCode, freq: []const f64) f64 {
        var total_freq: f64 = 0;
        var weighted: f64 = 0;
        for (0..self.n_symbols) |i| {
            total_freq += freq[i];
            weighted += freq[i] * @as(f64, @floatFromInt(self.lengths[i]));
        }
        if (total_freq == 0) return 0;
        return weighted / total_freq;
    }

    pub fn deinit(self: *HuffmanCode) void {
        self.allocator.free(self.codes);
        self.allocator.free(self.lengths);
    }
};

// --- Tests ---

test "build: uniform frequencies give equal-ish lengths" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 1.0, 1.0, 1.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    try std.testing.expectEqual(@as(usize, 4), hc.n_symbols);
    // All lengths should be 2 for 4 uniform symbols
    for (hc.lengths) |l| {
        try std.testing.expectEqual(@as(u5, 2), l);
    }
}

test "build: skewed frequencies give different lengths" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 100.0, 1.0, 1.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    // Most frequent symbol should have shortest code
    try std.testing.expect(hc.lengths[0] <= hc.lengths[1]);
    try std.testing.expect(hc.lengths[0] <= hc.lengths[2]);
    try std.testing.expect(hc.lengths[0] <= hc.lengths[3]);
}

test "build: 2 symbols" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 3.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    try std.testing.expectEqual(@as(u5, 1), hc.lengths[0]);
    try std.testing.expectEqual(@as(u5, 1), hc.lengths[1]);
    // Codes should be 0 and 1
    try std.testing.expect(hc.codes[0] != hc.codes[1]);
}

test "avgLength: uniform 4 symbols = 2 bits" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 1.0, 1.0, 1.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, 2.0), hc.avgLength(&freq), 1e-9);
}

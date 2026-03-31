/// Histogram with fixed-width bins for score distribution fitting.
///
/// Supports accumulating values into bins, computing basic statistics
/// (mean, variance, stddev), and chi-squared goodness-of-fit testing.
const std = @import("std");
const Allocator = std.mem.Allocator;
const math = std.math;

pub const Histogram = struct {
    /// Per-bin counts; counts[i] covers [edges[i], edges[i+1]).
    counts: []u64,
    /// Bin edges; length is n_bins + 1.
    edges: []f64,
    n_bins: usize,
    /// Running statistics (include out-of-range values).
    total: u64,
    min_val: f64,
    max_val: f64,
    sum: f64,
    sum_sq: f64,
    allocator: Allocator,

    /// Create a histogram with n_bins fixed-width bins covering [min_val, max_val).
    pub fn init(allocator: Allocator, min_val: f64, max_val: f64, n_bins: usize) !Histogram {
        const counts = try allocator.alloc(u64, n_bins);
        @memset(counts, 0);
        const edges = try allocator.alloc(f64, n_bins + 1);
        const width = (max_val - min_val) / @as(f64, @floatFromInt(n_bins));
        for (0..n_bins + 1) |i| {
            edges[i] = min_val + @as(f64, @floatFromInt(i)) * width;
        }
        return Histogram{
            .counts = counts,
            .edges = edges,
            .n_bins = n_bins,
            .total = 0,
            .min_val = math.inf(f64),
            .max_val = -math.inf(f64),
            .sum = 0,
            .sum_sq = 0,
            .allocator = allocator,
        };
    }

    /// Add a single value. Updates running statistics unconditionally;
    /// only increments a bin count if x is within [edges[0], edges[n_bins]).
    pub fn add(self: *Histogram, x: f64) void {
        self.total += 1;
        self.sum += x;
        self.sum_sq += x * x;
        if (x < self.min_val) self.min_val = x;
        if (x > self.max_val) self.max_val = x;

        const bin_width = (self.edges[self.n_bins] - self.edges[0]) / @as(f64, @floatFromInt(self.n_bins));
        if (bin_width <= 0) return;
        const idx_f = (x - self.edges[0]) / bin_width;
        if (idx_f < 0) return;
        const idx: usize = @intFromFloat(idx_f);
        if (idx >= self.n_bins) return;
        self.counts[idx] += 1;
    }

    /// Add a slice of values.
    pub fn addAll(self: *Histogram, values: []const f64) void {
        for (values) |v| self.add(v);
    }

    /// Mean of all added values. Returns 0 when no values have been added.
    pub fn mean(self: Histogram) f64 {
        if (self.total == 0) return 0;
        return self.sum / @as(f64, @floatFromInt(self.total));
    }

    /// Sample variance (Bessel-corrected) of all added values.
    /// Returns 0 when fewer than 2 values have been added.
    pub fn variance(self: Histogram) f64 {
        if (self.total < 2) return 0;
        const n = @as(f64, @floatFromInt(self.total));
        const m = self.sum / n;
        return (self.sum_sq - n * m * m) / (n - 1.0);
    }

    /// Standard deviation of all added values.
    pub fn stddev(self: Histogram) f64 {
        return math.sqrt(self.variance());
    }

    /// Chi-squared goodness-of-fit statistic against expected bin frequencies.
    /// expected.len must equal n_bins. Bins with expected[i] == 0 are skipped.
    pub fn chiSquared(self: Histogram, expected: []const f64) f64 {
        std.debug.assert(expected.len == self.n_bins);
        var chi2: f64 = 0;
        for (0..self.n_bins) |i| {
            if (expected[i] > 0) {
                const diff = @as(f64, @floatFromInt(self.counts[i])) - expected[i];
                chi2 += (diff * diff) / expected[i];
            }
        }
        return chi2;
    }

    pub fn deinit(self: *Histogram) void {
        self.allocator.free(self.counts);
        self.allocator.free(self.edges);
    }
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "one value per bin" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // 0.5, 1.5, ..., 9.5 — each falls into a different bin.
    var i: usize = 0;
    while (i < 10) : (i += 1) {
        h.add(@as(f64, @floatFromInt(i)) + 0.5);
    }

    try std.testing.expectEqual(@as(u64, 10), h.total);
    for (h.counts) |c| {
        try std.testing.expectEqual(@as(u64, 1), c);
    }
}

test "all values in the same bin" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // All values in [2, 3): bin index 2.
    h.add(2.1);
    h.add(2.5);
    h.add(2.9);

    try std.testing.expectEqual(@as(u64, 3), h.total);
    try std.testing.expectEqual(@as(u64, 3), h.counts[2]);
    // Every other bin must be zero.
    for (h.counts, 0..) |c, idx| {
        if (idx != 2) try std.testing.expectEqual(@as(u64, 0), c);
    }
}

test "mean variance stddev of known dataset" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // Dataset: 2, 4, 4, 4, 5, 5, 7, 9
    // mean = 5, variance = 4 (population), sample variance ≈ 4.571
    const data = [_]f64{ 2, 4, 4, 4, 5, 5, 7, 9 };
    h.addAll(&data);

    const eps = 1e-10;
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), h.mean(), eps);

    const expected_var = 32.0 / 7.0; // sample variance
    try std.testing.expectApproxEqAbs(expected_var, h.variance(), eps);
    try std.testing.expectApproxEqAbs(math.sqrt(expected_var), h.stddev(), eps);
}

test "out-of-range values tracked in stats but not in bins" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    h.add(5.0); // in range
    h.add(-1.0); // below range
    h.add(10.0); // at upper edge — excluded by [edges[0], edges[n_bins])
    h.add(15.0); // above range

    try std.testing.expectEqual(@as(u64, 4), h.total);

    // Only the value 5.0 should be in a bin.
    var bin_total: u64 = 0;
    for (h.counts) |c| bin_total += c;
    try std.testing.expectEqual(@as(u64, 1), bin_total);
}

test "chiSquared with uniform expected" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // Perfectly uniform: one count per bin, expected = 1.0 each.
    var i: usize = 0;
    while (i < 10) : (i += 1) {
        h.add(@as(f64, @floatFromInt(i)) + 0.5);
    }

    var expected: [10]f64 = undefined;
    for (&expected) |*e| e.* = 1.0;

    try std.testing.expectApproxEqAbs(@as(f64, 0.0), h.chiSquared(&expected), 1e-10);
}

test "empty histogram" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    try std.testing.expectEqual(@as(u64, 0), h.total);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), h.mean(), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), h.variance(), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), h.stddev(), 1e-10);
}

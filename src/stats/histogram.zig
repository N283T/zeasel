/// Histogram with fixed-width bins for score distribution fitting.
///
/// Supports accumulating values into bins, computing basic statistics
/// (mean, variance, stddev), and chi-squared goodness-of-fit testing.
/// Bins are dynamically resized when values fall outside the current range.
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

    /// Tail fitting state. Once set, no more data can be added.
    is_done: bool,
    /// Tail cutoff: values <= phi are "unobserved" (virtual left-censoring).
    phi: f64,
    /// Number of censored observations (<= phi).
    n_censored: u64,
    /// Number of observed (tail) observations.
    n_observed: u64,
    /// Bin index of the censoring boundary.
    cmin: usize,

    /// Expected tail parameters.
    tail_lambda: f64,
    tail_mu: f64,
    has_expected_tail: bool,

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
            .is_done = false,
            .phi = -math.inf(f64),
            .n_censored = 0,
            .n_observed = 0,
            .cmin = 0,
            .tail_lambda = 0,
            .tail_mu = 0,
            .has_expected_tail = false,
        };
    }

    /// Bin width (constant across all bins).
    fn binWidth(self: Histogram) f64 {
        return (self.edges[self.n_bins] - self.edges[0]) / @as(f64, @floatFromInt(self.n_bins));
    }

    /// Convert a score to a bin index. May return negative (below range)
    /// or >= n_bins (above range).
    fn scoreToBin(self: Histogram, x: f64) i64 {
        const w = self.binWidth();
        if (w <= 0) return 0;
        return @intFromFloat(@floor((x - self.edges[0]) / w));
    }

    /// Grow the histogram bins so that bin index `target_bin` (possibly negative
    /// or >= n_bins) becomes valid. Follows Easel's strategy of 2x overallocation.
    fn growTo(self: *Histogram, target_bin: i64) !void {
        if (target_bin < 0) {
            // Expand below: add bins at the front.
            const deficit: usize = @intCast(-target_bin);
            const nnew = deficit * 2; // 2x overalloc
            const new_n_bins = self.n_bins + nnew;
            const w = self.binWidth();

            const new_counts = try self.allocator.alloc(u64, new_n_bins);
            @memset(new_counts[0..nnew], 0);
            @memcpy(new_counts[nnew..], self.counts);
            self.allocator.free(self.counts);
            self.counts = new_counts;

            const new_edges = try self.allocator.alloc(f64, new_n_bins + 1);
            const new_lo = self.edges[0] - @as(f64, @floatFromInt(nnew)) * w;
            for (0..new_n_bins + 1) |i| {
                new_edges[i] = new_lo + @as(f64, @floatFromInt(i)) * w;
            }
            self.allocator.free(self.edges);
            self.edges = new_edges;
            self.n_bins = new_n_bins;

            // Adjust cmin if tail is set.
            if (self.is_done) {
                self.cmin += nnew;
            }
        } else {
            // Expand above: add bins at the end.
            const target: usize = @intCast(target_bin);
            const deficit = target - self.n_bins + 1;
            const nnew = deficit * 2; // 2x overalloc
            const new_n_bins = self.n_bins + nnew;
            const w = self.binWidth();

            const new_counts = try self.allocator.alloc(u64, new_n_bins);
            @memcpy(new_counts[0..self.n_bins], self.counts);
            @memset(new_counts[self.n_bins..], 0);
            self.allocator.free(self.counts);
            self.counts = new_counts;

            const new_edges = try self.allocator.alloc(f64, new_n_bins + 1);
            for (0..new_n_bins + 1) |i| {
                new_edges[i] = self.edges[0] + @as(f64, @floatFromInt(i)) * w;
            }
            self.allocator.free(self.edges);
            self.edges = new_edges;
            self.n_bins = new_n_bins;
        }
    }

    /// Add a single value. Updates running statistics unconditionally.
    /// Dynamically resizes bins if x is outside the current range.
    pub fn add(self: *Histogram, x: f64) !void {
        if (self.is_done) return error.HistogramFinished;

        self.total += 1;
        self.sum += x;
        self.sum_sq += x * x;
        if (x < self.min_val) self.min_val = x;
        if (x > self.max_val) self.max_val = x;

        const w = self.binWidth();
        if (w <= 0) return;

        const bin_i = self.scoreToBin(x);

        // Resize if out of range.
        if (bin_i < 0 or bin_i >= @as(i64, @intCast(self.n_bins))) {
            try self.growTo(bin_i);
        }

        // Recompute after possible resize.
        const idx: usize = @intCast(self.scoreToBin(x));
        if (idx < self.n_bins) {
            self.counts[idx] += 1;
        }
    }

    /// Add a slice of values.
    pub fn addAll(self: *Histogram, values: []const f64) !void {
        for (values) |v| try self.add(v);
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

    /// Mark the tail cutoff point for fitting.
    /// Data points with values <= phi are treated as "unobserved" (virtual
    /// left-censoring). The actual phi is snapped to the nearest bin lower
    /// bound. Returns the fraction of data in the observed right tail.
    pub fn setTail(self: *Histogram, phi: f64) f64 {
        const bin_i = self.scoreToBin(phi);
        // Clamp to valid range.
        const cmin_idx: usize = if (bin_i < 0)
            0
        else if (bin_i >= @as(i64, @intCast(self.n_bins)))
            self.n_bins
        else
            @intCast(bin_i);

        // Snap phi to the lower bound of the cmin bin.
        self.cmin = cmin_idx;
        self.phi = self.edges[cmin_idx];

        // Count censored observations.
        var z: u64 = 0;
        for (0..cmin_idx) |b| {
            z += self.counts[b];
        }
        self.n_censored = z;
        self.n_observed = self.total - z;
        self.is_done = true;

        if (self.total == 0) return 0;
        return @as(f64, @floatFromInt(self.n_observed)) / @as(f64, @floatFromInt(self.total));
    }

    /// Set expected distribution parameters (exponential tail) for
    /// goodness-of-fit testing. The parameters describe the tail
    /// distribution: P(x) ~ lambda * exp(-lambda * (x - mu)).
    pub fn setExpectedTail(self: *Histogram, lambda: f64, mu: f64) void {
        self.tail_lambda = lambda;
        self.tail_mu = mu;
        self.has_expected_tail = true;
        self.is_done = true;
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
        try h.add(@as(f64, @floatFromInt(i)) + 0.5);
    }

    try std.testing.expectEqual(@as(u64, 10), h.total);
    for (h.counts[0..10]) |c| {
        try std.testing.expectEqual(@as(u64, 1), c);
    }
}

test "all values in the same bin" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // All values in [2, 3): bin index 2.
    try h.add(2.1);
    try h.add(2.5);
    try h.add(2.9);

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
    try h.addAll(&data);

    const eps = 1e-10;
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), h.mean(), eps);

    const expected_var = 32.0 / 7.0; // sample variance
    try std.testing.expectApproxEqAbs(expected_var, h.variance(), eps);
    try std.testing.expectApproxEqAbs(math.sqrt(expected_var), h.stddev(), eps);
}

test "out-of-range values dynamically resize bins" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    try h.add(5.0); // in range
    try h.add(-5.0); // below range -> triggers resize
    try h.add(15.0); // above range -> triggers resize

    try std.testing.expectEqual(@as(u64, 3), h.total);

    // All three values should be in bins now (after resize).
    var bin_total: u64 = 0;
    for (h.counts) |c| bin_total += c;
    try std.testing.expectEqual(@as(u64, 3), bin_total);

    // The histogram range should now include -5 and 15.
    try std.testing.expect(h.edges[0] <= -5.0);
    try std.testing.expect(h.edges[h.n_bins] >= 15.0);
}

test "dynamic resize below preserves existing counts" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // Add values to bins 0-4.
    for (0..5) |i| {
        try h.add(@as(f64, @floatFromInt(i)) + 0.5);
    }

    // Now add a value below range.
    try h.add(-3.0);

    try std.testing.expectEqual(@as(u64, 6), h.total);

    // All 6 values should be in bins.
    var bin_total: u64 = 0;
    for (h.counts) |c| bin_total += c;
    try std.testing.expectEqual(@as(u64, 6), bin_total);
}

test "dynamic resize above preserves existing counts" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    try h.add(0.5);
    try h.add(5.5);
    try h.add(20.0); // above range

    try std.testing.expectEqual(@as(u64, 3), h.total);
    var bin_total: u64 = 0;
    for (h.counts) |c| bin_total += c;
    try std.testing.expectEqual(@as(u64, 3), bin_total);
}

test "chiSquared with uniform expected" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // Perfectly uniform: one count per bin, expected = 1.0 each.
    var i: usize = 0;
    while (i < 10) : (i += 1) {
        try h.add(@as(f64, @floatFromInt(i)) + 0.5);
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

test "setTail: basic tail marking" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // Add 10 values: 0.5, 1.5, ..., 9.5
    for (0..10) |i| {
        try h.add(@as(f64, @floatFromInt(i)) + 0.5);
    }

    // Set tail at 5.0 -> bins [0..4] are censored, [5..9] are observed.
    const tail_mass = h.setTail(5.0);

    try std.testing.expectEqual(true, h.is_done);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), h.phi, 1e-10);
    try std.testing.expectEqual(@as(u64, 5), h.n_censored);
    try std.testing.expectEqual(@as(u64, 5), h.n_observed);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), tail_mass, 1e-10);

    // Cannot add more data after setTail.
    try std.testing.expectError(error.HistogramFinished, h.add(3.0));
}

test "setExpectedTail: stores parameters" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    for (0..10) |i| {
        try h.add(@as(f64, @floatFromInt(i)) + 0.5);
    }

    h.setExpectedTail(0.5, 3.0);

    try std.testing.expectEqual(true, h.has_expected_tail);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), h.tail_lambda, 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), h.tail_mu, 1e-10);
    try std.testing.expectEqual(true, h.is_done);
}

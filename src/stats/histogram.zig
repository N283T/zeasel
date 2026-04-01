/// Histogram with fixed-width bins for score distribution fitting.
///
/// Supports accumulating values into bins, computing basic statistics
/// (mean, variance, stddev), and chi-squared goodness-of-fit testing.
/// Bins are dynamically resized when values fall outside the current range.
///
/// "Full" mode additionally stores all raw sample values, enabling
/// sorted data retrieval, rank queries, tail extraction by mass fraction,
/// and goodness-of-fit testing with intelligent re-binning.
const std = @import("std");
const Allocator = std.mem.Allocator;
const math = std.math;
const gamma = @import("gamma.zig");

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

    /// Raw data storage (full mode only).
    raw_data: ?std.ArrayList(f64),
    /// Whether this histogram is in full mode (stores raw values).
    is_full: bool,
    /// Whether raw_data is currently sorted.
    is_sorted: bool,
    /// Tail threshold set by setTailByMass.
    tail_threshold: ?f64,

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
            .raw_data = null,
            .is_full = false,
            .is_sorted = true,
            .tail_threshold = null,
        };
    }

    /// Create a histogram in "full" mode that stores both binned counts AND raw values.
    /// The bin_width and min_val define the initial bin layout; bins grow dynamically.
    pub fn initFull(allocator: Allocator, bin_width: f64, min_val: f64) !Histogram {
        // Start with a single bin; will grow as data arrives.
        const initial_bins: usize = 1;
        const counts = try allocator.alloc(u64, initial_bins);
        @memset(counts, 0);
        const edges = try allocator.alloc(f64, initial_bins + 1);
        edges[0] = min_val;
        edges[1] = min_val + bin_width;
        return Histogram{
            .counts = counts,
            .edges = edges,
            .n_bins = initial_bins,
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
            .raw_data = .empty,
            .is_full = true,
            .is_sorted = true,
            .tail_threshold = null,
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
    /// In full mode, also stores the raw value.
    pub fn add(self: *Histogram, x: f64) !void {
        if (self.is_done) return error.HistogramFinished;

        self.total += 1;
        self.sum += x;
        self.sum_sq += x * x;
        if (x < self.min_val) self.min_val = x;
        if (x > self.max_val) self.max_val = x;

        // Store raw value in full mode.
        if (self.raw_data) |*rd| {
            try rd.append(self.allocator, x);
            self.is_sorted = false;
        }

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
    pub fn setTail(self: *Histogram, phi_val: f64) f64 {
        const bin_i = self.scoreToBin(phi_val);
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

    /// Ensure raw_data is sorted. No-op if already sorted.
    fn ensureSorted(self: *Histogram) void {
        if (!self.is_sorted) {
            if (self.raw_data) |*rd| {
                std.mem.sort(f64, rd.items, {}, std.sort.asc(f64));
                self.is_sorted = true;
            }
        }
    }

    /// Get sorted copy of raw data. Only available in full mode.
    /// Caller owns the returned slice and must free it with allocator.
    pub fn getData(self: *Histogram, allocator: Allocator) ![]f64 {
        const rd = self.raw_data orelse return error.NotFullMode;
        if (rd.items.len == 0) return error.NoData;
        const result = try allocator.alloc(f64, rd.items.len);
        @memcpy(result, rd.items);
        std.mem.sort(f64, result, {}, std.sort.asc(f64));
        return result;
    }

    /// Get the n-th highest value (0-indexed, 0 = highest).
    /// Only available in full mode.
    pub fn getRank(self: *Histogram, rank: usize) !f64 {
        const rd = self.raw_data orelse return error.NotFullMode;
        if (rd.items.len == 0) return error.NoData;
        if (rank >= rd.items.len) return error.RankOutOfRange;
        self.ensureSorted();
        return rd.items[rd.items.len - 1 - rank];
    }

    /// Set censoring threshold by mass fraction.
    /// tailp = fraction of data to include in the tail (e.g., 0.01 = top 1%).
    /// Returns the threshold value and the number of samples in the tail.
    /// Only available in full mode.
    pub fn setTailByMass(self: *Histogram, tailp: f64) !struct { threshold: f64, n_tail: usize } {
        const rd = self.raw_data orelse return error.NotFullMode;
        if (rd.items.len == 0) return error.NoData;
        if (tailp <= 0.0 or tailp > 1.0) return error.InvalidParameter;

        self.ensureSorted();

        const n = rd.items.len;
        // Number of samples in the tail: ceil(n * tailp), at least 1.
        const n_tail_f = @ceil(@as(f64, @floatFromInt(n)) * tailp);
        const n_tail: usize = @intFromFloat(@min(n_tail_f, @as(f64, @floatFromInt(n))));

        // Threshold is the value at the boundary.
        const threshold = rd.items[n - n_tail];

        self.tail_threshold = threshold;
        // n_censored here refers to the mass-based censoring count.
        self.n_censored = @intCast(n - n_tail);

        return .{ .threshold = threshold, .n_tail = n_tail };
    }

    /// Get raw tail values above the censoring threshold for a given mass fraction.
    /// Returns a sorted owned copy of tail values. Caller must free with allocator.
    /// Only available in full mode.
    pub fn getTailByMass(self: *Histogram, allocator: Allocator, tailp: f64) ![]f64 {
        const rd = self.raw_data orelse return error.NotFullMode;
        if (rd.items.len == 0) return error.NoData;
        if (tailp <= 0.0 or tailp > 1.0) return error.InvalidParameter;

        self.ensureSorted();

        const n = rd.items.len;
        const n_tail_f = @ceil(@as(f64, @floatFromInt(n)) * tailp);
        const n_tail: usize = @intFromFloat(@min(n_tail_f, @as(f64, @floatFromInt(n))));

        const start = n - n_tail;
        const result = try allocator.alloc(f64, n_tail);
        @memcpy(result, rd.items[start..n]);
        return result;
    }

    /// Result of goodness-of-fit testing.
    pub const GoodnessResult = struct {
        /// G-test statistic (log-likelihood ratio).
        g_stat: f64,
        /// Chi-squared statistic.
        chi2_stat: f64,
        /// P-value for G-test (from chi-squared distribution).
        g_pvalue: f64,
        /// P-value for chi-squared test.
        chi2_pvalue: f64,
        /// Degrees of freedom (number of merged bins - 1).
        ndof: usize,
    };

    /// Goodness-of-fit test with intelligent re-binning.
    /// Merges adjacent bins until each has at least `min_expected` expected count,
    /// then computes both G-test and chi-squared statistics with p-values.
    ///
    /// `expected_cdf` is the CDF of the hypothesized distribution.
    /// The test uses all bins (or the tail region if setTail was called).
    pub fn goodness(
        self: Histogram,
        allocator: Allocator,
        expected_cdf: *const fn (f64) f64,
        min_expected: f64,
    ) !GoodnessResult {
        // Determine the range of bins to test.
        const start_bin: usize = if (self.is_done) self.cmin else 0;
        const end_bin: usize = self.n_bins;
        const n_test_bins = end_bin - start_bin;
        if (n_test_bins == 0) return error.NoData;

        // Total observed in the test range.
        const n_total: f64 = blk: {
            if (self.is_done) {
                break :blk @floatFromInt(self.n_observed);
            } else {
                break :blk @floatFromInt(self.total);
            }
        };
        if (n_total == 0) return error.NoData;

        // Compute expected counts per bin: E[b] = N * (CDF(upper) - CDF(lower)).
        const expected = try allocator.alloc(f64, n_test_bins);
        defer allocator.free(expected);
        const observed = try allocator.alloc(f64, n_test_bins);
        defer allocator.free(observed);

        for (0..n_test_bins) |i| {
            const lo = self.edges[start_bin + i];
            const hi = self.edges[start_bin + i + 1];
            const cdf_lo = expected_cdf(lo);
            const cdf_hi = expected_cdf(hi);
            expected[i] = n_total * (cdf_hi - cdf_lo);
            observed[i] = @floatFromInt(self.counts[start_bin + i]);
        }

        // Merge adjacent bins until each merged bin has >= min_expected expected count.
        // We accumulate into a dynamic list of merged bins.
        var merged_obs: std.ArrayList(f64) = .empty;
        defer merged_obs.deinit(allocator);
        var merged_exp: std.ArrayList(f64) = .empty;
        defer merged_exp.deinit(allocator);

        var acc_obs: f64 = 0;
        var acc_exp: f64 = 0;
        for (0..n_test_bins) |i| {
            acc_obs += observed[i];
            acc_exp += expected[i];
            if (acc_exp >= min_expected) {
                try merged_obs.append(allocator, acc_obs);
                try merged_exp.append(allocator, acc_exp);
                acc_obs = 0;
                acc_exp = 0;
            }
        }
        // Merge leftover into the last bin.
        if (acc_exp > 0) {
            if (merged_exp.items.len > 0) {
                merged_obs.items[merged_obs.items.len - 1] += acc_obs;
                merged_exp.items[merged_exp.items.len - 1] += acc_exp;
            } else {
                // All bins merged into one — not enough data for meaningful test.
                try merged_obs.append(allocator, acc_obs);
                try merged_exp.append(allocator, acc_exp);
            }
        }

        const n_merged = merged_obs.items.len;
        if (n_merged < 2) return error.InsufficientBins;

        // Compute G-test and chi-squared statistics.
        var g_stat: f64 = 0;
        var chi2_stat: f64 = 0;
        for (0..n_merged) |i| {
            const o = merged_obs.items[i];
            const e = merged_exp.items[i];
            if (e > 0) {
                if (o > 0) {
                    g_stat += o * @log(o / e);
                }
                const diff = o - e;
                chi2_stat += (diff * diff) / e;
            }
        }
        g_stat *= 2.0;

        // Degrees of freedom = number of merged bins - 1.
        const ndof = n_merged - 1;

        // P-values from chi-squared survival function.
        // Chi-squared(k) = Gamma(k/2, 1/2) with mu=0.
        const k_half = @as(f64, @floatFromInt(ndof)) / 2.0;
        const g_pvalue = gamma.surv(g_stat, 0.0, 0.5, k_half);
        const chi2_pvalue = gamma.surv(chi2_stat, 0.0, 0.5, k_half);

        return GoodnessResult{
            .g_stat = g_stat,
            .chi2_stat = chi2_stat,
            .g_pvalue = g_pvalue,
            .chi2_pvalue = chi2_pvalue,
            .ndof = ndof,
        };
    }

    pub fn deinit(self: *Histogram) void {
        if (self.raw_data) |*rd| {
            rd.deinit(self.allocator);
        }
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

test "initFull: stores raw data and bins correctly" {
    const allocator = std.testing.allocator;
    var h = try Histogram.initFull(allocator, 1.0, 0.0);
    defer h.deinit();

    try h.add(1.5);
    try h.add(2.5);
    try h.add(0.5);

    try std.testing.expectEqual(@as(u64, 3), h.total);
    try std.testing.expect(h.is_full);
    try std.testing.expect(h.raw_data != null);
    try std.testing.expectEqual(@as(usize, 3), h.raw_data.?.items.len);

    const data = try h.getData(allocator);
    defer allocator.free(data);

    // Should be sorted: 0.5, 1.5, 2.5.
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), data[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), data[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), data[2], 1e-10);
}

test "getData: returns error for non-full histogram" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    try h.add(5.0);
    try std.testing.expectError(error.NotFullMode, h.getData(allocator));
}

test "getRank: returns n-th highest value" {
    const allocator = std.testing.allocator;
    var h = try Histogram.initFull(allocator, 1.0, 0.0);
    defer h.deinit();

    try h.add(3.0);
    try h.add(1.0);
    try h.add(4.0);
    try h.add(1.5);
    try h.add(9.0);

    // Sorted: 1.0, 1.5, 3.0, 4.0, 9.0
    // rank 0 = highest = 9.0
    // rank 1 = 4.0
    // rank 4 = lowest = 1.0
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), try h.getRank(0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), try h.getRank(1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), try h.getRank(2), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), try h.getRank(4), 1e-10);

    // Out of range.
    try std.testing.expectError(error.RankOutOfRange, h.getRank(5));
}

test "setTailByMass: extracts top 50% of 100 values" {
    const allocator = std.testing.allocator;
    var h = try Histogram.initFull(allocator, 1.0, 0.0);
    defer h.deinit();

    // Add 100 values: 1.0, 2.0, ..., 100.0
    for (1..101) |i| {
        try h.add(@as(f64, @floatFromInt(i)));
    }

    const result = try h.setTailByMass(0.5);
    // Top 50% = 50 values, threshold = value at index 50 = 51.0
    try std.testing.expectEqual(@as(usize, 50), result.n_tail);
    try std.testing.expectApproxEqAbs(@as(f64, 51.0), result.threshold, 1e-10);
    try std.testing.expect(h.tail_threshold != null);
}

test "getTailByMass: returns correct tail values" {
    const allocator = std.testing.allocator;
    var h = try Histogram.initFull(allocator, 1.0, 0.0);
    defer h.deinit();

    // Add 10 values: 1.0, 2.0, ..., 10.0
    for (1..11) |i| {
        try h.add(@as(f64, @floatFromInt(i)));
    }

    // Get top 30% = ceil(10 * 0.3) = 3 values: 8.0, 9.0, 10.0
    const tail = try h.getTailByMass(allocator, 0.3);
    defer allocator.free(tail);

    try std.testing.expectEqual(@as(usize, 3), tail.len);
    try std.testing.expectApproxEqAbs(@as(f64, 8.0), tail[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), tail[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), tail[2], 1e-10);
}

test "getTailByMass: full mass returns all values" {
    const allocator = std.testing.allocator;
    var h = try Histogram.initFull(allocator, 1.0, 0.0);
    defer h.deinit();

    try h.add(5.0);
    try h.add(3.0);
    try h.add(7.0);

    const tail = try h.getTailByMass(allocator, 1.0);
    defer allocator.free(tail);

    try std.testing.expectEqual(@as(usize, 3), tail.len);
}

test "setTailByMass: invalid parameters" {
    const allocator = std.testing.allocator;
    var h = try Histogram.initFull(allocator, 1.0, 0.0);
    defer h.deinit();

    try h.add(1.0);
    try std.testing.expectError(error.InvalidParameter, h.setTailByMass(0.0));
    try std.testing.expectError(error.InvalidParameter, h.setTailByMass(-0.5));
}

test "goodness: uniform distribution passes" {
    const allocator = std.testing.allocator;
    // Create a histogram with 10 bins over [0, 10).
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // Add 100 perfectly uniform values (10 per bin).
    for (0..10) |bin| {
        for (0..10) |_| {
            try h.add(@as(f64, @floatFromInt(bin)) + 0.5);
        }
    }

    // CDF for uniform distribution on [0, 10].
    const uniform_cdf = struct {
        fn f(x: f64) f64 {
            if (x <= 0) return 0;
            if (x >= 10) return 1;
            return x / 10.0;
        }
    }.f;

    const result = try h.goodness(allocator, &uniform_cdf, 5.0);

    // Perfect fit: both statistics should be ~0.
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.chi2_stat, 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.g_stat, 1e-10);
    // P-value should be 1.0 (or very close) for perfect fit.
    try std.testing.expect(result.g_pvalue > 0.99);
    try std.testing.expect(result.chi2_pvalue > 0.99);
    try std.testing.expectEqual(@as(usize, 9), result.ndof);
}

test "goodness: poor fit has low p-value" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // All 100 values in one bin.
    for (0..100) |_| {
        try h.add(0.5);
    }

    // Test against uniform CDF — should fail badly.
    const uniform_cdf = struct {
        fn f(x: f64) f64 {
            if (x <= 0) return 0;
            if (x >= 10) return 1;
            return x / 10.0;
        }
    }.f;

    const result = try h.goodness(allocator, &uniform_cdf, 5.0);

    // P-value should be very low for such a bad fit.
    try std.testing.expect(result.chi2_pvalue < 0.001);
    try std.testing.expect(result.g_pvalue < 0.001);
}

test "goodness: with tail censoring" {
    const allocator = std.testing.allocator;
    var h = try Histogram.init(allocator, 0, 10, 10);
    defer h.deinit();

    // 10 values per bin = 100 total.
    for (0..10) |bin| {
        for (0..10) |_| {
            try h.add(@as(f64, @floatFromInt(bin)) + 0.5);
        }
    }

    // Set tail at 5.0: only bins [5..9] are tested.
    // n_observed = 50 (bins 5-9).
    _ = h.setTail(5.0);

    // CDF for uniform on [5, 10] — correctly conditioned for the tail.
    // Expected = 50 * (CDF(hi) - CDF(lo)) = 50 * 0.2 = 10 per bin.
    const tail_uniform_cdf = struct {
        fn f(x: f64) f64 {
            if (x <= 5) return 0;
            if (x >= 10) return 1;
            return (x - 5.0) / 5.0;
        }
    }.f;

    const result = try h.goodness(allocator, &tail_uniform_cdf, 5.0);

    // Uniform data with matching conditioned CDF — perfect fit.
    try std.testing.expect(result.chi2_pvalue > 0.99);
    try std.testing.expectEqual(@as(usize, 4), result.ndof); // 5 bins - 1
}

test "initFull: bins grow dynamically" {
    const allocator = std.testing.allocator;
    var h = try Histogram.initFull(allocator, 1.0, 0.0);
    defer h.deinit();

    // Values outside initial range should trigger bin growth.
    try h.add(-5.0);
    try h.add(50.0);
    try h.add(0.5);

    try std.testing.expectEqual(@as(u64, 3), h.total);
    try std.testing.expect(h.edges[0] <= -5.0);
    try std.testing.expect(h.edges[h.n_bins] >= 50.0);

    // All values should be in bins.
    var bin_total: u64 = 0;
    for (h.counts) |c| bin_total += c;
    try std.testing.expectEqual(@as(u64, 3), bin_total);
}

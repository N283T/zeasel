// Mixture Dirichlet distributions for profile HMM priors.
//
// A mixture Dirichlet is a weighted mixture of K Dirichlet components:
//   P(p | alpha) = sum_k q_k * Dir(p | alpha_k)
// where q_k are mixture weights and alpha_k are Dirichlet parameters.
//
// Used in hmmbuild as the prior over amino acid/nucleotide emission
// distributions when estimating profile HMM parameters from small
// count vectors.

const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;
const stats_fn = @import("stats/functions.zig");

pub const MixDirichlet = struct {
    /// Number of mixture components.
    n_components: usize,
    /// Alphabet size (number of categories in each Dirichlet).
    n_categories: usize,
    /// Mixture weights q[0..n_components-1], sum to 1.
    weights: []f64,
    /// Dirichlet alpha parameters: alpha[k][i] for component k, category i.
    /// Stored as flat array: alpha[k * n_categories + i].
    alpha: []f64,
    allocator: Allocator,

    /// Create a mixture Dirichlet with given parameters.
    pub fn init(allocator: Allocator, n_components: usize, n_categories: usize) !MixDirichlet {
        const weights = try allocator.alloc(f64, n_components);
        @memset(weights, 1.0 / @as(f64, @floatFromInt(n_components)));

        const alpha_data = try allocator.alloc(f64, n_components * n_categories);
        @memset(alpha_data, 1.0); // uniform prior

        return .{
            .n_components = n_components,
            .n_categories = n_categories,
            .weights = weights,
            .alpha = alpha_data,
            .allocator = allocator,
        };
    }

    /// Get alpha[k][i].
    pub fn getAlpha(self: MixDirichlet, k: usize, i: usize) f64 {
        return self.alpha[k * self.n_categories + i];
    }

    /// Set alpha[k][i].
    pub fn setAlpha(self: *MixDirichlet, k: usize, i: usize, val: f64) void {
        self.alpha[k * self.n_categories + i] = val;
    }

    /// Log probability of a count vector c[0..n_categories-1] under component k.
    /// logP(c | alpha_k) = logGamma(sum(alpha_k)) - logGamma(sum(alpha_k + c))
    ///                    + sum_i [ logGamma(alpha_k[i] + c[i]) - logGamma(alpha_k[i]) ]
    pub fn logProbComponent(self: MixDirichlet, k: usize, counts: []const f64) f64 {
        const n = self.n_categories;
        var sum_alpha: f64 = 0;
        var sum_alpha_c: f64 = 0;
        var log_num: f64 = 0;
        var log_den: f64 = 0;

        for (0..n) |i| {
            const a = self.getAlpha(k, i);
            const c = counts[i];
            sum_alpha += a;
            sum_alpha_c += a + c;
            log_num += stats_fn.logGamma(a + c);
            log_den += stats_fn.logGamma(a);
        }

        return stats_fn.logGamma(sum_alpha) - stats_fn.logGamma(sum_alpha_c) + log_num - log_den;
    }

    /// Log probability of a count vector under the full mixture.
    /// logP(c) = log( sum_k q_k * P(c | alpha_k) )
    pub fn logProb(self: MixDirichlet, counts: []const f64) f64 {
        var max_log: f64 = -math.inf(f64);
        const lp = self.allocator.alloc(f64, self.n_components) catch return -math.inf(f64);
        defer self.allocator.free(lp);

        for (0..self.n_components) |k| {
            lp[k] = @log(self.weights[k]) + self.logProbComponent(k, counts);
            if (lp[k] > max_log) max_log = lp[k];
        }

        // Log-sum-exp for numerical stability
        var sum_exp: f64 = 0;
        for (lp[0..self.n_components]) |v| {
            sum_exp += @exp(v - max_log);
        }
        return max_log + @log(sum_exp);
    }

    /// Compute posterior component probabilities given observed counts.
    /// Returns P(k | counts) for each component. Caller owns the returned slice.
    pub fn posteriorWeights(self: MixDirichlet, counts: []const f64) ![]f64 {
        const post = try self.allocator.alloc(f64, self.n_components);
        errdefer self.allocator.free(post);

        var max_log: f64 = -math.inf(f64);
        for (0..self.n_components) |k| {
            post[k] = @log(self.weights[k]) + self.logProbComponent(k, counts);
            if (post[k] > max_log) max_log = post[k];
        }

        // Normalize
        var total: f64 = 0;
        for (post) |*v| {
            v.* = @exp(v.* - max_log);
            total += v.*;
        }
        for (post) |*v| v.* /= total;

        return post;
    }

    /// Compute the mean posterior probability vector given observed counts.
    /// This is the Dirichlet mixture estimate of the probabilities.
    /// result[i] = sum_k P(k|c) * (alpha_k[i] + c[i]) / (sum_j(alpha_k[j] + c[j]))
    pub fn meanPosterior(self: MixDirichlet, counts: []const f64) ![]f64 {
        const n = self.n_categories;
        const result = try self.allocator.alloc(f64, n);
        @memset(result, 0);
        errdefer self.allocator.free(result);

        const post = try self.posteriorWeights(counts);
        defer self.allocator.free(post);

        for (0..self.n_components) |k| {
            var sum_alpha_c: f64 = 0;
            for (0..n) |i| {
                sum_alpha_c += self.getAlpha(k, i) + counts[i];
            }
            for (0..n) |i| {
                result[i] += post[k] * (self.getAlpha(k, i) + counts[i]) / sum_alpha_c;
            }
        }

        return result;
    }

    pub fn deinit(self: *MixDirichlet) void {
        self.allocator.free(self.weights);
        self.allocator.free(self.alpha);
    }
};

// --- Tests ---

test "init: uniform mixture" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 3, 4);
    defer md.deinit();

    try std.testing.expectEqual(@as(usize, 3), md.n_components);
    try std.testing.expectEqual(@as(usize, 4), md.n_categories);
    // Weights should be uniform
    for (md.weights) |w| {
        try std.testing.expectApproxEqAbs(1.0 / 3.0, w, 1e-10);
    }
}

test "logProbComponent: uniform Dirichlet" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 1, 4);
    defer md.deinit();

    const counts = [_]f64{ 5, 3, 2, 1 };
    const lp = md.logProbComponent(0, &counts);
    // Should be finite
    try std.testing.expect(math.isFinite(lp));
    try std.testing.expect(lp < 0); // probability < 1
}

test "logProb: single component equals logProbComponent" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 1, 4);
    defer md.deinit();

    const counts = [_]f64{ 5, 3, 2, 1 };
    const lp_mix = md.logProb(&counts);
    const lp_comp = md.logProbComponent(0, &counts);

    try std.testing.expectApproxEqAbs(lp_comp, lp_mix, 1e-10);
}

test "posteriorWeights: single component is 1.0" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 1, 4);
    defer md.deinit();

    const counts = [_]f64{ 5, 3, 2, 1 };
    const post = try md.posteriorWeights(&counts);
    defer allocator.free(post);

    try std.testing.expectApproxEqAbs(@as(f64, 1.0), post[0], 1e-10);
}

test "meanPosterior: sums to 1" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 2, 4);
    defer md.deinit();

    const counts = [_]f64{ 10, 5, 3, 2 };
    const mean = try md.meanPosterior(&counts);
    defer allocator.free(mean);

    var total: f64 = 0;
    for (mean) |v| total += v;
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), total, 1e-10);
}

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

    /// Parse a mixture Dirichlet from text in the Easel file format.
    ///
    /// Format:
    /// ```
    /// # comment lines (optional)
    /// K Q          # K = n_categories, Q = n_components
    /// weight_1
    /// alpha_1_1 alpha_1_2 ... alpha_1_K
    /// weight_2
    /// alpha_2_1 alpha_2_2 ... alpha_2_K
    /// ...
    /// ```
    pub fn read(allocator: Allocator, data: []const u8) !MixDirichlet {
        var line_iter = std.mem.splitScalar(u8, data, '\n');
        var tokens: std.ArrayList(f64) = .empty;
        defer tokens.deinit(allocator);

        while (line_iter.next()) |line| {
            // Strip trailing CR for Windows line endings
            const trimmed = std.mem.trimRight(u8, line, " \t\r");
            // Skip empty lines and comment lines
            if (trimmed.len == 0) continue;
            if (trimmed[0] == '#') continue;

            var tok_iter = std.mem.tokenizeAny(u8, trimmed, " \t");
            while (tok_iter.next()) |tok| {
                const val = std.fmt.parseFloat(f64, tok) catch return error.InvalidFormat;
                try tokens.append(allocator, val);
            }
        }

        const items = tokens.items;
        if (items.len < 2) return error.InvalidFormat;

        const n_categories: usize = @intFromFloat(items[0]);
        const n_components: usize = @intFromFloat(items[1]);

        if (n_categories < 1 or n_components < 1) return error.InvalidFormat;

        // Each component needs 1 weight + n_categories alphas
        const expected = 2 + n_components * (1 + n_categories);
        if (items.len != expected) return error.InvalidFormat;

        var md = try MixDirichlet.init(allocator, n_components, n_categories);
        errdefer md.deinit();

        var idx: usize = 2;
        var weight_sum: f64 = 0;
        for (0..n_components) |k| {
            const w = items[idx];
            if (w < 0.0 or w > 1.0) return error.InvalidFormat;
            md.weights[k] = w;
            weight_sum += w;
            idx += 1;

            for (0..n_categories) |i| {
                const a = items[idx];
                if (a <= 0.0) return error.InvalidFormat;
                md.setAlpha(k, i, a);
                idx += 1;
            }
        }

        // Normalize weights (as Easel does)
        if (weight_sum <= 0.0) return error.InvalidFormat;
        for (md.weights[0..n_components]) |*w| w.* /= weight_sum;

        return md;
    }

    /// Write the mixture Dirichlet in the Easel text format.
    pub fn write(self: MixDirichlet, dest: std.io.AnyWriter) !void {
        try dest.print("{d} {d}\n", .{ self.n_categories, self.n_components });
        for (0..self.n_components) |k| {
            try dest.print("{d:.4}\n", .{self.weights[k]});
            for (0..self.n_categories) |i| {
                if (i > 0) try dest.writeAll(" ");
                try dest.print("{d:.4}", .{self.getAlpha(k, i)});
            }
            try dest.writeAll("\n");
        }
    }

    /// Validate that the mixture Dirichlet is well-formed:
    /// - weights sum to 1 (within tolerance)
    /// - all alpha parameters are positive
    pub fn validate(self: MixDirichlet, tol: f64) bool {
        var weight_sum: f64 = 0;
        for (self.weights[0..self.n_components]) |w| {
            if (w < 0.0) return false;
            weight_sum += w;
        }
        if (@abs(weight_sum - 1.0) > tol) return false;

        for (self.alpha[0..self.n_components * self.n_categories]) |a| {
            if (a <= 0.0) return false;
        }
        return true;
    }

    /// Compare two mixture Dirichlet instances within tolerance.
    /// Returns true if they have the same dimensions and all weights/alphas
    /// match within the given tolerance.
    pub fn compare(self: MixDirichlet, other: MixDirichlet, tol: f64) bool {
        if (self.n_components != other.n_components) return false;
        if (self.n_categories != other.n_categories) return false;

        for (0..self.n_components) |k| {
            if (@abs(self.weights[k] - other.weights[k]) > tol) return false;
        }
        for (0..self.n_components * self.n_categories) |idx| {
            if (@abs(self.alpha[idx] - other.alpha[idx]) > tol) return false;
        }
        return true;
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

test "read: parse simple mixture Dirichlet" {
    const allocator = std.testing.allocator;
    const data =
        \\# Test mixture Dirichlet
        \\4 2
        \\0.6
        \\1.0 2.0 3.0 4.0
        \\0.4
        \\4.0 3.0 2.0 1.0
    ;

    var md = try MixDirichlet.read(allocator, data);
    defer md.deinit();

    try std.testing.expectEqual(@as(usize, 2), md.n_components);
    try std.testing.expectEqual(@as(usize, 4), md.n_categories);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), md.weights[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), md.weights[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), md.getAlpha(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), md.getAlpha(0, 3), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), md.getAlpha(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), md.getAlpha(1, 3), 1e-10);
}

test "read: reject invalid format" {
    const allocator = std.testing.allocator;

    // Too few tokens
    try std.testing.expectError(error.InvalidFormat, MixDirichlet.read(allocator, "4"));
    // Negative alpha
    try std.testing.expectError(error.InvalidFormat, MixDirichlet.read(allocator, "2 1\n1.0\n-1.0 2.0"));
    // Wrong number of tokens
    try std.testing.expectError(error.InvalidFormat, MixDirichlet.read(allocator, "2 1\n1.0\n1.0"));
}

test "write: output format" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 2, 3);
    defer md.deinit();

    md.weights[0] = 0.7;
    md.weights[1] = 0.3;
    md.setAlpha(0, 0, 1.5);
    md.setAlpha(0, 1, 2.5);
    md.setAlpha(0, 2, 3.5);
    md.setAlpha(1, 0, 4.0);
    md.setAlpha(1, 1, 5.0);
    md.setAlpha(1, 2, 6.0);

    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);

    try md.write(buf.writer(allocator).any());

    const output = buf.items;
    // Should start with dimensions
    try std.testing.expect(std.mem.startsWith(u8, output, "3 2\n"));
}

test "write then read: round-trip" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 3, 4);
    defer md.deinit();

    md.weights[0] = 0.5;
    md.weights[1] = 0.3;
    md.weights[2] = 0.2;
    md.setAlpha(0, 0, 1.1);
    md.setAlpha(0, 1, 2.2);
    md.setAlpha(0, 2, 3.3);
    md.setAlpha(0, 3, 4.4);
    md.setAlpha(1, 0, 5.5);
    md.setAlpha(1, 1, 6.6);
    md.setAlpha(1, 2, 7.7);
    md.setAlpha(1, 3, 8.8);
    md.setAlpha(2, 0, 0.1);
    md.setAlpha(2, 1, 0.2);
    md.setAlpha(2, 2, 0.3);
    md.setAlpha(2, 3, 0.4);

    // Write
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try md.write(buf.writer(allocator).any());

    // Read back
    var md2 = try MixDirichlet.read(allocator, buf.items);
    defer md2.deinit();

    try std.testing.expect(md.compare(md2, 1e-4));
}

test "validate: valid mixture" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 2, 3);
    defer md.deinit();

    md.weights[0] = 0.6;
    md.weights[1] = 0.4;
    try std.testing.expect(md.validate(1e-10));
}

test "validate: weights do not sum to 1" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 2, 3);
    defer md.deinit();

    md.weights[0] = 0.6;
    md.weights[1] = 0.6; // sum = 1.2
    try std.testing.expect(!md.validate(1e-10));
}

test "validate: zero alpha rejected" {
    const allocator = std.testing.allocator;
    var md = try MixDirichlet.init(allocator, 1, 3);
    defer md.deinit();

    md.setAlpha(0, 1, 0.0);
    try std.testing.expect(!md.validate(1e-10));
}

test "compare: identical instances" {
    const allocator = std.testing.allocator;
    var md1 = try MixDirichlet.init(allocator, 2, 3);
    defer md1.deinit();
    var md2 = try MixDirichlet.init(allocator, 2, 3);
    defer md2.deinit();

    try std.testing.expect(md1.compare(md2, 1e-10));
}

test "compare: different dimensions" {
    const allocator = std.testing.allocator;
    var md1 = try MixDirichlet.init(allocator, 2, 3);
    defer md1.deinit();
    var md2 = try MixDirichlet.init(allocator, 3, 3);
    defer md2.deinit();

    try std.testing.expect(!md1.compare(md2, 1e-10));
}

test "compare: different weights" {
    const allocator = std.testing.allocator;
    var md1 = try MixDirichlet.init(allocator, 2, 3);
    defer md1.deinit();
    var md2 = try MixDirichlet.init(allocator, 2, 3);
    defer md2.deinit();

    md1.weights[0] = 0.8;
    md1.weights[1] = 0.2;
    try std.testing.expect(!md1.compare(md2, 1e-10));
}

/// Dirichlet distribution and Dirichlet mixture priors.
///
/// Used in HMMER's Bayesian estimation of emission and transition
/// probabilities from observed counts.
const std = @import("std");
const Allocator = std.mem.Allocator;
const math = std.math;

// ---------------------------------------------------------------------------
// Dirichlet
// ---------------------------------------------------------------------------

/// A single Dirichlet component with K concentration parameters.
pub const Dirichlet = struct {
    /// Concentration parameters alpha[0..K-1].
    alpha: []const f64,

    /// Compute the posterior mean probability vector given observed counts.
    ///
    /// result[i] = (counts[i] + alpha[i]) / sum_j(counts[j] + alpha[j])
    ///
    /// Caller owns the returned slice and must free it with the same allocator.
    pub fn posteriorMean(self: Dirichlet, allocator: Allocator, counts: []const f64) ![]f64 {
        const k = self.alpha.len;
        std.debug.assert(counts.len == k);

        const result = try allocator.alloc(f64, k);
        errdefer allocator.free(result);

        var total: f64 = 0;
        for (0..k) |i| {
            result[i] = counts[i] + self.alpha[i];
            total += result[i];
        }
        if (total > 0) {
            for (result) |*r| r.* /= total;
        }
        return result;
    }

    /// Log probability density of probability vector p under this Dirichlet.
    ///
    /// log P(p | alpha) = log Gamma(sum(alpha))
    ///                  - sum_i(log Gamma(alpha[i]))
    ///                  + sum_i((alpha[i] - 1) * log(p[i]))
    ///
    /// p[i] == 0 terms with alpha[i] > 1 contribute -inf; p[i] == 0 with
    /// alpha[i] == 1 contributes 0 (consistent with Laplace prior).
    pub fn logPdf(self: Dirichlet, p: []const f64) f64 {
        const k = self.alpha.len;
        std.debug.assert(p.len == k);

        var sum_alpha: f64 = 0;
        var sum_lgamma: f64 = 0;
        var sum_log: f64 = 0;

        for (0..k) |i| {
            sum_alpha += self.alpha[i];
            sum_lgamma += math.lgamma(f64, self.alpha[i]);
            if (p[i] > 0) {
                sum_log += (self.alpha[i] - 1.0) * @log(p[i]);
            }
        }
        return math.lgamma(f64, sum_alpha) - sum_lgamma + sum_log;
    }
};

// ---------------------------------------------------------------------------
// MixtureDirichlet
// ---------------------------------------------------------------------------

/// Mixture Dirichlet prior with Q components.
///
/// The prior is a convex combination of Q Dirichlet distributions:
///   P(p) = sum_q mixture_coeffs[q] * Dir(p | alphas[q])
pub const MixtureDirichlet = struct {
    /// Number of mixture components.
    q: usize,
    /// Prior mixture weights; must sum to 1.0. Length == q.
    mixture_coeffs: []const f64,
    /// Concentration parameters per component. Outer length == q, inner length == k.
    alphas: []const []const f64,
    /// Alphabet size K.
    k: usize,

    /// Compute the posterior mean probability vector given observed counts.
    ///
    /// Algorithm:
    ///   1. Compute posterior component weight for each component:
    ///      log w_q = log(mixture_coeffs[q]) + logDataLikelihood(counts, alphas[q])
    ///   2. Normalise using the log-sum-exp trick.
    ///   3. Posterior mean = sum_q w_q * component_posterior_mean_q
    ///
    /// Caller owns the returned slice and must free it with the same allocator.
    pub fn posteriorMean(self: MixtureDirichlet, allocator: Allocator, counts: []const f64) ![]f64 {
        std.debug.assert(counts.len == self.k);

        // Step 1: log posterior weight for each component.
        const log_weights = try allocator.alloc(f64, self.q);
        defer allocator.free(log_weights);

        for (0..self.q) |qi| {
            log_weights[qi] = @log(self.mixture_coeffs[qi]) +
                logDataLikelihood(counts, self.alphas[qi]);
        }

        // Step 2: normalise via log-sum-exp.
        var max_lw = log_weights[0];
        for (log_weights[1..]) |w| {
            if (w > max_lw) max_lw = w;
        }

        const weights = try allocator.alloc(f64, self.q);
        defer allocator.free(weights);

        var total_w: f64 = 0;
        for (0..self.q) |qi| {
            weights[qi] = @exp(log_weights[qi] - max_lw);
            total_w += weights[qi];
        }
        for (weights) |*w| w.* /= total_w;

        // Step 3: weighted average of component posterior means.
        const result = try allocator.alloc(f64, self.k);
        errdefer allocator.free(result);
        @memset(result, 0.0);

        for (0..self.q) |qi| {
            var comp_total: f64 = 0;
            for (0..self.k) |i| {
                comp_total += counts[i] + self.alphas[qi][i];
            }
            if (comp_total > 0) {
                for (0..self.k) |i| {
                    result[i] += weights[qi] * (counts[i] + self.alphas[qi][i]) / comp_total;
                }
            }
        }

        return result;
    }
};

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Log marginal likelihood of observed counts under a single Dirichlet:
///
///   log P(c | alpha) = log Gamma(sum(alpha)) - log Gamma(sum(alpha + c))
///                    + sum_i(log Gamma(alpha[i] + c[i]) - log Gamma(alpha[i]))
fn logDataLikelihood(counts: []const f64, alpha: []const f64) f64 {
    const k = counts.len;
    std.debug.assert(alpha.len == k);

    var sum_alpha: f64 = 0;
    var sum_alpha_c: f64 = 0;
    var sum_lgamma_diff: f64 = 0;

    for (0..k) |i| {
        sum_alpha += alpha[i];
        sum_alpha_c += alpha[i] + counts[i];
        sum_lgamma_diff += math.lgamma(f64, alpha[i] + counts[i]) - math.lgamma(f64, alpha[i]);
    }

    return math.lgamma(f64, sum_alpha) - math.lgamma(f64, sum_alpha_c) + sum_lgamma_diff;
}

// ---------------------------------------------------------------------------
// Built-in default priors
// ---------------------------------------------------------------------------

const laplace_amino_alpha = [20]f64{
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
};

/// Laplace (uniform) Dirichlet prior for 20 amino acids.
/// All alpha = 1.0; acts as additive smoothing / pseudocount of 1.
pub const default_amino_prior = Dirichlet{
    .alpha = &laplace_amino_alpha,
};

const laplace_dna_alpha = [4]f64{ 1.0, 1.0, 1.0, 1.0 };

/// Laplace (uniform) Dirichlet prior for the 4 DNA nucleotides.
pub const default_dna_prior = Dirichlet{
    .alpha = &laplace_dna_alpha,
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "Dirichlet posteriorMean: uniform alpha, skewed counts" {
    // alpha = [1,1,1,1], counts = [10,0,0,0]
    // result[i] = (counts[i] + 1) / (10+1+1+1+1) = [11/14, 1/14, 1/14, 1/14]
    const allocator = std.testing.allocator;
    const d = Dirichlet{ .alpha = &[_]f64{ 1, 1, 1, 1 } };
    const counts = [_]f64{ 10, 0, 0, 0 };

    const result = try d.posteriorMean(allocator, &counts);
    defer allocator.free(result);

    try std.testing.expect(math.approxEqAbs(f64, result[0], 11.0 / 14.0, 1e-10));
    try std.testing.expect(math.approxEqAbs(f64, result[1], 1.0 / 14.0, 1e-10));
    try std.testing.expect(math.approxEqAbs(f64, result[2], 1.0 / 14.0, 1e-10));
    try std.testing.expect(math.approxEqAbs(f64, result[3], 1.0 / 14.0, 1e-10));
}

test "Dirichlet posteriorMean: zero counts yields alpha-normalised" {
    // With zero counts the posterior mean is alpha / sum(alpha).
    const allocator = std.testing.allocator;
    const alpha = [_]f64{ 2.0, 4.0, 1.0, 3.0 };
    const d = Dirichlet{ .alpha = &alpha };
    const counts = [_]f64{ 0, 0, 0, 0 };

    const result = try d.posteriorMean(allocator, &counts);
    defer allocator.free(result);

    const sum_alpha: f64 = 2 + 4 + 1 + 3; // 10
    try std.testing.expect(math.approxEqAbs(f64, result[0], 2.0 / sum_alpha, 1e-10));
    try std.testing.expect(math.approxEqAbs(f64, result[1], 4.0 / sum_alpha, 1e-10));
    try std.testing.expect(math.approxEqAbs(f64, result[2], 1.0 / sum_alpha, 1e-10));
    try std.testing.expect(math.approxEqAbs(f64, result[3], 3.0 / sum_alpha, 1e-10));
}

test "Dirichlet logPdf: symmetric Dirichlet at uniform p" {
    // For Dir(alpha=1, K=4) at p = [0.25, 0.25, 0.25, 0.25]:
    // log P = lgamma(4) - 4*lgamma(1) + sum(0 * log(0.25)) = lgamma(4) - 0 + 0 = log(6)
    const d = Dirichlet{ .alpha = &[_]f64{ 1, 1, 1, 1 } };
    const p = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    const result = d.logPdf(&p);
    const expected = math.lgamma(f64, 4.0); // lgamma(4) = log(3!) = log(6)
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-10));
}

test "logDataLikelihood: uniform alpha [1,1,1,1], counts [1,1,1,1]" {
    // log P(c|alpha) = lgamma(4) - lgamma(8)
    //               + 4*(lgamma(2) - lgamma(1))
    //             = log(6) - log(5040) + 4*0
    //             = log(6/5040)
    const counts = [_]f64{ 1, 1, 1, 1 };
    const alpha = [_]f64{ 1, 1, 1, 1 };
    const result = logDataLikelihood(&counts, &alpha);
    const expected = @log(6.0 / 5040.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-10));
}

test "MixtureDirichlet posteriorMean: two components, skewed counts select component" {
    // Component 0: favours position 0 (alpha [10, 1, 1, 1])
    // Component 1: favours position 3 (alpha [1, 1, 1, 10])
    // Equal prior weights.
    const allocator = std.testing.allocator;

    const alpha0 = [_]f64{ 10, 1, 1, 1 };
    const alpha1 = [_]f64{ 1, 1, 1, 10 };
    const alphas = [_][]const f64{ &alpha0, &alpha1 };
    const coeffs = [_]f64{ 0.5, 0.5 };

    const mix = MixtureDirichlet{
        .q = 2,
        .mixture_coeffs = &coeffs,
        .alphas = &alphas,
        .k = 4,
    };

    // Counts heavily favouring position 0 → result[0] should dominate.
    {
        const counts_a = [_]f64{ 100, 0, 0, 0 };
        const result = try mix.posteriorMean(allocator, &counts_a);
        defer allocator.free(result);
        // result[0] should be much larger than result[3]
        try std.testing.expect(result[0] > result[3]);
    }

    // Counts heavily favouring position 3 → result[3] should dominate.
    {
        const counts_c = [_]f64{ 0, 0, 0, 100 };
        const result = try mix.posteriorMean(allocator, &counts_c);
        defer allocator.free(result);
        try std.testing.expect(result[3] > result[0]);
    }
}

test "MixtureDirichlet posteriorMean: probabilities sum to 1" {
    const allocator = std.testing.allocator;

    const alpha0 = [_]f64{ 2, 1, 1, 1 };
    const alpha1 = [_]f64{ 1, 2, 1, 1 };
    const alphas = [_][]const f64{ &alpha0, &alpha1 };
    const coeffs = [_]f64{ 0.6, 0.4 };

    const mix = MixtureDirichlet{
        .q = 2,
        .mixture_coeffs = &coeffs,
        .alphas = &alphas,
        .k = 4,
    };

    const counts = [_]f64{ 3, 7, 2, 1 };
    const result = try mix.posteriorMean(allocator, &counts);
    defer allocator.free(result);

    var total: f64 = 0;
    for (result) |v| total += v;
    try std.testing.expect(math.approxEqAbs(f64, total, 1.0, 1e-10));
}

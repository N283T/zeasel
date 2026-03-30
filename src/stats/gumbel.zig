/// Gumbel (Type I extreme value) distribution.
///
/// Parameters:
///   mu     - location parameter
///   lambda - scale parameter (must be > 0)
///
/// The standardized variable is z = lambda * (x - mu).
/// CDF: P(X <= x) = exp(-exp(-z))
const std = @import("std");
const math = std.math;

/// PDF: lambda * exp(-z - exp(-z)), where z = lambda*(x-mu)
pub fn pdf(x: f64, mu: f64, lambda: f64) f64 {
    const z = lambda * (x - mu);
    return lambda * @exp(-z - @exp(-z));
}

/// Log PDF: log(lambda) - z - exp(-z)
pub fn logPdf(x: f64, mu: f64, lambda: f64) f64 {
    const z = lambda * (x - mu);
    return @log(lambda) - z - @exp(-z);
}

/// CDF: exp(-exp(-z))
pub fn cdf(x: f64, mu: f64, lambda: f64) f64 {
    const z = lambda * (x - mu);
    return @exp(-@exp(-z));
}

/// Survival function: 1 - CDF
pub fn surv(x: f64, mu: f64, lambda: f64) f64 {
    return 1.0 - cdf(x, mu, lambda);
}

/// Log survival: log(1 - exp(-exp(-z)))
/// Uses log1p for numerical stability when survival is near 1 (z is very negative).
pub fn logSurv(x: f64, mu: f64, lambda: f64) f64 {
    const z = lambda * (x - mu);
    // log(1 - exp(-exp(-z))) = log1p(-exp(-exp(-z)))
    return math.log1p(-@exp(-@exp(-z)));
}

/// Inverse CDF (quantile): mu - (1/lambda) * log(-log(p))
pub fn invcdf(p: f64, mu: f64, lambda: f64) f64 {
    return mu - (1.0 / lambda) * @log(-@log(p));
}

/// Inverse survival: mu - (1/lambda) * log(-log(1-p))
pub fn invsurv(p: f64, mu: f64, lambda: f64) f64 {
    return invcdf(1.0 - p, mu, lambda);
}

/// Maximum likelihood estimation of Gumbel parameters from complete data.
///
/// Uses Newton-Raphson iteration to find lambda, then derives mu.
/// Returns an error if the input has fewer than 2 samples or zero variance.
pub fn fitComplete(x: []const f64) !struct { mu: f64, lambda: f64 } {
    if (x.len < 2) return error.InsufficientData;

    const n: f64 = @floatFromInt(x.len);

    // Compute sample mean and variance.
    var sum: f64 = 0.0;
    for (x) |xi| sum += xi;
    const mean = sum / n;

    var var_sum: f64 = 0.0;
    for (x) |xi| {
        const d = xi - mean;
        var_sum += d * d;
    }
    const variance = var_sum / n;
    if (variance == 0.0) return error.ZeroVariance;
    const stddev = @sqrt(variance);

    // Initial lambda estimate from method of moments.
    // For Gumbel, variance = pi^2 / (6 * lambda^2)
    // => lambda = pi / (sqrt(6) * stddev)
    var lambda = math.pi / (@sqrt(6.0) * stddev);

    // Newton-Raphson iteration for lambda.
    // Log-likelihood: n*log(lambda) - lambda*sum(xi) - sum(exp(-lambda*xi)) * ...
    // The score equation for lambda reduces to:
    //   1/lambda - mean(xi) + sum(xi * exp(-lambda*xi)) / sum(exp(-lambda*xi)) = 0
    for (0..100) |_| {
        var s0: f64 = 0.0; // sum of exp(-lambda * xi)
        var s1: f64 = 0.0; // sum of xi * exp(-lambda * xi)
        var s2: f64 = 0.0; // sum of xi^2 * exp(-lambda * xi)
        for (x) |xi| {
            const e = @exp(-lambda * xi);
            s0 += e;
            s1 += xi * e;
            s2 += xi * xi * e;
        }

        // Gradient of log-likelihood w.r.t. lambda:
        // dL/dlambda = n/lambda - sum(xi) + sum(xi*exp(-lambda*xi))/sum(exp(-lambda*xi))
        // Wait — the full Gumbel log-likelihood includes the (x-mu) terms.
        // Easel parametrizes as: score = 1/lambda + mean(xi*exp(-lambda*xi))/mean(exp(-lambda*xi)) - mean(xi)
        // where we've already profiled out mu via mu = -(1/lambda)*log(mean(exp(-lambda*xi))).
        const w_mean_x = s1 / s0; // weighted mean of x
        const gradient = (1.0 / lambda) + w_mean_x - mean;

        // Second derivative (negative Hessian diagonal):
        // d2L/dlambda2 = -n/lambda^2 - [ sum(x^2*e)/sum(e) - (sum(x*e)/sum(e))^2 ]
        const w_mean_x2 = s2 / s0;
        const hessian = -(1.0 / (lambda * lambda)) - (w_mean_x2 - w_mean_x * w_mean_x);

        const delta = gradient / hessian;
        lambda -= delta;
        if (lambda <= 0.0) lambda = 1e-10; // guard

        if (@abs(delta) < 1e-8 * lambda) break;
    }

    // Derive mu from the MLE relationship:
    // mu = -(1/lambda) * log( (1/n) * sum(exp(-lambda*xi)) )
    var s0: f64 = 0.0;
    for (x) |xi| s0 += @exp(-lambda * xi);
    const mu = -(1.0 / lambda) * @log(s0 / n);

    return .{ .mu = mu, .lambda = lambda };
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "gumbel pdf at x=0, mu=0, lambda=1" {
    // pdf = 1 * exp(-0 - exp(0)) = exp(-1) ≈ 0.36787944
    const result = pdf(0.0, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, @exp(-1.0), 1e-6));
}

test "gumbel cdf at x=0, mu=0, lambda=1" {
    // cdf = exp(-exp(0)) = exp(-1) ≈ 0.36787944
    const result = cdf(0.0, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, @exp(-1.0), 1e-6));
}

test "gumbel surv at x=0, mu=0, lambda=1" {
    // surv = 1 - exp(-1) ≈ 0.63212056
    const result = surv(0.0, 0.0, 1.0);
    const expected = 1.0 - @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "gumbel cdf + surv = 1" {
    const x: f64 = 2.5;
    const c = cdf(x, 1.0, 0.5);
    const s = surv(x, 1.0, 0.5);
    try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-14));
}

test "gumbel invcdf round-trip" {
    const mu: f64 = 3.0;
    const lambda: f64 = 0.7;
    const xs = [_]f64{ -2.0, 0.0, 1.0, 5.0, 10.0 };
    for (xs) |x| {
        const p = cdf(x, mu, lambda);
        const x2 = invcdf(p, mu, lambda);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-6));
    }
}

test "gumbel invsurv round-trip" {
    const mu: f64 = 3.0;
    const lambda: f64 = 0.7;
    // Avoid extreme tails where surv approaches 0 or 1 (numerical issues near boundaries).
    const xs = [_]f64{ 1.0, 2.0, 3.0, 5.0 };
    for (xs) |x| {
        const s = surv(x, mu, lambda);
        const x2 = invsurv(s, mu, lambda);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-6));
    }
}

test "gumbel fitComplete recovers parameters" {
    // Generate 2000 samples from Gumbel(mu=10, lambda=0.5) using invcdf on LCG uniform.
    const true_mu: f64 = 10.0;
    const true_lambda: f64 = 0.5;
    const n = 2000;

    var samples: [n]f64 = undefined;
    // Simple LCG: a=1664525, c=1013904223, m=2^32
    var state: u64 = 12345;
    for (&samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_lambda);
    }

    const fit = try fitComplete(&samples);
    // Tolerance is generous because 2000 samples has sampling variance.
    try std.testing.expect(math.approxEqAbs(f64, fit.mu, true_mu, 0.3));
    try std.testing.expect(math.approxEqAbs(f64, fit.lambda, true_lambda, 0.05));
}

test "gumbel fitComplete error on small input" {
    const x = [_]f64{1.0};
    try std.testing.expectError(error.InsufficientData, fitComplete(&x));
}

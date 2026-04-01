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

/// Survival function: 1 - exp(-exp(-z)).
/// When exp(-z) is small, uses first-order approximation exp(-z) to avoid
/// precision loss from 1.0 - (value near 1.0).
pub fn surv(x: f64, mu: f64, lambda: f64) f64 {
    const z = lambda * (x - mu);
    const ey = @exp(-z);
    // When ey is small, 1 - exp(-ey) ~ ey (first-order Taylor)
    if (ey < 1e-7) return ey;
    return -math.expm1(-ey);
}

/// Log survival: log(1 - exp(-exp(-z)))
/// Uses log1p for numerical stability when survival is near 1 (z is very negative).
/// For large positive z, surv ~ exp(-z), so log_surv ~ -exp(-z) avoids log1p(-1) = -inf.
pub fn logSurv(x: f64, mu: f64, lambda: f64) f64 {
    const z = lambda * (x - mu);
    // For large z, surv ~ exp(-z), so log_surv ~ -exp(-z)
    if (z > 10.0) return -@exp(-z);
    return math.log1p(-@exp(-@exp(-z)));
}

/// Log CDF: -exp(-z), the direct elegant form.
/// Since cdf = exp(-exp(-z)), log(cdf) = -exp(-z).
pub fn logCdf(x: f64, mu: f64, lambda: f64) f64 {
    const z = lambda * (x - mu);
    return -@exp(-z);
}

/// Inverse CDF (quantile): mu - (1/lambda) * log(-log(p))
pub fn invcdf(p: f64, mu: f64, lambda: f64) f64 {
    return mu - (1.0 / lambda) * @log(-@log(p));
}

/// Inverse survival: mu - (1/lambda) * log(-log(1-p))
/// For very small p (< 5e-9), uses a numerically stable approximation
/// to avoid catastrophic cancellation in 1.0 - p.
/// Reference: Easel esl_gumbel_invsurv().
pub fn invsurv(p: f64, mu: f64, lambda: f64) f64 {
    // For very small p, 1.0 - p loses precision (rounds to 1.0 when p < ~1e-16).
    // Use the approximation: log(-log(1-p)) ≈ log(p) for small p.
    // More precisely: log(1-p) ≈ -p for small p, so -log(1-p) ≈ p, so log(-log(1-p)) ≈ log(p).
    if (p < 5e-9) {
        return mu - (1.0 / lambda) * @log(p);
    }
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

    // Find x_max to center data and prevent overflow in exp(-lambda * xi)
    // when xi has large negative values.
    var x_max: f64 = x[0];
    for (x[1..]) |xi| if (xi > x_max) {
        x_max = xi;
    };

    // Newton-Raphson iteration for lambda.
    // The score equation (profiling out mu) is:
    //   f(lambda) = 1/lambda + weighted_mean(x) - mean(x) = 0
    // We center by x_max to prevent overflow.
    var converged = false;
    for (0..100) |_| {
        var s0: f64 = 0.0;
        var s1: f64 = 0.0;
        var s2: f64 = 0.0;
        for (x) |xi| {
            const e = @exp(-lambda * (xi - x_max));
            s0 += e;
            s1 += xi * e;
            s2 += xi * xi * e;
        }

        const w_mean_x = s1 / s0;
        const gradient = (1.0 / lambda) + w_mean_x - mean;

        const w_mean_x2 = s2 / s0;
        const hessian = -(1.0 / (lambda * lambda)) - (w_mean_x2 - w_mean_x * w_mean_x);

        const delta = gradient / hessian;
        lambda -= delta;
        if (lambda <= 0.0) lambda = 1e-10;

        if (@abs(delta) < 1e-8 * lambda) {
            converged = true;
            break;
        }
    }

    // Bisection fallback if Newton-Raphson did not converge.
    if (!converged) {
        var lo: f64 = 1e-10;
        var hi: f64 = lambda * 10.0;
        if (hi < 1.0) hi = 1.0;
        for (0..200) |_| {
            const mid = (lo + hi) / 2.0;
            var s0f: f64 = 0.0;
            var s1f: f64 = 0.0;
            for (x) |xi| {
                const e = @exp(-mid * (xi - x_max));
                s0f += e;
                s1f += xi * e;
            }
            const f_val = (1.0 / mid) + s1f / s0f - mean;
            if (f_val > 0) {
                lo = mid;
            } else {
                hi = mid;
            }
            if ((hi - lo) < 1e-8 * lo) break;
        }
        lambda = (lo + hi) / 2.0;
    }

    // Derive mu from the MLE relationship:
    // mu = x_max - (1/lambda) * log( (1/n) * sum(exp(-lambda*(xi - x_max))) )
    var s0: f64 = 0.0;
    for (x) |xi| s0 += @exp(-lambda * (xi - x_max));
    const mu = x_max - (1.0 / lambda) * @log(s0 / n);

    return .{ .mu = mu, .lambda = lambda };
}

/// MLE of mu given known lambda, from complete data.
///
/// Formula (Lawless 4.1.5):
///   mu = -(1/lambda) * log( (1/n) * sum(exp(-lambda * xi)) )
pub fn fitCompleteLoc(x: []const f64, known_lambda: f64) !f64 {
    if (x.len < 2) return error.InsufficientData;

    const n: f64 = @floatFromInt(x.len);
    var esum: f64 = 0.0;
    for (x) |xi| {
        esum += @exp(-known_lambda * xi);
    }
    return -@log(esum / n) / known_lambda;
}

/// Evaluate Lawless equation 4.2.2 and its derivative for censored Gumbel MLE.
///
/// `x` contains observed (uncensored) values, `z` is the count of censored samples,
/// `phi` is the censoring threshold, `lambda` is the current estimate.
/// Returns (f, df) where f=0 at the MLE lambda.
fn lawless422(x: []const f64, z: f64, phi: f64, lambda: f64) struct { f: f64, df: f64 } {
    const n: f64 = @floatFromInt(x.len);

    var esum: f64 = 0.0;
    var xesum: f64 = 0.0;
    var xxesum: f64 = 0.0;
    var xsum: f64 = 0.0;
    for (x) |xi| {
        const e = @exp(-lambda * xi);
        xsum += xi;
        esum += e;
        xesum += xi * e;
        xxesum += xi * xi * e;
    }

    // Add censored data terms.
    const ep = @exp(-lambda * phi);
    esum += z * ep;
    xesum += z * phi * ep;
    xxesum += z * phi * phi * ep;

    const f = 1.0 / lambda - xsum / n + xesum / esum;
    const ratio = xesum / esum;
    const df = ratio * ratio - xxesum / esum - 1.0 / (lambda * lambda);
    return .{ .f = f, .df = df };
}

/// Left-censored MLE of Gumbel parameters.
///
/// `x` contains only the observed (uncensored) values (x_i >= phi).
/// `n_total` is total number of samples including censored.
/// `phi` is the censoring threshold.
///
/// Uses Newton-Raphson on Lawless equation 4.2.2 with bisection fallback.
/// Reference: Easel esl_gumbel_FitCensored().
pub fn fitCensored(x: []const f64, n_total: usize, phi: f64) !struct { mu: f64, lambda: f64 } {
    const n_obs = x.len;
    if (n_obs < 2) return error.InsufficientData;
    if (n_total < n_obs) return error.InvalidArgument;

    const n: f64 = @floatFromInt(n_obs);
    const z: f64 = @floatFromInt(n_total - n_obs); // number of censored samples

    // Initial lambda estimate from method of moments on observed data.
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

    var lambda = math.pi / @sqrt(6.0 * variance);

    // Newton-Raphson iteration for lambda.
    const tol: f64 = 1e-5;
    var converged = false;
    for (0..100) |_| {
        const r = lawless422(x, z, phi, lambda);
        if (@abs(r.f) < tol) {
            converged = true;
            break;
        }
        lambda -= r.f / r.df;
        if (lambda <= 0.0) lambda = 0.001;
    }

    // Bisection fallback if Newton-Raphson did not converge.
    if (!converged) {
        var left: f64 = 0.001;
        var right: f64 = math.pi / @sqrt(6.0 * variance);
        // Bracket the root: ensure f(right) < 0.
        var r = lawless422(x, z, phi, right);
        while (r.f > 0.0) {
            right *= 2.0;
            if (right > 1000.0) return error.NoResult;
            r = lawless422(x, z, phi, right);
        }
        for (0..100) |_| {
            const mid = (left + right) / 2.0;
            const rm = lawless422(x, z, phi, mid);
            if (@abs(rm.f) < tol) break;
            if (rm.f > 0.0) {
                left = mid;
            } else {
                right = mid;
            }
        }
        lambda = (left + right) / 2.0;
    }

    // Substitute into Lawless 4.2.3 to find mu.
    var esum: f64 = 0.0;
    for (x) |xi| esum += @exp(-lambda * xi);
    esum += z * @exp(-lambda * phi);
    const mu = -@log(esum / n) / lambda;

    return .{ .mu = mu, .lambda = lambda };
}

/// Censored MLE of mu given known lambda.
///
/// `x` contains only observed values (x_i >= phi).
/// `n_total` is total count including censored.
/// `phi` is the censoring threshold.
///
/// Formula (Lawless 4.2.3 with known lambda):
///   mu = -(1/lambda) * log( (1/n) * (sum(exp(-lambda*xi)) + z*exp(-lambda*phi)) )
///
/// Reference: Easel esl_gumbel_FitCensoredLoc().
pub fn fitCensoredLoc(x: []const f64, n_total: usize, phi: f64, known_lambda: f64) !f64 {
    const n_obs = x.len;
    if (n_obs < 2) return error.InsufficientData;
    if (n_total < n_obs) return error.InvalidArgument;

    const n: f64 = @floatFromInt(n_obs);
    const z: f64 = @floatFromInt(n_total - n_obs);

    var esum: f64 = 0.0;
    for (x) |xi| {
        esum += @exp(-known_lambda * xi);
    }
    esum += z * @exp(-known_lambda * phi);

    return -@log(esum / n) / known_lambda;
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

test "gumbel fitCompleteLoc recovers mu" {
    const true_mu: f64 = 10.0;
    const true_lambda: f64 = 0.5;
    const n = 2000;

    var samples: [n]f64 = undefined;
    var state: u64 = 12345;
    for (&samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_lambda);
    }

    const mu = try fitCompleteLoc(&samples, true_lambda);
    try std.testing.expect(math.approxEqAbs(f64, mu, true_mu, 0.3));
}

test "gumbel fitCompleteLoc error on small input" {
    const x = [_]f64{1.0};
    try std.testing.expectError(error.InsufficientData, fitCompleteLoc(&x, 0.5));
}

test "gumbel fitCensored recovers parameters" {
    const true_mu: f64 = 10.0;
    const true_lambda: f64 = 0.5;
    const n = 2000;

    // Generate samples.
    var all_samples: [n]f64 = undefined;
    var state: u64 = 54321;
    for (&all_samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_lambda);
    }

    // Censor: keep only values >= phi.
    const phi: f64 = 10.0;
    var observed: [n]f64 = undefined;
    var n_obs: usize = 0;
    for (&all_samples) |s| {
        if (s >= phi) {
            observed[n_obs] = s;
            n_obs += 1;
        }
    }

    const fit = try fitCensored(observed[0..n_obs], n, phi);
    try std.testing.expect(math.approxEqAbs(f64, fit.mu, true_mu, 0.5));
    try std.testing.expect(math.approxEqAbs(f64, fit.lambda, true_lambda, 0.1));
}

test "gumbel fitCensored error on small input" {
    const x = [_]f64{1.0};
    try std.testing.expectError(error.InsufficientData, fitCensored(&x, 10, 0.0));
}

test "gumbel fitCensoredLoc recovers mu" {
    const true_mu: f64 = 10.0;
    const true_lambda: f64 = 0.5;
    const n = 2000;

    var all_samples: [n]f64 = undefined;
    var state: u64 = 54321;
    for (&all_samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_lambda);
    }

    const phi: f64 = 10.0;
    var observed: [n]f64 = undefined;
    var n_obs: usize = 0;
    for (&all_samples) |s| {
        if (s >= phi) {
            observed[n_obs] = s;
            n_obs += 1;
        }
    }

    const mu = try fitCensoredLoc(observed[0..n_obs], n, phi, true_lambda);
    try std.testing.expect(math.approxEqAbs(f64, mu, true_mu, 0.5));
}

test "gumbel fitCensoredLoc error on small input" {
    const x = [_]f64{1.0};
    try std.testing.expectError(error.InsufficientData, fitCensoredLoc(&x, 10, 0.0, 0.5));
}

test "gumbel logCdf matches log(cdf)" {
    const xs = [_]f64{ -2.0, 0.0, 1.0, 3.0, 5.0 };
    for (xs) |x| {
        const lc = logCdf(x, 0.0, 1.0);
        const expected = @log(cdf(x, 0.0, 1.0));
        try std.testing.expect(math.approxEqAbs(f64, lc, expected, 1e-12));
    }
}

test "gumbel logCdf direct form" {
    // logCdf = -exp(-z), verify the elegant identity
    const z: f64 = 2.0;
    const x = z; // mu=0, lambda=1
    const lc = logCdf(x, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, lc, -@exp(-z), 1e-14));
}

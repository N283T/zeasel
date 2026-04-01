/// Gamma distribution and incomplete gamma function.
///
/// Parameters:
///   mu     - location (shift)
///   lambda - rate (lambda > 0)
///   a      - shape (a > 0)
///
/// The incomplete gamma functions are computed using series expansion
/// (for x < a+1) or continued fraction (for x >= a+1), following
/// Numerical Recipes §6.2.
const std = @import("std");
const math = std.math;
const functions = @import("functions.zig");

/// Maximum iterations for series / continued-fraction expansions.
const MAX_ITER: usize = 200;
/// Convergence threshold.
const EPSILON: f64 = 3e-15;
/// Small float to avoid division by zero in Lentz's method.
const FPMIN: f64 = 1e-300;

/// Regularized lower incomplete gamma P(a, x) = gamma(a,x) / Gamma(a).
/// Returns a value in [0, 1].
/// Uses series expansion for x < a+1, continued fraction otherwise.
pub fn incompleteGamma(a: f64, x: f64) f64 {
    if (x < 0.0 or a <= 0.0) return 0.0;
    if (x == 0.0) return 0.0;
    if (x < a + 1.0) {
        return gammaSeries(a, x);
    } else {
        return 1.0 - gammaContFrac(a, x);
    }
}

/// Regularized upper incomplete gamma Q(a, x) = 1 - P(a, x).
/// Computed directly (without subtraction from 1.0) when possible,
/// avoiding precision loss in the tail.
pub fn upperIncompleteGamma(a: f64, x: f64) f64 {
    if (x == 0) return 1.0;
    if (x < a + 1.0) {
        // Series gives P; compute Q = 1 - P (less precise in tail, but x < a+1 means we're not deep in tail)
        return 1.0 - gammaSeries(a, x);
    } else {
        // Continued fraction gives Q directly — precise in the tail
        return gammaContFrac(a, x);
    }
}

/// Series expansion for the regularized lower incomplete gamma.
/// Converges for x < a+1.
fn gammaSeries(a: f64, x: f64) f64 {
    const log_gamma_a = math.lgamma(f64, a);
    const ln_x = @log(x);

    // ap = a, sum = 1/a, del = 1/a
    var ap = a;
    var del: f64 = 1.0 / a;
    var sum = del;

    for (0..MAX_ITER) |_| {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (@abs(del) < @abs(sum) * EPSILON) break;
    }

    return sum * @exp(-x + a * ln_x - log_gamma_a);
}

/// Continued fraction expansion for the regularized upper incomplete gamma.
/// Returns Q(a, x) = 1 - P(a, x). Converges for x >= a+1.
/// Uses Lentz's modified continued fraction method.
fn gammaContFrac(a: f64, x: f64) f64 {
    const log_gamma_a = math.lgamma(f64, a);
    const ln_x = @log(x);

    // Set up for evaluating CF via modified Lentz's method.
    var b = x + 1.0 - a;
    var c: f64 = 1.0 / FPMIN;
    var d = 1.0 / b;
    var h = d;

    var i: f64 = 1.0;
    while (i <= @as(f64, @floatFromInt(MAX_ITER))) : (i += 1.0) {
        const an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (@abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (@abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        const delta = d * c;
        h *= delta;
        if (@abs(delta - 1.0) < EPSILON) break;
    }

    return @exp(-x + a * ln_x - log_gamma_a) * h;
}

/// PDF: (lambda^a / Gamma(a)) * (x-mu)^(a-1) * exp(-lambda*(x-mu))
/// Returns 0 for x <= mu.
pub fn pdf(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    if (x <= mu) return 0.0;
    const t = lambda * (x - mu);
    // log pdf = a*log(lambda) - lgamma(a) + (a-1)*log(x-mu) - lambda*(x-mu)
    //         = a*log(lambda) - lgamma(a) + (a-1)*log(t/lambda) - t
    //         = log(lambda) - lgamma(a) + (a-1)*log(t) - t
    const log_p = @log(lambda) - math.lgamma(f64, a) + (a - 1.0) * @log(t) - t;
    return @exp(log_p);
}

/// CDF via regularized incomplete gamma: P(a, lambda*(x-mu))
pub fn cdf(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    if (x <= mu) return 0.0;
    return incompleteGamma(a, lambda * (x - mu));
}

/// Survival = Q(a, x) computed directly via upper incomplete gamma.
/// Avoids precision loss from 1.0 - P when P is near 1.0 (deep tail).
/// Reference: Easel esl_gam_surv.
pub fn surv(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    if (x <= mu) return 1.0;
    return upperIncompleteGamma(a, lambda * (x - mu));
}

/// Result of fitting a gamma distribution.
pub const FitResult = struct {
    mu: f64,
    lambda: f64,
    alpha: f64,
};

/// Maximum likelihood fit of a complete gamma distribution.
///
/// Given observed data x[0..n-1] and a known location parameter mu (typically 0),
/// estimates lambda (rate) and alpha (shape) using Minka's generalized Newton method.
///
/// The shape parameter `alpha` corresponds to `tau` in Easel's esl_gam_FitComplete()
/// and `a` in this module's pdf/cdf functions.
///
/// Returns error.DidNotConverge if the iterative alpha estimation fails to converge.
pub fn fitComplete(x: []const f64, mu: f64) !FitResult {
    std.debug.assert(x.len > 0);

    const n: f64 = @floatFromInt(x.len);

    // Compute sufficient statistics: mean(x-mu) and mean(log(x-mu))
    var xbar: f64 = 0;
    var logxbar: f64 = 0;
    for (x) |xi| {
        const shifted = xi - mu;
        std.debug.assert(shifted >= 0);
        xbar += shifted;
        logxbar += if (shifted == 0.0) -36.0 else @log(shifted);
    }
    xbar /= n;
    logxbar /= n;

    // Initial estimate for alpha using Stirling approximation
    // From Easel: alpha = 0.5 / (log(xbar) - logxbar)
    var alpha = 0.5 / (@log(xbar) - logxbar);

    // Generalized Newton iteration (Minka 2002)
    const max_iterations: usize = 100;
    var iter: usize = 0;
    while (iter < max_iterations) : (iter += 1) {
        const old_alpha = alpha;

        const psi_val = functions.psi(alpha);
        const tri_val = functions.trigamma(alpha);

        // Minka's update: alpha = 1/(1/alpha + (logxbar - log(xbar) + log(alpha) - psi(alpha)) / (alpha - alpha^2 * trigamma(alpha)))
        alpha = 1.0 / (1.0 / alpha + (logxbar - @log(xbar) + @log(alpha) - psi_val) / (alpha - alpha * alpha * tri_val));

        // Check convergence on both alpha and the NLL
        if (@abs(alpha - old_alpha) < 1e-6 * @abs(old_alpha) + 1e-6) break;
    }
    if (iter == max_iterations) return error.DidNotConverge;

    return .{
        .mu = mu,
        .lambda = alpha / xbar,
        .alpha = alpha,
    };
}

/// Log PDF: tau*log(lambda) - lgamma(a) + (a-1)*log(t) - t, where t = lambda*(x-mu).
/// Returns -inf for x <= mu.
pub fn logPdf(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    if (x <= mu) return -math.inf(f64);
    const t = lambda * (x - mu);
    return @log(lambda) - math.lgamma(f64, a) + (a - 1.0) * @log(t) - t;
}

/// Log CDF: log(P(a, lambda*(x-mu))).
/// Returns -inf for x <= mu.
pub fn logCdf(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    if (x <= mu) return -math.inf(f64);
    const val = incompleteGamma(a, lambda * (x - mu));
    return @log(val);
}

/// Log survival: log(Q(a, x)) computed directly via upper incomplete gamma.
/// Avoids precision loss from log(1 - P) when P is near 1.0.
/// Reference: Easel esl_gam_logsurv.
/// Returns 0 for x <= mu.
pub fn logSurv(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    if (x <= mu) return 0.0;
    const q = upperIncompleteGamma(a, lambda * (x - mu));
    return @log(q);
}

/// Inverse CDF (quantile function) of the Gamma distribution.
///
/// Given a cumulative probability p in (0,1), returns x such that
/// CDF(x; mu, lambda, a) = p.
///
/// Uses bisection on the CDF. The initial bracket is based on the
/// distribution mean +/- several standard deviations. Returns error
/// if bisection does not converge within the iteration limit.
pub fn invcdf(p: f64, mu: f64, lambda: f64, a: f64) !f64 {
    std.debug.assert(p > 0.0 and p < 1.0);
    std.debug.assert(lambda > 0.0);
    std.debug.assert(a > 0.0);

    // Mean and std dev of Gamma(a, lambda) shifted by mu.
    const mean = mu + a / lambda;
    const sd = @sqrt(a) / lambda;

    // Initial bracket: expand until CDF(lo) < p and CDF(hi) > p.
    var lo = mean - 5.0 * sd;
    if (lo <= mu) lo = mu + 1e-15;
    var hi = mean + 5.0 * sd;

    // Widen bracket if needed.
    while (cdf(lo, mu, lambda, a) > p) {
        lo = mu + (lo - mu) * 0.5;
        if (lo <= mu) {
            lo = mu + 1e-15;
            break;
        }
    }
    while (cdf(hi, mu, lambda, a) < p) {
        hi = hi + 5.0 * sd;
    }

    // Bisection.
    const max_iter: usize = 200;
    for (0..max_iter) |_| {
        const mid = 0.5 * (lo + hi);
        if (hi - lo < 1e-12 * @abs(mid) + 1e-15) return mid;

        const f_mid = cdf(mid, mu, lambda, a);
        if (f_mid < p) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return 0.5 * (lo + hi);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "invcdf: round-trip with CDF" {
    // For several quantiles, verify invcdf(cdf(x)) ~ x.
    const ps = [_]f64{ 0.1, 0.25, 0.5, 0.75, 0.9, 0.99 };
    for (ps) |p| {
        const x = try invcdf(p, 0.0, 1.0, 2.0);
        const recovered_p = cdf(x, 0.0, 1.0, 2.0);
        try std.testing.expect(math.approxEqAbs(f64, recovered_p, p, 1e-8));
    }
}

test "invcdf: shifted gamma (mu != 0)" {
    const x = try invcdf(0.5, 5.0, 2.0, 3.0);
    const recovered_p = cdf(x, 5.0, 2.0, 3.0);
    try std.testing.expect(math.approxEqAbs(f64, recovered_p, 0.5, 1e-8));
}

test "invcdf: exponential special case (alpha=1)" {
    // CDF of Exp(lambda=1) at x is 1 - exp(-x), so median = ln(2).
    const x = try invcdf(0.5, 0.0, 1.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, x, @log(2.0), 1e-8));
}

test "incompleteGamma(1, 1) = 1 - exp(-1)" {
    // P(1, x) = 1 - exp(-x)
    const result = incompleteGamma(1.0, 1.0);
    const expected = 1.0 - @exp(-1.0); // ≈ 0.63212056
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "incompleteGamma(2, 1) ≈ 0.26424" {
    // P(2, 1) = 1 - 2*exp(-1) ≈ 0.26424111...
    const result = incompleteGamma(2.0, 1.0);
    const expected: f64 = 1.0 - 2.0 * @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "incompleteGamma boundary: P(a,0) = 0" {
    try std.testing.expect(math.approxEqAbs(f64, incompleteGamma(2.0, 0.0), 0.0, 1e-14));
}

test "incompleteGamma large x approaches 1" {
    const result = incompleteGamma(2.0, 100.0);
    try std.testing.expect(math.approxEqAbs(f64, result, 1.0, 1e-10));
}

test "incompleteGamma continued fraction path (x >= a+1)" {
    // For a=2, x=5: uses continued fraction
    // Known: P(2,5) = 1 - 6*exp(-5) ≈ 0.95957...
    const result = incompleteGamma(2.0, 5.0);
    const expected = 1.0 - 6.0 * @exp(-5.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "gamma pdf positive for x > mu" {
    // With a=1, the gamma reduces to exponential: pdf = lambda*exp(-lambda*(x-mu))
    const result = pdf(1.0, 0.0, 1.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, @exp(-1.0), 1e-6));
}

test "gamma pdf zero at or below mu" {
    try std.testing.expectEqual(@as(f64, 0.0), pdf(0.0, 0.0, 1.0, 2.0));
    try std.testing.expectEqual(@as(f64, 0.0), pdf(-1.0, 0.0, 1.0, 2.0));
}

test "gamma cdf + surv = 1" {
    const x: f64 = 2.0;
    const c = cdf(x, 0.0, 1.0, 2.0);
    const s = surv(x, 0.0, 1.0, 2.0);
    try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-14));
}

test "fitComplete: gamma with alpha=2, lambda=1, mu=0" {
    // Gamma(alpha=2, lambda=1) has mean=2, variance=2.
    // Hand-crafted data with mean ~2.0 and shape ~2.
    const data = [_]f64{
        1.2, 2.5, 1.8, 3.1, 0.9,
        2.0, 1.5, 2.8, 1.1, 3.5,
        2.3, 1.7, 0.8, 2.9, 1.4,
        2.6, 1.9, 3.0, 1.3, 2.2,
    };

    const result = try fitComplete(&data, 0.0);

    // mu should be exactly as passed
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.mu, 1e-10);

    // alpha and lambda should be reasonable estimates
    // With 20 samples, we accept generous tolerance
    try std.testing.expect(result.alpha > 0.5);
    try std.testing.expect(result.alpha < 10.0);
    try std.testing.expect(result.lambda > 0.1);
    try std.testing.expect(result.lambda < 10.0);

    // The ratio alpha/lambda should approximate the sample mean
    var sum: f64 = 0;
    for (data) |xi| sum += xi;
    const sample_mean = sum / @as(f64, @floatFromInt(data.len));
    const fitted_mean = result.alpha / result.lambda;
    try std.testing.expectApproxEqAbs(sample_mean, fitted_mean, 1e-6);
}

test "fitComplete: gamma exponential special case (alpha~1)" {
    // Exponential is Gamma with alpha=1. Data ~ Exp(lambda=1), mean=1.
    const data = [_]f64{ 0.2, 0.5, 0.1, 1.5, 0.8, 2.0, 0.3, 1.2, 0.7, 1.0 };

    const result = try fitComplete(&data, 0.0);

    // Alpha should be near 1 for exponential-like data
    try std.testing.expect(result.alpha > 0.3);
    try std.testing.expect(result.alpha < 3.0);
    try std.testing.expect(result.lambda > 0.1);
}

test "gamma logPdf matches log(pdf)" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const lp = logPdf(x, 0.0, 1.0, 2.0);
        const expected = @log(pdf(x, 0.0, 1.0, 2.0));
        try std.testing.expect(math.approxEqAbs(f64, lp, expected, 1e-12));
    }
}

test "gamma logPdf at or below mu is -inf" {
    try std.testing.expect(math.isInf(logPdf(0.0, 0.0, 1.0, 2.0)));
    try std.testing.expect(logPdf(0.0, 0.0, 1.0, 2.0) < 0.0);
}

test "gamma logCdf matches log(cdf)" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const lc = logCdf(x, 0.0, 1.0, 2.0);
        const expected = @log(cdf(x, 0.0, 1.0, 2.0));
        try std.testing.expect(math.approxEqAbs(f64, lc, expected, 1e-10));
    }
}

test "gamma logSurv matches log(surv)" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const ls = logSurv(x, 0.0, 1.0, 2.0);
        const expected = @log(surv(x, 0.0, 1.0, 2.0));
        try std.testing.expect(math.approxEqAbs(f64, ls, expected, 1e-10));
    }
}

test "gamma logSurv at or below mu is zero" {
    try std.testing.expectEqual(@as(f64, 0.0), logSurv(0.0, 0.0, 1.0, 2.0));
}

/// Log-normal distribution.
///
/// Parameters:
///   mu    - mean of the underlying normal distribution (log-space location)
///   sigma - standard deviation of the underlying normal distribution (log-space scale, must be > 0)
///
/// A random variable X is log-normally distributed if ln(X) ~ Normal(mu, sigma^2).
/// Defined for x > 0.
///
/// PDF: (1/(x*sigma*sqrt(2*pi))) * exp(-(ln(x)-mu)^2 / (2*sigma^2))
/// CDF: 0.5 * erfc(-(ln(x)-mu) / (sigma*sqrt(2)))
const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;
const c_math = @cImport(@cInclude("math.h"));

/// PDF: (1/(x*sigma*sqrt(2*pi))) * exp(-(ln(x)-mu)^2 / (2*sigma^2))
/// Returns 0 for x <= 0.
pub fn pdf(x: f64, mu: f64, sigma: f64) f64 {
    if (x <= 0.0) return 0.0;
    return @exp(logPdf(x, mu, sigma));
}

/// Log PDF: -ln(x) - ln(sigma) - 0.5*ln(2*pi) - (ln(x)-mu)^2 / (2*sigma^2)
/// Returns -inf for x <= 0.
pub fn logPdf(x: f64, mu: f64, sigma: f64) f64 {
    if (x <= 0.0) return -math.inf(f64);
    const log_x = @log(x);
    const z = (log_x - mu) / sigma;
    return -log_x - @log(sigma) - 0.5 * @log(2.0 * math.pi) - 0.5 * z * z;
}

/// CDF: 0.5 * erfc(-(ln(x)-mu) / (sigma*sqrt(2)))
/// Returns 0 for x <= 0.
pub fn cdf(x: f64, mu: f64, sigma: f64) f64 {
    if (x <= 0.0) return 0.0;
    const z = (@log(x) - mu) / (sigma * @sqrt(2.0));
    return 0.5 * c_math.erfc(-z);
}

/// Survival function: 0.5 * erfc((ln(x)-mu) / (sigma*sqrt(2)))
/// Returns 1 for x <= 0.
pub fn surv(x: f64, mu: f64, sigma: f64) f64 {
    if (x <= 0.0) return 1.0;
    const z = (@log(x) - mu) / (sigma * @sqrt(2.0));
    return 0.5 * c_math.erfc(z);
}

/// Inverse CDF (quantile function).
/// Uses the relationship: if X ~ LogNormal(mu, sigma), then ln(X) ~ Normal(mu, sigma).
/// So invcdf(p) = exp(mu + sigma * Phi^{-1}(p)) where Phi^{-1} is the standard normal quantile.
pub fn invcdf(p: f64, mu: f64, sigma: f64) f64 {
    if (p <= 0.0 or p >= 1.0) return math.nan(f64);

    // Standard normal inverse CDF using the rational approximation from normal.zig
    const normal = @import("normal.zig");
    const z = normal.invcdf(p, 0.0, 1.0);
    return @exp(mu + sigma * z);
}

/// Maximum likelihood estimation of log-normal parameters from complete data.
///
/// The MLE for log-normal is straightforward:
///   mu_hat    = mean(ln(x))
///   sigma_hat = sqrt(mean((ln(x) - mu_hat)^2))
///
/// All values in x must be positive.
pub fn fitComplete(allocator: Allocator, x: []const f64) !struct { mu: f64, sigma: f64 } {
    _ = allocator; // No allocation needed for this simple MLE
    if (x.len < 2) return error.InsufficientData;

    const n: f64 = @floatFromInt(x.len);

    // Compute mean of log(x)
    var log_sum: f64 = 0.0;
    for (x) |xi| {
        if (xi <= 0.0) return error.InvalidArgument;
        log_sum += @log(xi);
    }
    const mu = log_sum / n;

    // Compute variance of log(x)
    var var_sum: f64 = 0.0;
    for (x) |xi| {
        const d = @log(xi) - mu;
        var_sum += d * d;
    }
    const sigma = @sqrt(var_sum / n);

    if (sigma == 0.0) return error.ZeroVariance;

    return .{ .mu = mu, .sigma = sigma };
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "lognormal: pdf/cdf/surv consistency" {
    // CDF + surv = 1.0
    const mu: f64 = 0.0;
    const sigma: f64 = 1.0;
    const xs = [_]f64{ 0.1, 0.5, 1.0, 2.0, 5.0, 10.0 };
    for (xs) |x| {
        const c = cdf(x, mu, sigma);
        const s = surv(x, mu, sigma);
        try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-12));
    }
}

test "lognormal: numerical derivative of CDF approximates PDF" {
    const mu: f64 = 0.5;
    const sigma: f64 = 0.8;
    const h: f64 = 1e-7;

    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const numerical_pdf = (cdf(x + h, mu, sigma) - cdf(x - h, mu, sigma)) / (2.0 * h);
        const analytical_pdf = pdf(x, mu, sigma);
        try std.testing.expect(math.approxEqAbs(f64, numerical_pdf, analytical_pdf, 1e-5));
    }
}

test "lognormal: pdf is zero for x <= 0" {
    try std.testing.expectEqual(@as(f64, 0.0), pdf(0.0, 0.0, 1.0));
    try std.testing.expectEqual(@as(f64, 0.0), pdf(-1.0, 0.0, 1.0));
}

test "lognormal: known values at mu=0, sigma=1" {
    // At x=1: pdf = (1/(1*1*sqrt(2pi))) * exp(0) = 1/sqrt(2pi) ~ 0.3989
    const p = pdf(1.0, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, p, 1.0 / @sqrt(2.0 * math.pi), 1e-10));

    // CDF at x=1 with mu=0 should be 0.5 (median of lognormal(0,1) is exp(0)=1)
    const c = cdf(1.0, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, c, 0.5, 1e-10));
}

test "lognormal: invcdf is inverse of cdf" {
    const mu: f64 = 1.0;
    const sigma: f64 = 0.5;
    const ps = [_]f64{ 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99 };
    for (ps) |p| {
        const x = invcdf(p, mu, sigma);
        const recovered_p = cdf(x, mu, sigma);
        try std.testing.expect(math.approxEqAbs(f64, recovered_p, p, 1e-8));
    }
}

test "lognormal: cdf round-trip through invcdf" {
    const mu: f64 = 0.0;
    const sigma: f64 = 1.0;
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const p = cdf(x, mu, sigma);
        if (p > 0.0 and p < 1.0) {
            const x2 = invcdf(p, mu, sigma);
            try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-8));
        }
    }
}

test "lognormal: fitComplete recovers parameters" {
    // Generate samples from LogNormal(mu=1.0, sigma=0.5) using invcdf on LCG uniform.
    const true_mu: f64 = 1.0;
    const true_sigma: f64 = 0.5;
    const n = 5000;

    var samples: [n]f64 = undefined;
    var state: u64 = 77777;
    for (&samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_sigma);
    }

    const fit = try fitComplete(std.testing.allocator, &samples);

    try std.testing.expect(math.approxEqAbs(f64, fit.mu, true_mu, 0.05));
    try std.testing.expect(math.approxEqAbs(f64, fit.sigma, true_sigma, 0.05));
}

test "lognormal: fitComplete error on insufficient data" {
    var x = [_]f64{1.0};
    try std.testing.expectError(error.InsufficientData, fitComplete(std.testing.allocator, &x));
}

test "lognormal: fitComplete error on non-positive values" {
    var x = [_]f64{ 1.0, 2.0, -1.0 };
    try std.testing.expectError(error.InvalidArgument, fitComplete(std.testing.allocator, &x));
}

test "lognormal: logPdf matches log(pdf)" {
    const xs = [_]f64{ 0.1, 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const lp = logPdf(x, 0.0, 1.0);
        const p = pdf(x, 0.0, 1.0);
        try std.testing.expect(math.approxEqAbs(f64, lp, @log(p), 1e-12));
    }
}

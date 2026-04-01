/// Exponential distribution with location (mu) and rate (lambda).
///
/// Defined for x >= mu.
/// PDF: lambda * exp(-lambda * (x - mu))
const std = @import("std");
const math = std.math;

/// PDF: lambda * exp(-lambda * (x - mu)) for x >= mu, else 0.
pub fn pdf(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return 0.0;
    return lambda * @exp(-lambda * (x - mu));
}

/// CDF: 1 - exp(-lambda * (x - mu)) for x >= mu, else 0.
/// Uses Taylor approximation (1 - exp(-y) ≈ y) for small y to avoid
/// catastrophic cancellation near mu.
pub fn cdf(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return 0.0;
    const y = lambda * (x - mu);
    if (y < 5e-9) return y; // Taylor: 1-exp(-y) ≈ y for small y
    return 1.0 - @exp(-y);
}

/// Survival function: exp(-lambda * (x - mu)) for x >= mu, else 1.
pub fn surv(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return 1.0;
    return @exp(-lambda * (x - mu));
}

/// Log PDF: log(lambda) - lambda*(x - mu) for x >= mu, else -inf.
/// Guards against lambda = infinity to prevent NaN.
pub fn logPdf(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return -math.inf(f64);
    if (math.isPositiveInf(lambda)) {
        return if (x == mu) math.inf(f64) else -math.inf(f64);
    }
    return @log(lambda) - lambda * (x - mu);
}

/// Log CDF: log(1 - exp(-lambda*(x - mu))) for x >= mu, else -inf.
/// Uses log1p and direct approximations for numerical stability.
pub fn logCdf(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return -math.inf(f64);
    const y = lambda * (x - mu);
    if (y == 0.0) return -math.inf(f64);
    const ey = @exp(-y);
    // When y is small, 1-e^{-y} ~ y, so log(cdf) ~ log(y).
    if (y < 1e-7) return @log(y);
    // When y is large, e^{-y} is tiny, so log(1-e^{-y}) ~ -e^{-y}.
    if (ey < 1e-7) return -ey;
    return math.log1p(-ey);
}

/// Log survival: -lambda * (x - mu) for x >= mu, else 0.
pub fn logSurv(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return 0.0;
    return -lambda * (x - mu);
}

/// Inverse CDF: mu - (1/lambda) * log(1 - p)
pub fn invcdf(p: f64, mu: f64, lambda: f64) f64 {
    return mu - (1.0 / lambda) * @log(1.0 - p);
}

/// Inverse survival: mu - (1/lambda) * log(p).
/// Reference: Easel esl_exp_invsurv().
pub fn invsurv(p: f64, mu: f64, lambda: f64) f64 {
    return mu - (1.0 / lambda) * @log(p);
}

/// Result of fitting an exponential distribution.
pub const FitResult = struct {
    mu: f64,
    lambda: f64,
};

/// Maximum likelihood fit of a complete exponential distribution.
///
/// Given observed data x[0..n-1], estimates mu (location) and lambda (rate).
/// mu is set to the minimum observed value (left truncation point),
/// and lambda = n / sum(x_i - mu), the analytic MLE.
///
/// Corresponds to Easel's esl_exp_FitComplete().
pub fn fitComplete(x: []const f64) FitResult {
    std.debug.assert(x.len > 0);

    // ML mu is the minimum observation
    var mu = x[0];
    for (x[1..]) |xi| {
        if (xi < mu) mu = xi;
    }

    var sum: f64 = 0;
    for (x) |xi| {
        sum += xi - mu;
    }
    const mean = sum / @as(f64, @floatFromInt(x.len));

    return .{
        .mu = mu,
        .lambda = 1.0 / mean,
    };
}

/// Maximum likelihood fit with a known (fixed) location parameter mu.
///
/// Corresponds to Easel's esl_exp_FitCompleteScale().
pub fn fitCompleteScale(x: []const f64, mu: f64) f64 {
    std.debug.assert(x.len > 0);

    var sum: f64 = 0;
    for (x) |xi| {
        sum += xi - mu;
    }
    const mean = sum / @as(f64, @floatFromInt(x.len));
    return 1.0 / mean;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "exponential pdf at x=1, mu=0, lambda=1" {
    // pdf = exp(-1) ≈ 0.36787944
    const result = pdf(1.0, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, @exp(-1.0), 1e-6));
}

test "exponential cdf at x=1, mu=0, lambda=1" {
    // cdf = 1 - exp(-1) ≈ 0.63212056
    const result = cdf(1.0, 0.0, 1.0);
    const expected = 1.0 - @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "exponential pdf below mu is zero" {
    try std.testing.expectEqual(@as(f64, 0.0), pdf(-1.0, 0.0, 1.0));
}

test "exponential cdf below mu is zero" {
    try std.testing.expectEqual(@as(f64, 0.0), cdf(-1.0, 0.0, 1.0));
}

test "exponential cdf small-value guard near mu" {
    // For very small y = lambda*(x-mu), cdf should use Taylor approx: 1-exp(-y) ≈ y
    const tiny = 1e-12;
    const result = cdf(tiny, 0.0, 1.0);
    // Should be approximately equal to y = tiny
    try std.testing.expectApproxEqAbs(tiny, result, 1e-20);
    // At mu exactly, should be 0
    try std.testing.expectEqual(@as(f64, 0.0), cdf(0.0, 0.0, 1.0));
}

test "exponential surv below mu is one" {
    try std.testing.expectEqual(@as(f64, 1.0), surv(-1.0, 0.0, 1.0));
}

test "exponential cdf + surv = 1" {
    const xs = [_]f64{ 0.0, 0.5, 1.0, 3.0 };
    for (xs) |x| {
        const c = cdf(x, 0.0, 1.0);
        const s = surv(x, 0.0, 1.0);
        try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-14));
    }
}

test "exponential invcdf round-trip" {
    const mu: f64 = 2.0;
    const lambda: f64 = 1.5;
    const xs = [_]f64{ 2.0, 2.5, 3.0, 5.0 };
    for (xs) |x| {
        const p = cdf(x, mu, lambda);
        const x2 = invcdf(p, mu, lambda);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-6));
    }
}

test "fitComplete: exponential samples with known parameters" {
    // Samples from Exp(mu=1.0, lambda=2.0): mean should be mu + 1/lambda = 1.5
    // Hand-crafted data shifted by mu=1.0, with mean(x-mu)=0.5 -> lambda=2.0
    const data = [_]f64{ 1.1, 1.05, 1.2, 1.8, 1.3, 1.15, 1.5, 1.25, 1.4, 1.35 };
    const result = fitComplete(&data);

    // mu should be the minimum value
    try std.testing.expectApproxEqAbs(@as(f64, 1.05), result.mu, 1e-10);

    // lambda = 1 / mean(x - mu)
    var sum: f64 = 0;
    for (data) |xi| sum += xi - 1.05;
    const expected_lambda = @as(f64, @floatFromInt(data.len)) / sum;
    try std.testing.expectApproxEqAbs(expected_lambda, result.lambda, 1e-10);
}

test "fitComplete: single sample" {
    const data = [_]f64{5.0};
    const result = fitComplete(&data);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result.mu, 1e-10);
    // mean(x - mu) = 0, lambda = inf
    try std.testing.expect(math.isInf(result.lambda));
}

test "fitCompleteScale: known mu" {
    // With mu=0, lambda should be 1/mean(x)
    const data = [_]f64{ 0.5, 1.0, 1.5, 2.0, 0.8, 1.2 };
    const lambda = fitCompleteScale(&data, 0.0);
    var sum: f64 = 0;
    for (data) |xi| sum += xi;
    const expected = @as(f64, @floatFromInt(data.len)) / sum;
    try std.testing.expectApproxEqAbs(expected, lambda, 1e-10);
}

test "exponential logPdf matches log(pdf)" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const lp = logPdf(x, 0.0, 1.0);
        const expected = @log(pdf(x, 0.0, 1.0));
        try std.testing.expect(math.approxEqAbs(f64, lp, expected, 1e-12));
    }
}

test "exponential logPdf below mu is -inf" {
    try std.testing.expect(math.isInf(logPdf(-1.0, 0.0, 1.0)));
    try std.testing.expect(logPdf(-1.0, 0.0, 1.0) < 0.0);
}

test "exponential logCdf matches log(cdf)" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const lc = logCdf(x, 0.0, 1.0);
        const expected = @log(cdf(x, 0.0, 1.0));
        try std.testing.expect(math.approxEqAbs(f64, lc, expected, 1e-10));
    }
}

test "exponential logSurv matches log(surv)" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const ls = logSurv(x, 0.0, 1.0);
        const expected = @log(surv(x, 0.0, 1.0));
        try std.testing.expect(math.approxEqAbs(f64, ls, expected, 1e-12));
    }
}

test "exponential logSurv below mu is zero" {
    try std.testing.expectEqual(@as(f64, 0.0), logSurv(-1.0, 0.0, 1.0));
}

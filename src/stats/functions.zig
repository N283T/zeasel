// Foundation statistical functions shared across distribution modules.
//
// Provides: logGamma, psi (digamma), trigamma, erfc, chiSquaredTest,
// gTest, linearRegression, mean/variance.

const std = @import("std");
const math = std.math;
const gamma_mod = @import("gamma.zig");
const c_math = @cImport(@cInclude("math.h"));

/// Compute log(Gamma(x)) for x > 0.
/// Uses Lanczos approximation with 11 coefficients.
pub fn logGamma(x: f64) f64 {
    std.debug.assert(x > 0);

    const cof = [_]f64{
        4.694580336184385e+04,
        -1.560605207784446e+05,
        2.065049568014106e+05,
        -1.388934775095388e+05,
        5.031796415085709e+04,
        -9.601592329182778e+03,
        8.785855930895250e+02,
        -3.155153906098611e+01,
        2.908143421162229e-01,
        -2.319827630494973e-04,
        1.251639670050933e-10,
    };

    const xx = x - 1.0;
    var tx = xx + 11.0;
    var tmp = tx;
    var value: f64 = 1.0;
    var i: usize = 10;
    while (true) {
        value += cof[i] / tmp;
        tmp -= 1.0;
        if (i == 0) break;
        i -= 1;
    }
    value = @log(value);
    tx += 0.5;
    value += 0.918938533 + (xx + 0.5) * @log(tx) - tx;
    return value;
}

/// Compute Psi(x) (digamma function), the derivative of log(Gamma(x)).
/// Implements Bernardo's "Algorithm AS 103".
pub fn psi(x_in: f64) f64 {
    std.debug.assert(x_in > 0);

    if (x_in <= 1e-5) return -0.5772156649015329 - 1.0 / x_in;

    var result: f64 = 0;
    var x = x_in;

    // Use Psi(1+x) = Psi(x) + 1/x to shift x up
    while (x < 8.5) {
        result -= 1.0 / x;
        x += 1.0;
    }

    // Stirling approximation for large x
    const x2 = 1.0 / x;
    result += @log(x) - 0.5 * x2;
    const x4 = x2 * x2;
    result += (-1.0 / 12.0) * x4 + (1.0 / 120.0) * x4 * x4 - (1.0 / 252.0) * x4 * x4 * x4;

    return result;
}

/// Compute Psi'(x) (trigamma function), the second derivative of log(Gamma(x)).
/// Implements Schneider's "Algorithm AS 121".
pub fn trigamma(x_in: f64) f64 {
    std.debug.assert(x_in > 0);

    if (x_in <= 1e-4) return 1.0 / (x_in * x_in);

    var result: f64 = 0;
    var x = x_in;

    while (x < 5.0) {
        result += 1.0 / (x * x);
        x += 1.0;
    }

    const y = 1.0 / (x * x);
    result += 0.5 * y + (1.0 + y * ((1.0 / 6.0) + y * ((1.0 / 30.0) + y * ((1.0 / 42.0) + y * (1.0 / 30.0))))) / x;

    return result;
}

/// Complementary error function erfc(x) = 1 - erf(x).
/// Uses C standard library for full double precision (~15 digits).
pub fn erfc(x: f64) f64 {
    return c_math.erfc(x);
}

/// Chi-squared test: P(chi^2 > x) with v degrees of freedom.
/// Returns the p-value (survival function of chi-squared distribution).
pub fn chiSquaredTest(v: u32, x: f64) f64 {
    // P(chi^2 > x) = Q(v/2, x/2) = 1 - P(v/2, x/2) where P is the
    // regularized incomplete gamma function.
    return 1.0 - gamma_mod.incompleteGamma(@as(f64, @floatFromInt(v)) / 2.0, x / 2.0);
}

/// G-test (log-likelihood ratio test) for comparing two proportions.
/// ca/na vs cb/nb. Returns the G statistic and p-value.
pub fn gTest(ca: u32, na: u32, cb: u32, nb: u32) struct { g: f64, p: f64 } {
    const a: f64 = @floatFromInt(ca);
    const b: f64 = @floatFromInt(na - ca);
    const c: f64 = @floatFromInt(cb);
    const d: f64 = @floatFromInt(nb - cb);
    const n: f64 = @as(f64, @floatFromInt(na)) + @as(f64, @floatFromInt(nb));

    var llr: f64 = 0;
    if (a > 0) llr += a * @log(a);
    if (b > 0) llr += b * @log(b);
    if (c > 0) llr += c * @log(c);
    if (d > 0) llr += d * @log(d);
    if (n > 0) llr += n * @log(n);
    if (a + b > 0) llr -= (a + b) * @log(a + b);
    if (c + d > 0) llr -= (c + d) * @log(c + d);
    if (a + c > 0) llr -= (a + c) * @log(a + c);
    if (b + d > 0) llr -= (b + d) * @log(b + d);

    // G = 2 * LLR (standard G-test statistic)
    const g = 2.0 * llr;

    // P-value from chi-squared with 1 df: P(chi^2 > G) = Q(0.5, G/2) = 1 - P(0.5, G/2)
    const p = 1.0 - gamma_mod.incompleteGamma(0.5, g / 2.0);

    return .{ .g = g, .p = p };
}

/// Result of linear regression y = a + bx.
pub const LinearRegressionResult = struct {
    a: f64, // intercept
    b: f64, // slope
    sigma_a: f64, // std error of a
    sigma_b: f64, // std error of b
    cov_ab: f64, // covariance of a and b
    cc: f64, // Pearson correlation coefficient
    q: ?f64, // chi-squared goodness-of-fit P-value (only when sigma is provided)
};

/// Fit n points (x,y) to a straight line y = a + bx by linear regression.
/// sigma is optional per-point standard deviations for y (null = uniform weight).
pub fn linearRegression(x: []const f64, y: []const f64, sigma: ?[]const f64) !LinearRegressionResult {
    const n = x.len;
    if (n <= 2) return error.InvalidInput;
    if (y.len != n) return error.InvalidInput;
    if (sigma) |s| {
        if (s.len != n) return error.InvalidInput;
    }

    // S = sum(1/sigma_i^2)
    var big_s: f64 = 0;
    for (0..n) |i| {
        if (sigma) |s| {
            big_s += 1.0 / (s[i] * s[i]);
        } else {
            big_s += 1.0;
        }
    }

    // Sx = sum(x_i / sigma_i^2)
    var sx: f64 = 0;
    for (0..n) |i| {
        if (sigma) |s| {
            sx += x[i] / (s[i] * s[i]);
        } else {
            sx += x[i];
        }
    }

    // Sy = sum(y_i / sigma_i^2)
    var sy: f64 = 0;
    for (0..n) |i| {
        if (sigma) |s| {
            sy += y[i] / (s[i] * s[i]);
        } else {
            sy += y[i];
        }
    }

    // Stt = sum(t_i^2) where t_i = (x_i - Sx/S) / sigma_i
    var stt: f64 = 0;
    var b_val: f64 = 0;
    for (0..n) |i| {
        var t = x[i] - sx / big_s;
        if (sigma) |s| t /= s[i];
        stt += t * t;
        if (sigma) |s| {
            b_val += t * y[i] / s[i];
        } else {
            b_val += t * y[i];
        }
    }
    if (stt == 0) return error.InvalidInput;

    b_val /= stt;
    const a_val = (sy - sx * b_val) / big_s;
    const sa = @sqrt((1.0 + (sx * sx) / (big_s * stt)) / big_s);
    const sb = @sqrt(1.0 / stt);
    const cov = -sx / (big_s * stt);

    // Pearson correlation coefficient
    var sxy: f64 = 0;
    var sxx: f64 = 0;
    var syy: f64 = 0;
    const mean_x = sx / big_s;
    const mean_y = sy / big_s;
    for (0..n) |i| {
        const xdev = x[i] - mean_x;
        const ydev = y[i] - mean_y;
        sxy += xdev * ydev;
        sxx += xdev * xdev;
        syy += ydev * ydev;
    }
    const cc = if (sxx > 0 and syy > 0) sxy / @sqrt(sxx * syy) else 0;

    // Chi-squared goodness-of-fit P-value (only meaningful with per-point sigma)
    const q_val: ?f64 = if (sigma) |s| blk: {
        var chi2: f64 = 0;
        for (0..n) |i| {
            const residual = y[i] - a_val - b_val * x[i];
            chi2 += (residual * residual) / (s[i] * s[i]);
        }
        const df: f64 = @floatFromInt(n - 2);
        break :blk 1.0 - gamma_mod.incompleteGamma(df / 2.0, chi2 / 2.0);
    } else null;

    return LinearRegressionResult{
        .a = a_val,
        .b = b_val,
        .sigma_a = sa,
        .sigma_b = sb,
        .cov_ab = cov,
        .cc = cc,
        .q = q_val,
    };
}

/// Compute mean and sample variance of a slice.
pub fn meanVariance(x: []const f64) struct { mean: f64, variance: f64 } {
    if (x.len == 0) return .{ .mean = 0, .variance = 0 };

    var sum: f64 = 0;
    var sqsum: f64 = 0;
    for (x) |v| {
        sum += v;
        sqsum += v * v;
    }
    const n: f64 = @floatFromInt(x.len);
    const mean_val = sum / n;
    const var_val = if (x.len > 1)
        @abs((sqsum - sum * sum / n) / (n - 1.0))
    else
        0.0;

    return .{ .mean = mean_val, .variance = var_val };
}

// --- Tests ---

test "logGamma: known values" {
    // log(Gamma(1)) = log(1) = 0
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), logGamma(1.0), 1e-6);
    // log(Gamma(2)) = log(1) = 0
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), logGamma(2.0), 1e-6);
    // log(Gamma(0.5)) = log(sqrt(pi)) ≈ 0.5723649
    try std.testing.expectApproxEqAbs(@as(f64, 0.5723649), logGamma(0.5), 1e-4);
    // log(Gamma(5)) = log(24) ≈ 3.1780538
    try std.testing.expectApproxEqAbs(@as(f64, 3.1780538), logGamma(5.0), 1e-4);
}

test "psi: known values" {
    // Psi(1) = -gamma ≈ -0.5772
    try std.testing.expectApproxEqAbs(@as(f64, -0.5772156649), psi(1.0), 1e-6);
    // Psi(2) = 1 - gamma ≈ 0.4228
    try std.testing.expectApproxEqAbs(@as(f64, 0.4227843351), psi(2.0), 1e-6);
}

test "trigamma: known values" {
    // Psi'(1) = pi^2/6 ≈ 1.6449
    try std.testing.expectApproxEqAbs(@as(f64, math.pi * math.pi / 6.0), trigamma(1.0), 1e-3);
}

test "erfc: known values" {
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), erfc(0.0), 1e-9);
    try std.testing.expect(erfc(3.0) < 0.01);
}

test "chiSquaredTest: 1 df" {
    // chi^2 = 3.841 at p=0.05 for 1 df
    const p = chiSquaredTest(1, 3.841);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), p, 0.01);
}

test "gTest: identical proportions" {
    // Same proportion: G should be ~0
    const result = gTest(50, 100, 50, 100);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.g, 1e-9);
}

test "gTest: different proportions" {
    const result = gTest(80, 100, 20, 100);
    try std.testing.expect(result.g > 0);
    try std.testing.expect(result.p < 0.01);
}

test "linearRegression: perfect line y = 2x + 1" {
    const x = [_]f64{ 1, 2, 3, 4, 5 };
    const y = [_]f64{ 3, 5, 7, 9, 11 };

    const result = try linearRegression(&x, &y, null);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result.a, 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), result.b, 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result.cc, 1e-9);
    // q is null when sigma is not provided
    try std.testing.expect(result.q == null);
}

test "linearRegression: chi-squared Q-value with sigma" {
    const x = [_]f64{ 1, 2, 3, 4, 5 };
    const y = [_]f64{ 3.1, 4.9, 7.2, 8.8, 11.1 };
    const sigma = [_]f64{ 0.5, 0.5, 0.5, 0.5, 0.5 };

    const result = try linearRegression(&x, &y, &sigma);
    // q should be present when sigma is provided
    try std.testing.expect(result.q != null);
    // For data close to a perfect line with reasonable sigma, q should be large (good fit)
    try std.testing.expect(result.q.? > 0.01);
    try std.testing.expect(result.q.? <= 1.0);
}

test "linearRegression: too few points" {
    const x = [_]f64{ 1, 2 };
    const y = [_]f64{ 3, 5 };
    try std.testing.expectError(error.InvalidInput, linearRegression(&x, &y, null));
}

test "meanVariance: basic" {
    const x = [_]f64{ 2, 4, 6, 8 };
    const result = meanVariance(&x);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result.mean, 1e-9);
    // Variance = sum((x-mean)^2)/(n-1) = (9+1+1+9)/3 = 20/3 ≈ 6.667
    try std.testing.expectApproxEqAbs(@as(f64, 20.0 / 3.0), result.variance, 1e-9);
}

test "meanVariance: single value" {
    const x = [_]f64{42.0};
    const result = meanVariance(&x);
    try std.testing.expectApproxEqAbs(@as(f64, 42.0), result.mean, 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.variance, 1e-9);
}

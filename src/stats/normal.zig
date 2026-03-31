/// Normal (Gaussian) distribution.
///
/// Parameters:
///   mu    - location (mean)
///   sigma - scale (standard deviation, must be > 0)
///
/// PDF: (1/(sigma*sqrt(2*pi))) * exp(-(x-mu)^2 / (2*sigma^2))
const std = @import("std");
const math = std.math;
const c_math = @cImport(@cInclude("math.h"));

/// PDF: (1/(sigma*sqrt(2*pi))) * exp(-(x-mu)^2 / (2*sigma^2))
pub fn pdf(x: f64, mu: f64, sigma: f64) f64 {
    const z = (x - mu) / sigma;
    return @exp(-0.5 * z * z) / (sigma * @sqrt(2.0 * math.pi));
}

/// Log PDF: -0.5*log(2*pi) - log(sigma) - (x-mu)^2/(2*sigma^2)
pub fn logPdf(x: f64, mu: f64, sigma: f64) f64 {
    const z = (x - mu) / sigma;
    return -0.5 * @log(2.0 * math.pi) - @log(sigma) - 0.5 * z * z;
}

/// CDF using erfc to avoid cancellation in tails.
/// CDF(x) = 0.5 * erfc(-(x-mu) / (sigma*sqrt(2)))
pub fn cdf(x: f64, mu: f64, sigma: f64) f64 {
    const z = (x - mu) / (sigma * @sqrt(2.0));
    return 0.5 * c_math.erfc(-z);
}

/// Survival function using erfc to avoid cancellation in tails.
/// surv(x) = 0.5 * erfc((x-mu) / (sigma*sqrt(2)))
pub fn surv(x: f64, mu: f64, sigma: f64) f64 {
    const z = (x - mu) / (sigma * @sqrt(2.0));
    return 0.5 * c_math.erfc(z);
}

// Rational approximation coefficients for the standard normal inverse CDF.
// Peter Acklam's algorithm; maximum relative error ~1.15e-9.
// Reference: https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/

// Numerator coefficients for the central region
const na1 = -3.969683028665376e+01;
const na2 =  2.209460984245205e+02;
const na3 = -2.759285104469687e+02;
const na4 =  1.383577518672690e+02;
const na5 = -3.066479806614716e+01;
const na6 =  2.506628277459239e+00;

// Denominator coefficients for the central region
const nb1 = -5.447609879822406e+01;
const nb2 =  1.615858368580409e+02;
const nb3 = -1.556989798598866e+02;
const nb4 =  6.680131188771972e+01;
const nb5 = -1.328068155288572e+01;

// Numerator coefficients for the tail region
const nc1 = -7.784894002430293e-03;
const nc2 = -3.223964580411365e-01;
const nc3 = -2.400758277161838e+00;
const nc4 = -2.549732539343734e+00;
const nc5 =  4.374664141464968e+00;
const nc6 =  2.938163982698783e+00;

// Denominator coefficients for the tail region
const nd1 =  7.784695709041462e-03;
const nd2 =  3.224671290700398e-01;
const nd3 =  2.445134137142996e+00;
const nd4 =  3.754408661907416e+00;

const p_low  = 0.02425;
const p_high = 1.0 - p_low;

/// Inverse CDF of the standard normal distribution (mean=0, sigma=1).
/// Uses Peter Acklam's rational approximation; max relative error ~1.15e-9.
fn stdNormalInv(p: f64) f64 {
    if (p <= 0.0) return -math.inf(f64);
    if (p >= 1.0) return math.inf(f64);

    if (p < p_low) {
        // Lower tail
        const q = @sqrt(-2.0 * @log(p));
        return (((((nc1*q+nc2)*q+nc3)*q+nc4)*q+nc5)*q+nc6) /
               ((((nd1*q+nd2)*q+nd3)*q+nd4)*q+1.0);
    } else if (p <= p_high) {
        // Central region
        const q = p - 0.5;
        const r = q * q;
        return (((((na1*r+na2)*r+na3)*r+na4)*r+na5)*r+na6)*q /
               (((((nb1*r+nb2)*r+nb3)*r+nb4)*r+nb5)*r+1.0);
    } else {
        // Upper tail: use symmetry
        const q = @sqrt(-2.0 * @log(1.0 - p));
        return -(((((nc1*q+nc2)*q+nc3)*q+nc4)*q+nc5)*q+nc6) /
                ((((nd1*q+nd2)*q+nd3)*q+nd4)*q+1.0);
    }
}

/// Inverse CDF (quantile): mu + sigma * Phi^{-1}(p)
pub fn invcdf(p: f64, mu: f64, sigma: f64) f64 {
    return mu + sigma * stdNormalInv(p);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "normal pdf at x=0, mu=0, sigma=1" {
    // pdf = 1/sqrt(2*pi) ≈ 0.39894228
    const result = pdf(0.0, 0.0, 1.0);
    const expected = 1.0 / @sqrt(2.0 * math.pi);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "normal logPdf at x=0, mu=0, sigma=1" {
    const result = logPdf(0.0, 0.0, 1.0);
    const expected = @log(pdf(0.0, 0.0, 1.0));
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-10));
}

test "normal cdf at mean is 0.5" {
    try std.testing.expect(math.approxEqAbs(f64, cdf(0.0, 0.0, 1.0), 0.5, 1e-7));
    try std.testing.expect(math.approxEqAbs(f64, cdf(5.0, 5.0, 2.0), 0.5, 1e-7));
}

test "normal cdf at 1.96 sigma is approx 0.975" {
    const result = cdf(1.96, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, 0.975, 1e-4));
}

test "normal surv at mean is 0.5" {
    try std.testing.expect(math.approxEqAbs(f64, surv(0.0, 0.0, 1.0), 0.5, 1e-7));
}

test "normal cdf + surv = 1" {
    const xs = [_]f64{ -3.0, -1.0, 0.0, 1.0, 3.0 };
    for (xs) |x| {
        const c = cdf(x, 0.0, 1.0);
        const s = surv(x, 0.0, 1.0);
        try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-14));
    }
}

test "normal invcdf at 0.5 is mean" {
    try std.testing.expect(math.approxEqAbs(f64, invcdf(0.5, 0.0, 1.0), 0.0, 1e-4));
    try std.testing.expect(math.approxEqAbs(f64, invcdf(0.5, 3.0, 2.0), 3.0, 1e-4));
}

test "normal invcdf at 0.975 is approx 1.96" {
    const result = invcdf(0.975, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, 1.96, 1e-2));
}

test "normal invcdf round-trip" {
    const mu: f64 = 1.0;
    const sigma: f64 = 2.0;
    const xs = [_]f64{ -2.0, 0.0, 1.0, 3.0, 5.0 };
    for (xs) |x| {
        const p = cdf(x, mu, sigma);
        const x2 = invcdf(p, mu, sigma);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-4));
    }
}

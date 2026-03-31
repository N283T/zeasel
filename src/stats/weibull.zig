/// Weibull distribution.
///
/// Parameters:
///   mu     - location parameter
///   lambda - scale parameter (must be > 0)
///   tau    - shape parameter (must be > 0)
///
/// Defined for x >= mu.
/// PDF: (tau * lambda) * (lambda*(x-mu))^(tau-1) * exp(-(lambda*(x-mu))^tau)
const std = @import("std");
const math = std.math;
const minimizer = @import("minimizer.zig");

/// PDF: (tau * lambda) * (lambda*(x-mu))^(tau-1) * exp(-(lambda*(x-mu))^tau)  for x >= mu, else 0.
pub fn pdf(x: f64, mu: f64, lambda: f64, tau: f64) f64 {
    if (x < mu) return 0.0;
    const u = lambda * (x - mu);
    if (u == 0.0) {
        // Limit depends on tau: pdf(mu) = tau*lambda when tau=1, 0 when tau>1, inf when tau<1
        if (tau < 1.0) return math.inf(f64);
        if (tau > 1.0) return 0.0;
        return lambda; // tau == 1: pure exponential
    }
    return (tau * lambda) * math.pow(f64, u, tau - 1.0) * @exp(-math.pow(f64, u, tau));
}

/// CDF: 1 - exp(-(lambda*(x-mu))^tau)  for x >= mu, else 0.
pub fn cdf(x: f64, mu: f64, lambda: f64, tau: f64) f64 {
    if (x < mu) return 0.0;
    const u = lambda * (x - mu);
    return 1.0 - @exp(-math.pow(f64, u, tau));
}

/// Survival function: exp(-(lambda*(x-mu))^tau)  for x >= mu, else 1.
pub fn surv(x: f64, mu: f64, lambda: f64, tau: f64) f64 {
    if (x < mu) return 1.0;
    const u = lambda * (x - mu);
    return @exp(-math.pow(f64, u, tau));
}

/// Inverse CDF: mu + (1/lambda) * (-log(1-p))^(1/tau)
pub fn invcdf(p: f64, mu: f64, lambda: f64, tau: f64) f64 {
    return mu + (1.0 / lambda) * math.pow(f64, -@log(1.0 - p), 1.0 / tau);
}

/// Log-PDF for fitting: log(tau*lambda) + (tau-1)*log(lambda*(x-mu)) - (lambda*(x-mu))^tau
/// Returns -inf for x at mu when tau < 1 (skip these in fitting, following Easel).
fn logPdf(x: f64, mu: f64, lambda: f64, tau: f64) f64 {
    if (x < mu) return -math.inf(f64);
    const u = lambda * (x - mu);
    if (u == 0.0) {
        if (tau <= 1.0) return -math.inf(f64);
        return -math.inf(f64); // log(0) for tau > 1
    }
    return @log(tau * lambda) + (tau - 1.0) * @log(u) - math.pow(f64, u, tau);
}

/// Data passed to the Weibull NLL objective/gradient functions.
const WeibullFitData = struct {
    x: []const f64,
    mu: f64,
};

/// Negative log-likelihood objective for Weibull fitting.
/// Parameters p[0] = log(lambda), p[1] = log(tau) (change of variables for positivity).
fn weibullObjective(p: []const f64, user_data: ?*anyopaque) f64 {
    const data: *const WeibullFitData = @ptrCast(@alignCast(user_data.?));
    const lambda = @exp(p[0]);
    const tau = @exp(p[1]);

    var nll: f64 = 0;
    for (data.x) |xi| {
        if (tau < 1.0 and xi == data.mu) continue; // skip infinity, following Easel
        nll -= logPdf(xi, data.mu, lambda, tau);
    }
    return nll;
}

/// Gradient of the negative log-likelihood w.r.t. log(lambda) and log(tau).
/// Uses chain rule: d/d(log(lambda)) = lambda * d/d(lambda), etc.
fn weibullGradient(p: []const f64, grad: []f64, user_data: ?*anyopaque) void {
    const data: *const WeibullFitData = @ptrCast(@alignCast(user_data.?));
    const lambda = @exp(p[0]);
    const tau = @exp(p[1]);

    var dl: f64 = 0; // d(NLL)/d(lambda)
    var dt: f64 = 0; // d(NLL)/d(tau)

    for (data.x) |xi| {
        if (tau < 1.0 and xi == data.mu) continue;
        const shifted = xi - data.mu;
        if (shifted == 0.0) continue;
        const u = lambda * shifted;
        const u_tau = math.pow(f64, u, tau);

        // d(logpdf)/d(lambda) = 1/lambda + (tau-1)/lambda - tau * u^(tau-1) * shifted
        //                     = tau/lambda * (1 - u^tau)
        dl += tau / lambda * (1.0 - u_tau);

        // d(logpdf)/d(tau) = 1/tau + log(u) - u^tau * log(u)
        dt += 1.0 / tau + @log(u) * (1.0 - u_tau);
    }

    // Chain rule: d/d(log(lambda)) = lambda * d/d(lambda)
    grad[0] = -dl * lambda;
    grad[1] = -dt * tau;
}

/// Result of fitting a Weibull distribution.
pub const FitResult = struct {
    mu: f64,
    lambda: f64,
    tau: f64,
};

/// Maximum likelihood fit of a complete Weibull distribution.
///
/// Given observed data x[0..n-1], estimates mu (location), lambda (scale),
/// and tau (shape). mu is set to the minimum observed value.
/// lambda and tau are fitted by minimizing the negative log-likelihood
/// using conjugate gradient descent.
///
/// Corresponds to Easel's esl_wei_FitComplete().
pub fn fitComplete(allocator: std.mem.Allocator, x: []const f64) !FitResult {
    std.debug.assert(x.len > 0);

    // mu = min(x)
    var mu = x[0];
    for (x[1..]) |xi| {
        if (xi < mu) mu = xi;
    }

    // Initial guesses based on exponential fit
    var mean_val: f64 = 0;
    for (x) |xi| mean_val += xi;
    mean_val /= @as(f64, @floatFromInt(x.len));
    const lambda_init = 1.0 / (mean_val - mu);
    const tau_init: f64 = 0.9;

    // Change of variables: p[0] = log(lambda), p[1] = log(tau)
    var p = [_]f64{ @log(lambda_init), @log(tau_init) };

    var data = WeibullFitData{
        .x = x,
        .mu = mu,
    };

    _ = try minimizer.minimize(
        allocator,
        &p,
        weibullObjective,
        weibullGradient,
        @ptrCast(&data),
        .{ .max_iter = 200, .tol = 1e-8 },
    );

    return .{
        .mu = mu,
        .lambda = @exp(p[0]),
        .tau = @exp(p[1]),
    };
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "weibull cdf at x=mu is 0" {
    try std.testing.expectEqual(@as(f64, 0.0), cdf(0.0, 0.0, 1.0, 1.0));
    try std.testing.expectEqual(@as(f64, 0.0), cdf(2.0, 2.0, 1.5, 2.0));
}

test "weibull cdf below mu is 0" {
    try std.testing.expectEqual(@as(f64, 0.0), cdf(-1.0, 0.0, 1.0, 2.0));
}

test "weibull surv below mu is 1" {
    try std.testing.expectEqual(@as(f64, 1.0), surv(-1.0, 0.0, 1.0, 2.0));
}

test "weibull cdf tau=1 is exponential special case" {
    // When tau=1, Weibull reduces to Exponential(lambda)
    // cdf(1, 0, 1, 1) = 1 - exp(-1) ≈ 0.6321
    const result = cdf(1.0, 0.0, 1.0, 1.0);
    const expected = 1.0 - @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-10));
}

test "weibull pdf at known value tau=1" {
    // pdf(1, 0, 1, 1) = lambda * exp(-1) = exp(-1) ≈ 0.36787944
    const result = pdf(1.0, 0.0, 1.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, @exp(-1.0), 1e-10));
}

test "weibull pdf at known value tau=2" {
    // pdf(1, 0, 1, 2) = 2*1*(1)^1 * exp(-1) = 2*exp(-1)
    const result = pdf(1.0, 0.0, 1.0, 2.0);
    const expected = 2.0 * @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-10));
}

test "weibull pdf below mu is zero" {
    try std.testing.expectEqual(@as(f64, 0.0), pdf(-1.0, 0.0, 1.0, 2.0));
}

test "weibull cdf + surv = 1" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const c = cdf(x, 0.0, 1.0, 2.0);
        const s = surv(x, 0.0, 1.0, 2.0);
        try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-14));
    }
}

test "weibull invcdf round-trip" {
    const mu: f64 = 1.0;
    const lambda: f64 = 0.5;
    const tau: f64 = 2.0;
    const xs = [_]f64{ 1.5, 2.0, 3.0, 5.0 };
    for (xs) |x| {
        const p = cdf(x, mu, lambda, tau);
        const x2 = invcdf(p, mu, lambda, tau);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-10));
    }
}

test "weibull invcdf tau=1 round-trip" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 4.0 };
    for (xs) |x| {
        const p = cdf(x, 0.0, 1.0, 1.0);
        const x2 = invcdf(p, 0.0, 1.0, 1.0);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-10));
    }
}

test "fitComplete: weibull with known parameters" {
    const allocator = std.testing.allocator;

    // Data drawn from Weibull(mu=0, lambda=1, tau=2) -- Rayleigh-like
    // Generated analytically: x = (-log(1-p))^(1/tau) / lambda
    // with p = 0.05, 0.15, ..., 0.95
    const data = [_]f64{
        invcdf(0.05, 0.0, 1.0, 2.0),
        invcdf(0.15, 0.0, 1.0, 2.0),
        invcdf(0.25, 0.0, 1.0, 2.0),
        invcdf(0.35, 0.0, 1.0, 2.0),
        invcdf(0.45, 0.0, 1.0, 2.0),
        invcdf(0.55, 0.0, 1.0, 2.0),
        invcdf(0.65, 0.0, 1.0, 2.0),
        invcdf(0.75, 0.0, 1.0, 2.0),
        invcdf(0.85, 0.0, 1.0, 2.0),
        invcdf(0.95, 0.0, 1.0, 2.0),
    };

    const result = try fitComplete(allocator, &data);

    // mu should be min(data) = invcdf(0.05, 0, 1, 2)
    try std.testing.expectApproxEqAbs(data[0], result.mu, 1e-10);

    // lambda and tau should be in reasonable range
    // With only 10 quantile-spaced points, we accept generous tolerance
    try std.testing.expect(result.lambda > 0.3);
    try std.testing.expect(result.lambda < 3.0);
    try std.testing.expect(result.tau > 0.5);
    try std.testing.expect(result.tau < 5.0);
}

test "fitComplete: weibull exponential special case (tau~1)" {
    const allocator = std.testing.allocator;

    // Weibull with tau=1 is exponential(lambda). Use lambda=0.5.
    const data = [_]f64{
        invcdf(0.05, 0.0, 0.5, 1.0),
        invcdf(0.15, 0.0, 0.5, 1.0),
        invcdf(0.25, 0.0, 0.5, 1.0),
        invcdf(0.35, 0.0, 0.5, 1.0),
        invcdf(0.45, 0.0, 0.5, 1.0),
        invcdf(0.55, 0.0, 0.5, 1.0),
        invcdf(0.65, 0.0, 0.5, 1.0),
        invcdf(0.75, 0.0, 0.5, 1.0),
        invcdf(0.85, 0.0, 0.5, 1.0),
        invcdf(0.95, 0.0, 0.5, 1.0),
    };

    const result = try fitComplete(allocator, &data);

    // tau should be close to 1
    try std.testing.expect(result.tau > 0.5);
    try std.testing.expect(result.tau < 2.0);
    // lambda should be close to 0.5
    try std.testing.expect(result.lambda > 0.1);
    try std.testing.expect(result.lambda < 2.0);
}

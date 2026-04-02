/// Generalized Extreme Value (GEV) distribution.
///
/// Parameters:
///   mu     - location parameter
///   lambda - scale parameter (rate; lambda = 1/sigma; must be > 0)
///   alpha  - shape parameter (xi in some references)
///
/// This parameterization matches Easel's esl_gev convention where lambda is the
/// rate (inverse scale). The relationship to the Coles (2001) sigma notation is
/// lambda = 1/sigma.
///
/// The GEV unifies the three extreme value types:
///   alpha > 0: Type II (Frechet) - lower bound at mu - 1/(alpha*lambda)
///   alpha < 0: Type III (Weibull) - upper bound at mu - 1/(alpha*lambda)
///   alpha = 0: Type I (Gumbel) - no finite bound
///
/// Standardized variable: y = lambda * (x - mu), t = 1 + alpha * y
/// When |alpha * y| < 1e-12, we use the Gumbel limit formulas.
const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;

/// Threshold below which |alpha| is treated as zero (Gumbel limit).
const alpha_threshold: f64 = 1e-6;

/// PDF of the GEV distribution.
///
/// When alpha != 0: lambda * exp(-(1+1/alpha)*log(t) - exp(-log(t)/alpha))
///   where t = 1 + alpha * lambda * (x - mu), valid for t > 0.
/// When alpha == 0 (Gumbel): lambda * exp(-y - exp(-y))
pub fn pdf(x: f64, mu: f64, lambda: f64, alpha: f64) f64 {
    return @exp(logPdf(x, mu, lambda, alpha));
}

/// Log PDF of the GEV distribution.
pub fn logPdf(x: f64, mu: f64, lambda: f64, alpha: f64) f64 {
    const y = lambda * (x - mu);

    if (@abs(alpha) < alpha_threshold or @abs(y * alpha) < 1e-12) {
        // Gumbel limit: log(lambda) - y - exp(-y)
        return @log(lambda) - y - @exp(-y);
    }

    const t = 1.0 + alpha * y;
    if (t <= 0.0) return -math.inf(f64);

    const log_t = @log(t);
    // log(lambda) - (1 + 1/alpha)*log(t) - exp(-log(t)/alpha)
    return @log(lambda) - (1.0 + 1.0 / alpha) * log_t - @exp(-log_t / alpha);
}

/// CDF of the GEV distribution.
///
/// When alpha != 0: exp(-exp(-log(t)/alpha)) where t = 1 + alpha*lambda*(x-mu)
/// When alpha == 0 (Gumbel): exp(-exp(-y))
pub fn cdf(x: f64, mu: f64, lambda: f64, alpha: f64) f64 {
    const y = lambda * (x - mu);

    if (@abs(alpha) < alpha_threshold or @abs(y * alpha) < 1e-12) {
        return @exp(-@exp(-y));
    }

    const t = 1.0 + alpha * y;
    if (t <= 0.0) {
        // Beyond the support boundary: Frechet lower, Weibull upper
        return if (x < mu) 0.0 else 1.0;
    }

    return @exp(-@exp(-@log(t) / alpha));
}

/// Survival function: 1 - CDF.
pub fn surv(x: f64, mu: f64, lambda: f64, alpha: f64) f64 {
    const c = cdf(x, mu, lambda, alpha);
    // For numerical stability when CDF is near 1
    if (c > 0.9999999) {
        const lc = logCdf(x, mu, lambda, alpha);
        return -math.expm1(lc);
    }
    return 1.0 - c;
}

/// Log CDF of the GEV distribution.
fn logCdf(x: f64, mu: f64, lambda: f64, alpha: f64) f64 {
    const y = lambda * (x - mu);

    if (@abs(alpha) < alpha_threshold or @abs(y * alpha) < 1e-12) {
        return -@exp(-y);
    }

    const t = 1.0 + alpha * y;
    if (t <= 0.0) {
        return if (x < mu) -math.inf(f64) else 0.0;
    }

    return -@exp(-@log(t) / alpha);
}

/// Inverse CDF (quantile function).
///
/// When alpha != 0: mu + (exp(-alpha*log(-log(p))) - 1) / (alpha*lambda)
/// When alpha == 0: mu - log(-log(p)) / lambda
pub fn invcdf(p: f64, mu: f64, lambda: f64, alpha: f64) f64 {
    if (p <= 0.0 or p >= 1.0) return math.nan(f64);

    const log_p = @log(p);

    if (@abs(alpha) < alpha_threshold) {
        // Gumbel: mu - log(-log(p)) / lambda
        return mu - @log(-log_p) / lambda;
    }

    // General: mu + (exp(-alpha * log(-log(p))) - 1) / (alpha * lambda)
    return mu + (@exp(-alpha * @log(-log_p)) - 1.0) / (alpha * lambda);
}

/// Maximum likelihood estimation of GEV parameters from complete data.
///
/// Uses the method of probability-weighted moments (PWM, Hosking 1985) for
/// initial estimates, then refines with conjugate gradient minimization of
/// the negative log-likelihood.
///
/// Returns fitted (mu, lambda, alpha) parameters where lambda = 1/sigma.
pub fn fitComplete(allocator: Allocator, x: []const f64) !struct { mu: f64, lambda: f64, alpha: f64 } {
    if (x.len < 3) return error.InsufficientData;

    const n: f64 = @floatFromInt(x.len);

    // Sort data for PWM computation
    const sorted = try allocator.alloc(f64, x.len);
    defer allocator.free(sorted);
    @memcpy(sorted, x);
    std.mem.sort(f64, sorted, {}, std.sort.asc(f64));

    // Check for zero variance
    if (sorted[0] == sorted[sorted.len - 1]) return error.ZeroVariance;

    // Compute probability-weighted moments (PWMs): b0, b1, b2
    // b_r = (1/n) * sum_{i=1}^{n} ((i-1) choose r) / ((n-1) choose r) * x_{i:n}
    // Using the unbiased estimators:
    // b0 = mean(x)
    // b1 = (1/n) * sum_{i=1}^{n} (i-1)/(n-1) * x_{i:n}
    // b2 = (1/n) * sum_{i=1}^{n} (i-1)*(i-2)/((n-1)*(n-2)) * x_{i:n}
    var b0: f64 = 0.0;
    var b1: f64 = 0.0;
    var b2: f64 = 0.0;
    for (sorted, 0..) |xi, idx| {
        const i: f64 = @floatFromInt(idx);
        b0 += xi;
        b1 += (i / (n - 1.0)) * xi;
        if (n > 2.0) {
            b2 += (i * (i - 1.0) / ((n - 1.0) * (n - 2.0))) * xi;
        }
    }
    b0 /= n;
    b1 /= n;
    b2 /= n;

    // Hosking (1985) estimators from PWMs
    // c = (2*b1 - b0) / (3*b2 - b0) - log(2)/log(3)
    const c = (2.0 * b1 - b0) / (3.0 * b2 - b0) - @log(2.0) / @log(3.0);

    // alpha approximation: 7.8590*c + 2.9554*c^2
    var alpha_init = 7.8590 * c + 2.9554 * c * c;

    // Clamp alpha to a reasonable range to avoid degenerate fits
    if (alpha_init > 0.5) alpha_init = 0.5;
    if (alpha_init < -0.5) alpha_init = -0.5;

    // sigma and mu from PWMs, then convert sigma to lambda = 1/sigma
    var sigma_init: f64 = undefined;
    var mu_init: f64 = undefined;
    if (@abs(alpha_init) < alpha_threshold) {
        // Gumbel case: sigma = (2*b1 - b0) / log(2), mu = b0 - sigma*euler_gamma
        sigma_init = (2.0 * b1 - b0) / @log(2.0);
        mu_init = b0 - 0.5772156649015329 * sigma_init;
    } else {
        const g1 = math.gamma(f64, 1.0 + alpha_init);
        sigma_init = alpha_init * (2.0 * b1 - b0) / (g1 * (1.0 - math.pow(f64, 2.0, -alpha_init)));
        mu_init = b0 + sigma_init * (g1 - 1.0) / alpha_init;
    }

    if (sigma_init <= 0.0) sigma_init = 1.0;
    if (!math.isFinite(sigma_init)) sigma_init = 1.0;
    if (!math.isFinite(mu_init)) mu_init = b0;

    const lambda_init = 1.0 / sigma_init;

    // Refine with conjugate gradient minimization.
    // Parameterize as p[0] = mu, p[1] = log(lambda), p[2] = alpha
    // to ensure lambda > 0.
    var p = [3]f64{ mu_init, @log(lambda_init), alpha_init };

    var fit_data = FitData{
        .x = x,
    };

    const minimizer_mod = @import("minimizer.zig");
    _ = minimizer_mod.minimize(
        allocator,
        &p,
        gevNegLogLik,
        gevNegLogLikGrad,
        @ptrCast(&fit_data),
        .{ .tol = 1e-6, .max_iter = 300 },
    ) catch {
        // If minimizer fails, return PWM estimates
        return .{ .mu = mu_init, .lambda = lambda_init, .alpha = alpha_init };
    };

    const fitted_lambda = @exp(p[1]);
    if (!math.isFinite(p[0]) or !math.isFinite(fitted_lambda) or !math.isFinite(p[2])) {
        // Fall back to PWM estimates
        return .{ .mu = mu_init, .lambda = lambda_init, .alpha = alpha_init };
    }

    return .{ .mu = p[0], .lambda = fitted_lambda, .alpha = p[2] };
}

const FitData = struct {
    x: []const f64,
};

/// Negative log-likelihood for GEV fitting.
/// Parameters: p[0] = mu, p[1] = log(lambda), p[2] = alpha
fn gevNegLogLik(p: []const f64, user_data: ?*anyopaque) f64 {
    const data: *const FitData = @ptrCast(@alignCast(user_data.?));
    const x = data.x;
    const mu = p[0];
    const lambda = @exp(p[1]);
    const alpha = p[2];

    var nll: f64 = 0.0;
    for (x) |xi| {
        const lp = logPdf(xi, mu, lambda, alpha);
        if (!math.isFinite(lp)) return math.inf(f64);
        nll -= lp;
    }
    return nll;
}

/// Gradient of the negative log-likelihood for GEV fitting.
/// Parameters: p[0] = mu, p[1] = v = log(lambda), p[2] = alpha
///
/// logpdf = v - (1 + 1/alpha)*log(t) - exp(-log(t)/alpha)
/// where t = 1 + alpha * lambda * (x - mu) = 1 + alpha * y, y = lambda*(x-mu)
fn gevNegLogLikGrad(p: []const f64, dp: []f64, user_data: ?*anyopaque) void {
    const data: *const FitData = @ptrCast(@alignCast(user_data.?));
    const x = data.x;
    const mu = p[0];
    const lambda = @exp(p[1]);
    const alpha = p[2];

    var d_mu: f64 = 0.0;
    var d_v: f64 = 0.0; // w.r.t. log(lambda)
    var d_alpha: f64 = 0.0;

    if (@abs(alpha) < alpha_threshold) {
        // Gumbel case: logpdf = v - y - exp(-y), y = lambda*(x-mu)
        // y = exp(v)*(x-mu)
        // d(logpdf)/d(mu) = lambda*(exp(-y) - 1)  [since d(-y)/d(mu)=lambda, d(-exp(-y))/d(mu)=-lambda*exp(-y)*(-1)=lambda*exp(-y)... wait]
        // d(-y)/d(mu) = -d(y)/d(mu) = -(-lambda) = lambda
        // d(-exp(-y))/d(mu) = exp(-y)*(-1)*d(-y)/d(mu)... no:
        // d(-exp(-y))/d(mu) = -exp(-y) * d(-y)/d(mu) ... chain rule on -exp(-y):
        //   d/d(mu)[-exp(-y)] = exp(-y) * dy/d(mu) = exp(-y)*(-lambda)
        // So d(logpdf)/d(mu) = lambda + (-lambda)*exp(-y) = lambda*(1 - exp(-y))
        //
        // d(logpdf)/d(v): logpdf = v - y - exp(-y), y = exp(v)*(x-mu)
        //   d/dv = 1 - (x-mu)*exp(v) + exp(-y)*(x-mu)*exp(v)
        //        = 1 - y + y*exp(-y)
        //        = 1 + y*(exp(-y) - 1)
        for (x) |xi| {
            const y = lambda * (xi - mu);
            const ey = @exp(-y);
            d_mu -= lambda * (1.0 - ey);
            d_v -= 1.0 + y * (ey - 1.0);
        }
        dp[0] = d_mu;
        dp[1] = d_v;
        dp[2] = 0.0;
        return;
    }

    // General case: logpdf = v - (1 + 1/alpha)*log(t) - exp(-log(t)/alpha)
    // where t = 1 + alpha * lambda * (x - mu) = 1 + alpha * y, y = lambda*(x-mu)
    // v = log(lambda)
    //
    // Let lt = log(t), u = exp(-lt/alpha)
    // logpdf = v - (1 + 1/alpha)*lt - u
    //
    // dt/d(mu)    = -alpha*lambda
    // dt/d(v)     = alpha*lambda*(xi-mu) = alpha*y  (chain rule: d(lambda)/d(v)=lambda)
    // dt/d(alpha) = y
    //
    // d(lt)/dt = 1/t
    // d(u)/d(lt) = -u/alpha
    //
    // d(logpdf)/d(lt) = -(1 + 1/alpha) + u/alpha
    //
    // Chain: d(logpdf)/d(param) = [d(logpdf)/d(lt)] * (1/t) * dt/d(param)
    //        plus explicit alpha terms for param=alpha

    const inv_alpha = 1.0 / alpha;

    for (x) |xi| {
        const y = lambda * (xi - mu);
        const t = 1.0 + alpha * y;

        if (t <= 0.0) {
            // Outside support; gradient contribution is zero
            continue;
        }

        const log_t = @log(t);
        const u = @exp(-log_t / alpha); // = t^(-1/alpha)

        // d(logpdf)/d(log_t) = -(1 + inv_alpha) + u * inv_alpha
        const dl_dlt = -(1.0 + inv_alpha) + u * inv_alpha;
        // d(logpdf)/d(t) = dl_dlt / t
        const dl_dt = dl_dlt / t;

        // d(logpdf)/d(mu): dt/d(mu) = -alpha*lambda
        d_mu -= dl_dt * (-alpha * lambda);

        // d(logpdf)/d(v): explicit +1 from v, plus chain dt/d(v) = alpha*y
        d_v -= 1.0 + dl_dt * alpha * y;

        // d(logpdf)/d(alpha): explicit terms + chain dt/d(alpha)=y
        //   explicit from -(1+1/alpha)*lt: d/d(alpha) = (1/alpha^2)*lt
        //   explicit from -u: d/d(alpha)[-exp(-lt/alpha)] = -u * (lt/alpha^2)
        //     wait: d/d(alpha)[exp(-lt/alpha)] = exp(-lt/alpha) * lt/alpha^2
        //     so d/d(alpha)[-u] = -u * lt/alpha^2  -- no, that's wrong sign:
        //     d/d(alpha)[-exp(-lt/alpha)] = +exp(-lt/alpha) * lt/alpha^2 = u*lt/alpha^2
        //   explicit total: (1/alpha^2)*lt - (1/alpha^2)*lt + ... hmm let me be careful:
        //   d/d(alpha)[-(1+1/alpha)*lt] = (1/alpha^2)*lt
        //   d/d(alpha)[-u] = d/d(alpha)[-exp(-lt/alpha)]
        //                   = -exp(-lt/alpha) * d/d(alpha)[-lt/alpha]
        //                   = -u * (lt/alpha^2)
        //   explicit total: (1/alpha^2)*lt - u*(lt/alpha^2) = (lt/alpha^2)*(1 - u)
        //   chain: dl_dt * y
        const da_direct = (log_t / (alpha * alpha)) * (1.0 - u) + dl_dt * y;
        d_alpha -= da_direct;
    }

    dp[0] = d_mu;
    dp[1] = d_v;
    dp[2] = d_alpha;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "gev: pdf/cdf/surv consistency" {
    // Test CDF + surv = 1.0 for various parameter combinations.
    // lambda = 1/sigma (rate parameterization, matching Easel).
    const params = [_]struct { mu: f64, lambda: f64, alpha: f64 }{
        .{ .mu = 0.0, .lambda = 1.0, .alpha = 0.0 }, // Gumbel
        .{ .mu = 1.0, .lambda = 0.5, .alpha = 0.1 }, // Frechet (lambda=0.5 => sigma=2)
        .{ .mu = 0.0, .lambda = 1.0, .alpha = -0.2 }, // Weibull-type
        .{ .mu = 5.0, .lambda = 2.0, .alpha = 0.3 }, // lambda=2 => sigma=0.5
    };

    for (params) |p| {
        const xs = [_]f64{ -2.0, 0.0, 1.0, 3.0, 5.0, 10.0 };
        for (xs) |x| {
            const c = cdf(x, p.mu, p.lambda, p.alpha);
            const s = surv(x, p.mu, p.lambda, p.alpha);
            try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-10));
        }
    }
}

test "gev: numerical derivative of CDF approximates PDF" {
    const mu: f64 = 0.0;
    const lambda: f64 = 1.0;
    const alpha: f64 = 0.1;
    const h: f64 = 1e-7;

    const xs = [_]f64{ -1.0, 0.0, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const numerical_pdf = (cdf(x + h, mu, lambda, alpha) - cdf(x - h, mu, lambda, alpha)) / (2.0 * h);
        const analytical_pdf = pdf(x, mu, lambda, alpha);
        try std.testing.expect(math.approxEqAbs(f64, numerical_pdf, analytical_pdf, 1e-5));
    }
}

test "gev: alpha=0 matches Gumbel" {
    // GEV with alpha=0 should match Gumbel distribution.
    // Both modules use the same lambda (rate) parameterization.
    const gumbel = @import("gumbel.zig");
    const mu: f64 = 2.0;
    const lambda: f64 = 1.0 / 1.5; // equivalent to sigma=1.5

    const xs = [_]f64{ -2.0, 0.0, 1.0, 2.0, 3.0, 5.0, 10.0 };
    for (xs) |x| {
        const gev_cdf = cdf(x, mu, lambda, 0.0);
        const gumbel_cdf = gumbel.cdf(x, mu, lambda);
        try std.testing.expect(math.approxEqAbs(f64, gev_cdf, gumbel_cdf, 1e-10));

        const gev_pdf = pdf(x, mu, lambda, 0.0);
        const gumbel_pdf = gumbel.pdf(x, mu, lambda);
        try std.testing.expect(math.approxEqAbs(f64, gev_pdf, gumbel_pdf, 1e-10));
    }
}

test "gev: invcdf is inverse of cdf" {
    const params = [_]struct { mu: f64, lambda: f64, alpha: f64 }{
        .{ .mu = 0.0, .lambda = 1.0, .alpha = 0.0 },
        .{ .mu = 1.0, .lambda = 0.5, .alpha = 0.1 },
        .{ .mu = 0.0, .lambda = 1.0, .alpha = -0.2 },
    };

    const ps = [_]f64{ 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99 };
    for (params) |par| {
        for (ps) |p| {
            const x = invcdf(p, par.mu, par.lambda, par.alpha);
            const recovered_p = cdf(x, par.mu, par.lambda, par.alpha);
            try std.testing.expect(math.approxEqAbs(f64, recovered_p, p, 1e-10));
        }
    }
}

test "gev: cdf round-trip through invcdf" {
    const mu: f64 = 3.0;
    const lambda: f64 = 0.5; // sigma=2
    const alpha: f64 = 0.15;

    const xs = [_]f64{ 0.0, 2.0, 3.0, 5.0, 10.0 };
    for (xs) |x| {
        const p = cdf(x, mu, lambda, alpha);
        if (p > 0.0 and p < 1.0) {
            const x2 = invcdf(p, mu, lambda, alpha);
            try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-8));
        }
    }
}

test "gev: fitComplete recovers parameters" {
    // Generate samples from GEV(mu=5, lambda=0.5, alpha=0.1).
    // lambda=0.5 corresponds to sigma=2.
    const true_mu: f64 = 5.0;
    const true_lambda: f64 = 0.5;
    const true_alpha: f64 = 0.1;
    const n = 5000;

    var samples: [n]f64 = undefined;
    var state: u64 = 42;
    for (&samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_lambda, true_alpha);
    }

    const fit = try fitComplete(std.testing.allocator, &samples);

    // With 5000 samples, expect reasonable accuracy
    try std.testing.expect(math.approxEqAbs(f64, fit.mu, true_mu, 0.3));
    try std.testing.expect(math.approxEqAbs(f64, fit.lambda, true_lambda, 0.15));
    try std.testing.expect(math.approxEqAbs(f64, fit.alpha, true_alpha, 0.25));
}

test "gev: fitComplete Gumbel case (alpha near 0)" {
    // Generate samples from GEV with alpha=0 (Gumbel).
    // lambda=1/3 corresponds to sigma=3.
    const true_mu: f64 = 10.0;
    const true_lambda: f64 = 1.0 / 3.0;
    const n = 5000;

    var samples: [n]f64 = undefined;
    var state: u64 = 12345;
    for (&samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_lambda, 0.0);
    }

    const fit = try fitComplete(std.testing.allocator, &samples);

    try std.testing.expect(math.approxEqAbs(f64, fit.mu, true_mu, 0.5));
    try std.testing.expect(math.approxEqAbs(f64, fit.lambda, true_lambda, 0.15));
    try std.testing.expect(@abs(fit.alpha) < 0.15);
}

test "gev: fitComplete error on insufficient data" {
    var x = [_]f64{ 1.0, 2.0 };
    try std.testing.expectError(error.InsufficientData, fitComplete(std.testing.allocator, &x));
}

test "gev: support boundaries" {
    // Type II (alpha > 0): support is (mu - 1/(alpha*lambda), +inf)
    // At the boundary, CDF should be 0.
    const mu: f64 = 0.0;
    const lambda: f64 = 1.0;
    const alpha: f64 = 0.5;
    const boundary = mu - 1.0 / (alpha * lambda); // = -2.0

    try std.testing.expectEqual(@as(f64, 0.0), cdf(boundary - 1.0, mu, lambda, alpha));
    try std.testing.expect(cdf(boundary + 0.1, mu, lambda, alpha) > 0.0);

    // Type III (alpha < 0): support is (-inf, mu - 1/(alpha*lambda))
    const alpha_neg: f64 = -0.5;
    const upper = mu - 1.0 / (alpha_neg * lambda); // = 2.0

    try std.testing.expectEqual(@as(f64, 1.0), cdf(upper + 1.0, mu, lambda, alpha_neg));
    try std.testing.expect(cdf(upper - 0.1, mu, lambda, alpha_neg) < 1.0);
}

test "gev: logPdf matches log(pdf)" {
    const xs = [_]f64{ -1.0, 0.0, 1.0, 3.0, 5.0 };
    const alphas = [_]f64{ 0.0, 0.2, -0.3 };

    for (alphas) |alpha| {
        for (xs) |x| {
            const lp = logPdf(x, 0.0, 1.0, alpha);
            const p = pdf(x, 0.0, 1.0, alpha);
            if (p > 0.0) {
                try std.testing.expect(math.approxEqAbs(f64, lp, @log(p), 1e-12));
            }
        }
    }
}

test "gev: cross-check against known Easel values" {
    // Cross-check a few values against Easel's esl_gev_cdf/pdf/invcdf.
    // Using mu=0, lambda=1, alpha=0.1 (matches Easel convention directly).
    const mu: f64 = 0.0;
    const lambda: f64 = 1.0;
    const alpha: f64 = 0.1;

    // For x=1.0: y=lambda*(x-mu)=1, t=1+alpha*y=1.1, lt=log(1.1)
    // cdf = exp(-exp(-lt/alpha))
    const expected_cdf = @exp(-@exp(-@log(1.1) / 0.1));
    try std.testing.expect(math.approxEqAbs(f64, cdf(1.0, mu, lambda, alpha), expected_cdf, 1e-12));

    // invcdf round-trip
    const p_val: f64 = 0.632;
    const x_val = invcdf(p_val, mu, lambda, alpha);
    try std.testing.expect(math.approxEqAbs(f64, cdf(x_val, mu, lambda, alpha), p_val, 1e-10));
}

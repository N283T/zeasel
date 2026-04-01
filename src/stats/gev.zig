/// Generalized Extreme Value (GEV) distribution.
///
/// Parameters:
///   mu    - location parameter
///   sigma - scale parameter (must be > 0)
///   alpha - shape parameter (xi in some references)
///
/// The GEV unifies the three extreme value types:
///   alpha > 0: Type II (Frechet) - upper bound at mu - sigma/alpha
///   alpha < 0: Type III (Weibull) - lower bound at mu - sigma/alpha
///   alpha = 0: Type I (Gumbel) - no finite bound
///
/// Standardized variable: y = (x - mu) / sigma, t = 1 + alpha * y
/// When |alpha| < 1e-6, we use the Gumbel limit formulas.
const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;

/// Threshold below which |alpha| is treated as zero (Gumbel limit).
const alpha_threshold: f64 = 1e-6;

/// PDF of the GEV distribution.
///
/// When alpha != 0: (1/sigma) * t^(-1/alpha - 1) * exp(-t^(-1/alpha))
///   where t = 1 + alpha * (x - mu) / sigma, valid for t > 0.
/// When alpha == 0 (Gumbel): (1/sigma) * exp(-(y + exp(-y)))
pub fn pdf(x: f64, mu: f64, sigma: f64, alpha: f64) f64 {
    return @exp(logPdf(x, mu, sigma, alpha));
}

/// Log PDF of the GEV distribution.
pub fn logPdf(x: f64, mu: f64, sigma: f64, alpha: f64) f64 {
    const y = (x - mu) / sigma;

    if (@abs(alpha) < alpha_threshold) {
        // Gumbel limit: log(1/sigma) - y - exp(-y)
        return -@log(sigma) - y - @exp(-y);
    }

    const t = 1.0 + alpha * y;
    if (t <= 0.0) return -math.inf(f64);

    const inv_alpha = 1.0 / alpha;
    // log(1/sigma) + (-1/alpha - 1)*log(t) - t^(-1/alpha)
    return -@log(sigma) + (-inv_alpha - 1.0) * @log(t) - math.pow(f64, t, -inv_alpha);
}

/// CDF of the GEV distribution.
///
/// When alpha != 0: exp(-t^(-1/alpha)) where t = 1 + alpha * y, for t > 0.
/// When alpha == 0 (Gumbel): exp(-exp(-y))
pub fn cdf(x: f64, mu: f64, sigma: f64, alpha: f64) f64 {
    const y = (x - mu) / sigma;

    if (@abs(alpha) < alpha_threshold) {
        return @exp(-@exp(-y));
    }

    const t = 1.0 + alpha * y;
    if (t <= 0.0) {
        // Beyond the support boundary
        return if (alpha > 0.0) 0.0 else 1.0;
    }

    return @exp(-math.pow(f64, t, -1.0 / alpha));
}

/// Survival function: 1 - CDF.
pub fn surv(x: f64, mu: f64, sigma: f64, alpha: f64) f64 {
    const c = cdf(x, mu, sigma, alpha);
    // For numerical stability when CDF is near 1
    if (c > 0.9999999) {
        // Use log-based computation
        const lc = logCdf(x, mu, sigma, alpha);
        return -math.expm1(lc);
    }
    return 1.0 - c;
}

/// Log CDF of the GEV distribution.
fn logCdf(x: f64, mu: f64, sigma: f64, alpha: f64) f64 {
    const y = (x - mu) / sigma;

    if (@abs(alpha) < alpha_threshold) {
        return -@exp(-y);
    }

    const t = 1.0 + alpha * y;
    if (t <= 0.0) {
        return if (alpha > 0.0) -math.inf(f64) else 0.0;
    }

    return -math.pow(f64, t, -1.0 / alpha);
}

/// Inverse CDF (quantile function).
///
/// When alpha != 0: mu + (sigma/alpha) * ((-log(p))^(-alpha) - 1)
/// When alpha == 0: mu - sigma * log(-log(p))
pub fn invcdf(p: f64, mu: f64, sigma: f64, alpha: f64) f64 {
    if (p <= 0.0 or p >= 1.0) return math.nan(f64);

    const log_p = @log(p);

    if (@abs(alpha) < alpha_threshold) {
        // Gumbel: mu - sigma * log(-log(p))
        return mu - sigma * @log(-log_p);
    }

    // General: mu + (sigma/alpha) * ((-log(p))^(-alpha) - 1)
    return mu + (sigma / alpha) * (math.pow(f64, -log_p, -alpha) - 1.0);
}

/// Maximum likelihood estimation of GEV parameters from complete data.
///
/// Uses the method of probability-weighted moments (PWM, Hosking 1985) for
/// initial estimates, then refines with conjugate gradient minimization of
/// the negative log-likelihood.
///
/// Returns fitted (mu, sigma, alpha) parameters.
pub fn fitComplete(allocator: Allocator, x: []const f64) !struct { mu: f64, sigma: f64, alpha: f64 } {
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

    // sigma and mu from PWMs
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

    // Refine with conjugate gradient minimization.
    // Parameterize as p[0] = mu, p[1] = log(sigma), p[2] = alpha
    // to ensure sigma > 0.
    var p = [3]f64{ mu_init, @log(sigma_init), alpha_init };

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
        return .{ .mu = mu_init, .sigma = sigma_init, .alpha = alpha_init };
    };

    const fitted_sigma = @exp(p[1]);
    if (!math.isFinite(p[0]) or !math.isFinite(fitted_sigma) or !math.isFinite(p[2])) {
        // Fall back to PWM estimates
        return .{ .mu = mu_init, .sigma = sigma_init, .alpha = alpha_init };
    }

    return .{ .mu = p[0], .sigma = fitted_sigma, .alpha = p[2] };
}

const FitData = struct {
    x: []const f64,
};

/// Negative log-likelihood for GEV fitting.
/// Parameters: p[0] = mu, p[1] = log(sigma), p[2] = alpha
fn gevNegLogLik(p: []const f64, user_data: ?*anyopaque) f64 {
    const data: *const FitData = @ptrCast(@alignCast(user_data.?));
    const x = data.x;
    const mu = p[0];
    const sigma = @exp(p[1]);
    const alpha = p[2];

    var nll: f64 = 0.0;
    for (x) |xi| {
        const lp = logPdf(xi, mu, sigma, alpha);
        if (!math.isFinite(lp)) return math.inf(f64);
        nll -= lp;
    }
    return nll;
}

/// Gradient of the negative log-likelihood for GEV fitting.
/// Parameters: p[0] = mu, p[1] = w = log(sigma), p[2] = alpha
fn gevNegLogLikGrad(p: []const f64, dp: []f64, user_data: ?*anyopaque) void {
    const data: *const FitData = @ptrCast(@alignCast(user_data.?));
    const x = data.x;
    const mu = p[0];
    const sigma = @exp(p[1]);
    const alpha = p[2];
    const n: f64 = @floatFromInt(x.len);

    var d_mu: f64 = 0.0;
    var d_w: f64 = 0.0; // w.r.t. log(sigma)
    var d_alpha: f64 = 0.0;

    if (@abs(alpha) < alpha_threshold) {
        // Gumbel case gradients
        // logpdf = -log(sigma) - y - exp(-y), where y = (x - mu) / sigma
        // d/d_mu of logpdf = (1/sigma)(1 - exp(-y))
        // d/d_w of logpdf = -1 - y*(-1) + y*exp(-y)*(-1) ... need chain rule through w=log(sigma)
        for (x) |xi| {
            const y = (xi - mu) / sigma;
            const ey = @exp(-y);
            // d(logpdf)/d(mu) = (1 - ey) / sigma
            d_mu -= (1.0 - ey) / sigma;
            // d(logpdf)/d(w) = d(logpdf)/d(sigma) * sigma = -1 + y - y*ey ... via chain rule
            // logpdf = -w - y - exp(-y), y = (x-mu)*exp(-w)
            // d/dw = -1 + (x-mu)/sigma - (x-mu)/sigma * exp(-y) = -1 + y(1 - ey)
            d_w -= -1.0 + y * (1.0 - ey);
        }
        // d_alpha = 0 in Gumbel limit, but compute numerical approx for smooth transition
        dp[0] = d_mu;
        dp[1] = d_w;
        dp[2] = 0.0;
        return;
    }

    // General case: logpdf = -log(sigma) + (-1/alpha - 1)*log(t) - t^(-1/alpha)
    // where t = 1 + alpha * (x - mu) / sigma
    const inv_alpha = 1.0 / alpha;

    for (x) |xi| {
        const y = (xi - mu) / sigma;
        const t = 1.0 + alpha * y;

        if (t <= 0.0) {
            // Outside support; gradient is zero (boundary)
            continue;
        }

        const log_t = @log(t);
        const t_neg_inv_alpha = math.pow(f64, t, -inv_alpha);

        // dt/d(mu) = -alpha/sigma
        // dt/d(sigma) = -alpha*(x-mu)/sigma^2 = -alpha*y/sigma
        // dt/d(alpha) = y = (x-mu)/sigma

        // d(logpdf)/dt = (-1/alpha - 1)/t + (1/alpha) * t^(-1/alpha - 1)
        //              = (-1/alpha - 1)/t + t_neg_inv_alpha / (alpha * t)
        const dldt = ((-inv_alpha - 1.0) + t_neg_inv_alpha * inv_alpha) / t;

        // d(logpdf)/d(alpha) from the explicit alpha terms:
        // d/d(alpha) [(-1/alpha - 1)*log(t)] = (1/alpha^2)*log(t) + (-1/alpha-1)*(y/t)... no.
        // Let's be precise:
        // logpdf = -log(sigma) + (-1/alpha - 1)*log(t) - t^(-1/alpha)
        // d/d(alpha):
        //   (1/alpha^2)*log(t) + (-1/alpha - 1)*(y/t) - t^(-1/alpha) * (log(t)/alpha^2) + t^(-1/alpha)*(y/(alpha*t))
        //   = (1/alpha^2)*log(t)*(1 - t_neg_inv_alpha) + (y/t)*(-1/alpha - 1 + t_neg_inv_alpha/alpha)

        const da_direct = (1.0 / (alpha * alpha)) * log_t * (1.0 - t_neg_inv_alpha) +
            (y / t) * (-inv_alpha - 1.0 + t_neg_inv_alpha * inv_alpha);

        d_mu -= dldt * (-alpha / sigma);
        // d(logpdf)/d(w) = d(logpdf)/d(sigma) * sigma = dldt * (-alpha*y/sigma) * sigma + (-1)
        // = dldt * (-alpha*y) - 1
        d_w -= dldt * (-alpha * y) - 1.0;
        d_alpha -= da_direct;
    }

    // The -1 per sample for d_w from -log(sigma) term: already accumulated as -n in sum
    // Actually, let me reconsider. logpdf has -log(sigma) term.
    // d(-log(sigma))/d(w) = -1. We need to add -n to d_w from the nll.
    // Wait, we already subtracted (dldt * (-alpha*y) - 1.0), the -1.0 IS the -log(sigma) term.
    // Actually the sign: we do d_w -= (stuff), where stuff = d(logpdf)/d(w).
    // d(logpdf)/d(w) = -1 + dldt * dt/dw * sigma ... hmm let me redo.

    // logpdf = -w + (-1/alpha - 1)*log(t) - t^(-1/alpha), with w = log(sigma)
    // t = 1 + alpha*(x-mu)*exp(-w)
    // dt/dw = -alpha*(x-mu)*exp(-w) = -alpha*y (since y=(x-mu)/sigma=(x-mu)*exp(-w))
    // d(logpdf)/dw = -1 + dldt * (-alpha*y)
    // That's what we have. But we need to correct: the loop adds n * (-(-1)) = n to d_w.
    // Let me just zero-base and redo the w accumulation.

    // Actually, reviewing: d_w -= (dldt * (-alpha * y) - 1.0)
    // = d_w -= dldt*(-alpha*y) + 1.0
    // For negloglik, gradient = -sum(d(logpdf)/dw) = -sum(-1 + dldt*(-alpha*y))
    // = sum(1) - sum(dldt*(-alpha*y)) = n - sum(dldt*(-alpha*y))
    // Our code: d_w starts at 0, then d_w -= (-1 + dldt*(-alpha*y)) for each sample
    //         = d_w = sum(1 - dldt*(-alpha*y)) = n - sum(dldt*(-alpha*y))
    // That gives the negloglik gradient directly. But we need to negate once more?
    // No. d_mu starts 0, we do d_mu -= d(logpdf)/d(mu), so d_mu = -sum(d(logpdf)/d(mu)) = d(negloglik)/d(mu).
    // Same for d_w. So this is correct.

    _ = n;

    dp[0] = d_mu;
    dp[1] = d_w;
    dp[2] = d_alpha;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "gev: pdf/cdf/surv consistency" {
    // Test CDF + surv = 1.0 for various parameter combinations
    const params = [_]struct { mu: f64, sigma: f64, alpha: f64 }{
        .{ .mu = 0.0, .sigma = 1.0, .alpha = 0.0 }, // Gumbel
        .{ .mu = 1.0, .sigma = 2.0, .alpha = 0.1 }, // Frechet
        .{ .mu = 0.0, .sigma = 1.0, .alpha = -0.2 }, // Weibull-type
        .{ .mu = 5.0, .sigma = 0.5, .alpha = 0.3 },
    };

    for (params) |p| {
        const xs = [_]f64{ -2.0, 0.0, 1.0, 3.0, 5.0, 10.0 };
        for (xs) |x| {
            const c = cdf(x, p.mu, p.sigma, p.alpha);
            const s = surv(x, p.mu, p.sigma, p.alpha);
            try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-10));
        }
    }
}

test "gev: numerical derivative of CDF approximates PDF" {
    const mu: f64 = 0.0;
    const sigma: f64 = 1.0;
    const alpha: f64 = 0.1;
    const h: f64 = 1e-7;

    const xs = [_]f64{ -1.0, 0.0, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const numerical_pdf = (cdf(x + h, mu, sigma, alpha) - cdf(x - h, mu, sigma, alpha)) / (2.0 * h);
        const analytical_pdf = pdf(x, mu, sigma, alpha);
        try std.testing.expect(math.approxEqAbs(f64, numerical_pdf, analytical_pdf, 1e-5));
    }
}

test "gev: alpha=0 matches Gumbel" {
    // GEV with alpha=0 should match Gumbel distribution.
    // Gumbel uses lambda = 1/sigma parameterization.
    const gumbel = @import("gumbel.zig");
    const mu: f64 = 2.0;
    const sigma: f64 = 1.5;
    const lambda: f64 = 1.0 / sigma;

    const xs = [_]f64{ -2.0, 0.0, 1.0, 2.0, 3.0, 5.0, 10.0 };
    for (xs) |x| {
        const gev_cdf = cdf(x, mu, sigma, 0.0);
        const gumbel_cdf = gumbel.cdf(x, mu, lambda);
        try std.testing.expect(math.approxEqAbs(f64, gev_cdf, gumbel_cdf, 1e-10));

        const gev_pdf = pdf(x, mu, sigma, 0.0);
        const gumbel_pdf = gumbel.pdf(x, mu, lambda);
        try std.testing.expect(math.approxEqAbs(f64, gev_pdf, gumbel_pdf, 1e-10));
    }
}

test "gev: invcdf is inverse of cdf" {
    const params = [_]struct { mu: f64, sigma: f64, alpha: f64 }{
        .{ .mu = 0.0, .sigma = 1.0, .alpha = 0.0 },
        .{ .mu = 1.0, .sigma = 2.0, .alpha = 0.1 },
        .{ .mu = 0.0, .sigma = 1.0, .alpha = -0.2 },
    };

    const ps = [_]f64{ 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99 };
    for (params) |par| {
        for (ps) |p| {
            const x = invcdf(p, par.mu, par.sigma, par.alpha);
            const recovered_p = cdf(x, par.mu, par.sigma, par.alpha);
            try std.testing.expect(math.approxEqAbs(f64, recovered_p, p, 1e-10));
        }
    }
}

test "gev: cdf round-trip through invcdf" {
    const mu: f64 = 3.0;
    const sigma: f64 = 2.0;
    const alpha: f64 = 0.15;

    const xs = [_]f64{ 0.0, 2.0, 3.0, 5.0, 10.0 };
    for (xs) |x| {
        const p = cdf(x, mu, sigma, alpha);
        if (p > 0.0 and p < 1.0) {
            const x2 = invcdf(p, mu, sigma, alpha);
            try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-8));
        }
    }
}

test "gev: fitComplete recovers parameters" {
    // Generate samples from GEV(mu=5, sigma=2, alpha=0.1) using invcdf on LCG uniform.
    const true_mu: f64 = 5.0;
    const true_sigma: f64 = 2.0;
    const true_alpha: f64 = 0.1;
    const n = 5000;

    var samples: [n]f64 = undefined;
    var state: u64 = 42;
    for (&samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_sigma, true_alpha);
    }

    const fit = try fitComplete(std.testing.allocator, &samples);

    // With 5000 samples, expect reasonable accuracy
    try std.testing.expect(math.approxEqAbs(f64, fit.mu, true_mu, 0.3));
    try std.testing.expect(math.approxEqAbs(f64, fit.sigma, true_sigma, 0.3));
    try std.testing.expect(math.approxEqAbs(f64, fit.alpha, true_alpha, 0.25));
}

test "gev: fitComplete Gumbel case (alpha near 0)" {
    // Generate samples from GEV with alpha=0 (Gumbel)
    const true_mu: f64 = 10.0;
    const true_sigma: f64 = 3.0;
    const n = 5000;

    var samples: [n]f64 = undefined;
    var state: u64 = 12345;
    for (&samples) |*s| {
        state = (state *% 1664525 +% 1013904223) & 0xFFFFFFFF;
        const u = @as(f64, @floatFromInt(state + 1)) / @as(f64, 0x100000000);
        s.* = invcdf(u, true_mu, true_sigma, 0.0);
    }

    const fit = try fitComplete(std.testing.allocator, &samples);

    try std.testing.expect(math.approxEqAbs(f64, fit.mu, true_mu, 0.5));
    try std.testing.expect(math.approxEqAbs(f64, fit.sigma, true_sigma, 0.5));
    try std.testing.expect(@abs(fit.alpha) < 0.15);
}

test "gev: fitComplete error on insufficient data" {
    var x = [_]f64{ 1.0, 2.0 };
    try std.testing.expectError(error.InsufficientData, fitComplete(std.testing.allocator, &x));
}

test "gev: support boundaries" {
    // Type II (alpha > 0): support is (mu - sigma/alpha, +inf)
    // At the boundary, CDF should be 0
    const mu: f64 = 0.0;
    const sigma: f64 = 1.0;
    const alpha: f64 = 0.5;
    const boundary = mu - sigma / alpha; // = -2.0

    try std.testing.expectEqual(@as(f64, 0.0), cdf(boundary - 1.0, mu, sigma, alpha));
    try std.testing.expect(cdf(boundary + 0.1, mu, sigma, alpha) > 0.0);

    // Type III (alpha < 0): support is (-inf, mu - sigma/alpha)
    const alpha_neg: f64 = -0.5;
    const upper = mu - sigma / alpha_neg; // = 2.0

    try std.testing.expectEqual(@as(f64, 1.0), cdf(upper + 1.0, mu, sigma, alpha_neg));
    try std.testing.expect(cdf(upper - 0.1, mu, sigma, alpha_neg) < 1.0);
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

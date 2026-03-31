/// Root-finding routines: bisection and Newton-Raphson.
const std = @import("std");
const math = std.math;

pub const RootfinderError = error{
    DidNotConverge,
    BracketNotFound,
};

/// Find root of f(x) = 0 by bisection on [a, b].
/// Requires f(a) and f(b) to have opposite signs.
pub fn bisect(
    f: *const fn (f64) f64,
    a: f64,
    b: f64,
    tol: f64,
    max_iter: usize,
) !f64 {
    var lo = a;
    var hi = b;
    var f_lo = f(lo);

    if (f_lo * f(hi) > 0) return error.BracketNotFound;

    for (0..max_iter) |_| {
        const mid = (lo + hi) / 2.0;
        const f_mid = f(mid);
        if (@abs(f_mid) < tol or (hi - lo) / 2.0 < tol) return mid;
        if (f_lo * f_mid < 0) {
            hi = mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    return error.DidNotConverge;
}

/// Find root of f(x) = 0 by Newton-Raphson, starting at x0.
/// Requires f and its derivative df.
pub fn newton(
    f: *const fn (f64) f64,
    df: *const fn (f64) f64,
    x0: f64,
    tol: f64,
    max_iter: usize,
) !f64 {
    var x = x0;
    for (0..max_iter) |_| {
        const fx = f(x);
        if (@abs(fx) < tol) return x;
        const dfx = df(x);
        if (dfx == 0.0) return error.DidNotConverge;
        x -= fx / dfx;
    }
    return error.DidNotConverge;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

fn x_sq_minus_4(x: f64) f64 {
    return x * x - 4.0;
}

fn d_x_sq_minus_4(x: f64) f64 {
    return 2.0 * x;
}

test "bisect: x^2 - 4 on [0, 4] → root at 2.0" {
    const root = try bisect(x_sq_minus_4, 0.0, 4.0, 1e-10, 200);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), root, 1e-9);
}

test "bisect: wrong bracket returns BracketNotFound" {
    // f(1) = -3, f(3) = 5 — bracket is valid; flip to invalid: [3, 5] → f=5, f=21
    const result = bisect(x_sq_minus_4, 3.0, 5.0, 1e-10, 200);
    try std.testing.expectError(error.BracketNotFound, result);
}

test "newton: x^2 - 4, x0=3 → root at 2.0" {
    const root = try newton(x_sq_minus_4, d_x_sq_minus_4, 3.0, 1e-10, 100);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), root, 1e-9);
}

fn cos_x(x: f64) f64 {
    return math.cos(x);
}

fn neg_sin_x(x: f64) f64 {
    return -math.sin(x);
}

test "newton: cos(x) with x0=1.5 → root at π/2" {
    const root = try newton(cos_x, neg_sin_x, 1.5, 1e-10, 100);
    try std.testing.expectApproxEqAbs(math.pi / 2.0, root, 1e-9);
}

/// Conjugate gradient minimizer for unconstrained optimization (Polak-Ribière).
///
/// Uses bracket+Brent line search following Easel's esl_minimizer.c approach:
/// golden-section bracketing followed by Brent's method (parabolic interpolation
/// + golden-section bisection) for robust 1D minimization along search directions.
const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;

pub const MinimizerError = error{
    DidNotConverge,
    BracketingFailed,
    OutOfMemory,
};

/// Function to minimize: takes parameter vector x, returns f(x).
pub const ObjectiveFn = *const fn (x: []const f64, user_data: ?*anyopaque) f64;

/// Gradient function: takes x, writes gradient into grad.
pub const GradientFn = *const fn (x: []const f64, grad: []f64, user_data: ?*anyopaque) void;

pub const Options = struct {
    max_iter: usize = 200,
    tol: f64 = 1e-8,
    /// Relative tolerance for Brent's method. Should not be less than sqrt(DBL_EPSILON).
    brent_rtol: f64 = 1e-3,
    /// Absolute tolerance for Brent's method (avoids issues when minimum is at x=0).
    brent_atol: f64 = 1e-8,
    /// Maximum iterations for bracketing.
    brack_max_iter: usize = 100,
};

const golden_ratio: f64 = 1.61803398874989484820458683437;

/// Result from bracket(): a triplet (a, b, c) with f(b) < f(a) and f(b) < f(c).
const Bracket = struct {
    ax: f64,
    bx: f64,
    cx: f64,
    fa: f64,
    fb: f64,
    fc: f64,
};

/// Result from brent(): the scalar multiplier and function value at the minimum.
const BrentResult = struct {
    x: f64,
    fx: f64,
};

/// Bracket a minimum along a 1D line through `ori` in direction `dir`.
///
/// Finds scalars ax < bx < cx such that f(ori + bx*dir) < f(ori + ax*dir)
/// and f(ori + bx*dir) < f(ori + cx*dir), using golden-section expansion.
///
/// Based on Easel's bracket() in esl_minimizer.c.
fn bracket(
    ori: []const f64,
    dir: []const f64,
    firststep: f64,
    objective: ObjectiveFn,
    user_data: ?*anyopaque,
    wrk: []f64,
    max_iterations: usize,
) MinimizerError!Bracket {
    const n = ori.len;

    // Evaluate f(a) at origin (ax = 0)
    var ax: f64 = 0.0;
    var fa: f64 = objective(ori, user_data);

    // Evaluate f(b) at firststep
    var bx: f64 = firststep;
    for (0..n) |i| wrk[i] = ori[i] + bx * dir[i];
    var fb: f64 = objective(wrk, user_data);

    // If fb > fa, swap a and b so that fb <= fa.
    // This lets us search in the "wrong" direction when the initial step overshoots.
    if (fb > fa) {
        var swapper = ax;
        ax = bx;
        bx = swapper;
        swapper = fa;
        fa = fb;
        fb = swapper;
    }

    // First guess at c using golden ratio expansion
    var cx: f64 = bx + (bx - ax) * golden_ratio;
    for (0..n) |i| wrk[i] = ori[i] + cx * dir[i];
    var fc: f64 = objective(wrk, user_data);

    // Expand until fb < fc (i.e., b is below both a and c)
    var niter: usize = 0;
    while (fc <= fb) {
        // Slide: discard a, shift b->a, c->b, pick new c further out
        ax = bx;
        bx = cx;
        fa = fb;
        fb = fc;
        cx = bx + (bx - ax) * golden_ratio;
        for (0..n) |i| wrk[i] = ori[i] + cx * dir[i];
        fc = objective(wrk, user_data);

        // Rare: all three points have the same value (already at minimum)
        if (ax != bx and bx != cx and fa == fb and fb == fc) break;

        niter += 1;
        if (niter > max_iterations) return MinimizerError.BracketingFailed;
    }

    // Ensure a < b < c ordering for the caller
    if (ax > cx) {
        return .{
            .ax = cx,
            .bx = bx,
            .cx = ax,
            .fa = fc,
            .fb = fb,
            .fc = fa,
        };
    }
    return .{
        .ax = ax,
        .bx = bx,
        .cx = cx,
        .fa = fa,
        .fb = fb,
        .fc = fc,
    };
}

/// Find the minimum of f(ori + x*dir) within bracket [a, b] using Brent's method.
///
/// Combines golden-section search with parabolic interpolation for superlinear
/// convergence. Based on Easel's brent() in esl_minimizer.c, which follows
/// R.P. Brent (1973).
fn brent(
    ori: []const f64,
    dir: []const f64,
    objective: ObjectiveFn,
    user_data: ?*anyopaque,
    a_init: f64,
    b_init: f64,
    xvec: []f64,
    rtol: f64,
    atol: f64,
) BrentResult {
    const n = ori.len;
    const c = 1.0 - (1.0 / golden_ratio); // 0.381966... golden section ratio

    var a = a_init;
    var b = b_init;

    // Initial guess by golden section
    var x: f64 = a + c * (b - a);
    for (0..n) |i| xvec[i] = ori[i] + x * dir[i];
    var fx: f64 = objective(xvec, user_data);

    var v: f64 = x;
    var w: f64 = x;
    var fv: f64 = fx;
    var fw: f64 = fx;

    var d: f64 = 0.0;
    var e: f64 = 0.0; // distance moved on the step before last

    // Brent's algorithm is guaranteed to converge; no max_iter needed.
    while (true) {
        const m = 0.5 * (a + b);
        const tol = rtol * @abs(x) + atol;
        if (@abs(x - m) <= 2.0 * tol - 0.5 * (b - a)) break; // convergence

        var p: f64 = 0.0;
        var q: f64 = 0.0;
        var r: f64 = 0.0;
        if (@abs(e) > tol) {
            // Attempt parabolic interpolation
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0) {
                p = -p;
            } else {
                q = -q;
            }
            r = e;
            e = d;
        }

        if (@abs(p) < @abs(0.5 * q * r) or p < q * (a - x) or p < q * (b - x)) {
            // Parabolic interpolation step
            d = p / q;
            const u_candidate = x + d;
            if (2.0 * (u_candidate - a) < tol or 2.0 * (b - u_candidate) < tol) {
                d = if (x < m) tol else -tol;
            }
        } else {
            // Golden section step
            e = if (x < m) b - x else a - x;
            d = c * e;
        }

        // Evaluate f(), but not too close to x
        const u = if (@abs(d) >= tol)
            x + d
        else if (d > 0)
            x + tol
        else
            x - tol;

        for (0..n) |i| xvec[i] = ori[i] + u * dir[i];
        const fu = objective(xvec, user_data);

        // Bookkeeping: update bracket and best points
        if (fu <= fx) {
            if (u < x) {
                b = x;
            } else {
                a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if (u < x) {
                a = u;
            } else {
                b = u;
            }
            if (fu <= fw or w == x) {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if (fu <= fv or v == x or v == w) {
                v = u;
                fv = fu;
            }
        }
    }

    // Build final xvec at the minimum
    for (0..n) |i| xvec[i] = ori[i] + x * dir[i];
    return .{ .x = x, .fx = fx };
}

/// Minimize f(x) using conjugate gradient descent (Polak-Ribiere).
/// x is modified in place to the minimum. Returns the minimum function value.
///
/// Line search uses bracket+Brent: first brackets the minimum along the
/// conjugate direction using golden-section expansion, then refines with
/// Brent's method (parabolic interpolation + golden section).
pub fn minimize(
    allocator: Allocator,
    x: []f64,
    objective: ObjectiveFn,
    gradient: GradientFn,
    user_data: ?*anyopaque,
    options: Options,
) !f64 {
    const n = x.len;

    const grad = try allocator.alloc(f64, n);
    defer allocator.free(grad);
    const prev_grad = try allocator.alloc(f64, n);
    defer allocator.free(prev_grad);
    const direction = try allocator.alloc(f64, n);
    defer allocator.free(direction);
    const wrk = try allocator.alloc(f64, n);
    defer allocator.free(wrk);
    const x_new = try allocator.alloc(f64, n);
    defer allocator.free(x_new);

    gradient(x, grad, user_data);

    // Initial direction = -gradient (steepest descent)
    for (0..n) |i| direction[i] = -grad[i];

    var f_val = objective(x, user_data);

    for (0..options.max_iter) |_| {
        // Compute initial step size: min_i |1/direction[i]|
        // This ensures the first step doesn't overshoot in any dimension.
        var first_step: f64 = @abs(1.0 / direction[0]);
        for (1..n) |i| {
            const candidate = @abs(1.0 / direction[i]);
            if (candidate < first_step) first_step = candidate;
        }

        // Bracket the minimum along the search direction
        const brk = bracket(
            x,
            direction,
            first_step,
            objective,
            user_data,
            wrk,
            options.brack_max_iter,
        ) catch {
            // Bracketing failed; reset to steepest descent
            gradient(x, grad, user_data);
            for (0..n) |i| direction[i] = -grad[i];
            continue;
        };

        // Find precise minimum within bracket using Brent's method
        const result = brent(
            x,
            direction,
            objective,
            user_data,
            brk.ax,
            brk.cx,
            x_new,
            options.brent_rtol,
            options.brent_atol,
        );

        // If no improvement, we're done
        if (result.fx >= f_val) break;

        @memcpy(x, x_new);
        f_val = result.fx;

        // Save old gradient, compute new gradient at updated x
        @memcpy(prev_grad, grad);
        gradient(x, grad, user_data);

        // Check convergence: |grad|^2 < tol^2
        var grad_norm_sq: f64 = 0;
        for (grad) |g| grad_norm_sq += g * g;
        if (grad_norm_sq < options.tol * options.tol) return f_val;

        // Polak-Ribiere beta: beta = (g_new^T (g_new - g_old)) / |g_old|^2
        var dot_old: f64 = 0;
        var dot_diff: f64 = 0;
        for (0..n) |i| {
            dot_old += prev_grad[i] * prev_grad[i];
            dot_diff += grad[i] * (grad[i] - prev_grad[i]);
        }
        const beta = if (dot_old > 0) @max(0.0, dot_diff / dot_old) else 0.0;

        // Update conjugate direction
        for (0..n) |i| direction[i] = -grad[i] + beta * direction[i];
    }

    return f_val;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

fn quadratic_obj(x: []const f64, _: ?*anyopaque) f64 {
    const dx = x[0] - 3.0;
    return dx * dx;
}

fn quadratic_grad(x: []const f64, grad: []f64, _: ?*anyopaque) void {
    grad[0] = 2.0 * (x[0] - 3.0);
}

test "minimize 1-d quadratic (x-3)^2" {
    const allocator = std.testing.allocator;
    var x = [_]f64{0.0};
    const f_min = try minimize(allocator, &x, quadratic_obj, quadratic_grad, null, .{});
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), x[0], 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), f_min, 1e-12);
}

fn rosenbrock_obj(x: []const f64, _: ?*anyopaque) f64 {
    const a = 1.0 - x[0];
    const b = x[1] - x[0] * x[0];
    return a * a + 100.0 * b * b;
}

fn rosenbrock_grad(x: []const f64, grad: []f64, _: ?*anyopaque) void {
    grad[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
    grad[1] = 200.0 * (x[1] - x[0] * x[0]);
}

test "minimize Rosenbrock function" {
    const allocator = std.testing.allocator;
    var x = [_]f64{ 0.0, 0.0 };
    const f_min = try minimize(
        allocator,
        &x,
        rosenbrock_obj,
        rosenbrock_grad,
        null,
        .{ .max_iter = 10_000, .tol = 1e-6 },
    );
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), x[0], 1e-3);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), x[1], 1e-3);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), f_min, 1e-4);
}

/// Ill-conditioned elongated quadratic: f(x) = sum_i w_i * (x_i - c_i)^2
/// with condition number = max(w)/min(w).
/// The naive backtracking line search struggles with high condition numbers.
fn elongated_quadratic_obj(x: []const f64, _: ?*anyopaque) f64 {
    // Weights: 1.0, 1000.0 -- condition number 1000
    const weights = [_]f64{ 1.0, 1000.0 };
    const centers = [_]f64{ 5.0, -3.0 };
    var sum: f64 = 0.0;
    for (0..x.len) |i| {
        const d = x[i] - centers[i];
        sum += weights[i] * d * d;
    }
    return sum;
}

fn elongated_quadratic_grad(x: []const f64, grad: []f64, _: ?*anyopaque) void {
    const weights = [_]f64{ 1.0, 1000.0 };
    const centers = [_]f64{ 5.0, -3.0 };
    for (0..x.len) |i| {
        grad[i] = 2.0 * weights[i] * (x[i] - centers[i]);
    }
}

test "minimize ill-conditioned elongated quadratic" {
    const allocator = std.testing.allocator;
    var x = [_]f64{ 0.0, 0.0 };
    const f_min = try minimize(
        allocator,
        &x,
        elongated_quadratic_obj,
        elongated_quadratic_grad,
        null,
        .{ .max_iter = 500, .tol = 1e-8 },
    );
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), x[0], 1e-4);
    try std.testing.expectApproxEqAbs(@as(f64, -3.0), x[1], 1e-4);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), f_min, 1e-4);
}

/// Higher-dimensional ill-conditioned problem: 4D with condition number 1000.
fn high_dim_ill_conditioned_obj(x: []const f64, _: ?*anyopaque) f64 {
    var sum: f64 = 0.0;
    for (0..x.len) |i| {
        // Weight grows exponentially: 1, 10, 100, 1000
        const fi: f64 = @floatFromInt(i);
        const w = math.pow(f64, 10.0, fi);
        const d = x[i] - @as(f64, @floatFromInt(i + 1));
        sum += w * d * d;
    }
    return sum;
}

fn high_dim_ill_conditioned_grad(x: []const f64, grad: []f64, _: ?*anyopaque) void {
    for (0..x.len) |i| {
        const fi: f64 = @floatFromInt(i);
        const w = math.pow(f64, 10.0, fi);
        grad[i] = 2.0 * w * (x[i] - @as(f64, @floatFromInt(i + 1)));
    }
}

test "minimize high-dimensional ill-conditioned quadratic" {
    const allocator = std.testing.allocator;
    var x = [_]f64{ 0.0, 0.0, 0.0, 0.0 };
    const f_min = try minimize(
        allocator,
        &x,
        high_dim_ill_conditioned_obj,
        high_dim_ill_conditioned_grad,
        null,
        .{ .max_iter = 1000, .tol = 1e-6 },
    );
    for (0..x.len) |i| {
        const expected: f64 = @floatFromInt(i + 1);
        try std.testing.expectApproxEqAbs(expected, x[i], 1e-3);
    }
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), f_min, 1e-2);
}

test "bracket finds valid triplet" {
    const allocator = std.testing.allocator;
    const wrk = try allocator.alloc(f64, 1);
    defer allocator.free(wrk);

    const ori = [_]f64{0.0};
    const dir = [_]f64{1.0};

    const brk = try bracket(
        &ori,
        &dir,
        1.0,
        quadratic_obj,
        null,
        wrk,
        100,
    );
    // b should be between a and c, and f(b) <= f(a) and f(b) <= f(c)
    try std.testing.expect(brk.ax <= brk.bx);
    try std.testing.expect(brk.bx <= brk.cx);
    try std.testing.expect(brk.fb <= brk.fa);
    try std.testing.expect(brk.fb <= brk.fc);
}

test "brent finds minimum of 1D quadratic" {
    const allocator = std.testing.allocator;
    const xvec = try allocator.alloc(f64, 1);
    defer allocator.free(xvec);

    const ori = [_]f64{0.0};
    const dir = [_]f64{1.0};

    // Bracket [0, 6] contains the minimum of (x-3)^2 at x=3
    const result = brent(
        &ori,
        &dir,
        quadratic_obj,
        null,
        0.0,
        6.0,
        xvec,
        1e-3,
        1e-8,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), result.x, 1e-3);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.fx, 1e-4);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), xvec[0], 1e-3);
}

/// Conjugate gradient minimizer for unconstrained optimization (Polak-Ribière).
const std = @import("std");
const Allocator = std.mem.Allocator;

pub const MinimizerError = error{
    DidNotConverge,
    OutOfMemory,
};

/// Function to minimize: takes parameter vector x, returns f(x).
pub const ObjectiveFn = *const fn (x: []const f64, user_data: ?*anyopaque) f64;

/// Gradient function: takes x, writes gradient into grad.
pub const GradientFn = *const fn (x: []const f64, grad: []f64, user_data: ?*anyopaque) void;

pub const Options = struct {
    max_iter: usize = 200,
    tol: f64 = 1e-8,
};

/// Minimize f(x) using conjugate gradient descent (Polak-Ribière).
/// x is modified in place to the minimum. Returns the minimum function value.
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
    const x_new = try allocator.alloc(f64, n);
    defer allocator.free(x_new);

    gradient(x, grad, user_data);

    // Initial direction = -gradient (steepest descent)
    for (0..n) |i| direction[i] = -grad[i];

    var f_val = objective(x, user_data);

    for (0..options.max_iter) |_| {
        // Backtracking line search along direction
        var step: f64 = 1.0;
        var f_new: f64 = undefined;
        var found = false;
        for (0..60) |_| {
            for (0..n) |i| x_new[i] = x[i] + step * direction[i];
            f_new = objective(x_new, user_data);
            if (f_new < f_val) {
                found = true;
                break;
            }
            step *= 0.5;
        }
        if (!found) {
            // Reset to steepest descent and try again with a tiny step
            for (0..n) |i| direction[i] = -grad[i];
            step = 1e-6;
            for (0..n) |i| x_new[i] = x[i] + step * direction[i];
            f_new = objective(x_new, user_data);
            if (f_new >= f_val) break; // Truly stuck
        }

        @memcpy(x, x_new);
        f_val = f_new;

        // Save old gradient, compute new gradient at updated x
        @memcpy(prev_grad, grad);
        gradient(x, grad, user_data);

        // Check convergence: |grad|^2 < tol^2
        var grad_norm_sq: f64 = 0;
        for (grad) |g| grad_norm_sq += g * g;
        if (grad_norm_sq < options.tol * options.tol) return f_val;

        // Polak-Ribière beta: beta = (g_new^T (g_new - g_old)) / |g_old|^2
        var dot_new: f64 = 0;
        var dot_diff: f64 = 0;
        for (0..n) |i| {
            dot_new += grad[i] * grad[i];
            dot_diff += grad[i] * (grad[i] - prev_grad[i]);
        }
        const beta = @max(0.0, dot_diff / dot_new); // Clamp to 0 to reset on bad directions

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

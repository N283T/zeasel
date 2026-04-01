const std = @import("std");
const math = std.math;

/// Returns a namespace with vector operations for the given element type T.
/// Reduction operations (sum, dot, max, min) use @Vector for SIMD acceleration.
/// Float-only operations (normalize, entropy, logSum) are gated with comptime checks.
pub fn VectorOps(comptime T: type) type {
    return struct {
        // ----------------------------------------------------------------
        // Compile-time helpers
        // ----------------------------------------------------------------

        const is_float = @typeInfo(T) == .float;
        const vec_len = std.simd.suggestVectorLength(T) orelse 4;
        const V = @Vector(vec_len, T);

        // ----------------------------------------------------------------
        // Element-wise operations
        // ----------------------------------------------------------------

        /// Set all elements of dst to value.
        pub fn set(dst: []T, value: T) void {
            const splat: V = @splat(value);
            var i: usize = 0;
            while (i + vec_len <= dst.len) : (i += vec_len) {
                dst[i..][0..vec_len].* = splat;
            }
            while (i < dst.len) : (i += 1) {
                dst[i] = value;
            }
        }

        /// Scale all elements of dst by factor (in-place).
        pub fn scale(dst: []T, factor: T) void {
            const splat: V = @splat(factor);
            var i: usize = 0;
            while (i + vec_len <= dst.len) : (i += vec_len) {
                const chunk: V = dst[i..][0..vec_len].*;
                dst[i..][0..vec_len].* = chunk * splat;
            }
            while (i < dst.len) : (i += 1) {
                dst[i] *= factor;
            }
        }

        /// Element-wise add: dst[i] += src[i]. Slices must have equal length.
        pub fn add(dst: []T, src: []const T) void {
            std.debug.assert(dst.len == src.len);
            var i: usize = 0;
            while (i + vec_len <= dst.len) : (i += vec_len) {
                const d: V = dst[i..][0..vec_len].*;
                const s: V = src[i..][0..vec_len].*;
                dst[i..][0..vec_len].* = d + s;
            }
            while (i < dst.len) : (i += 1) {
                dst[i] += src[i];
            }
        }

        /// Scaled add: dst[i] += a * src[i]. Slices must have equal length.
        pub fn addScaled(dst: []T, src: []const T, a: T) void {
            std.debug.assert(dst.len == src.len);
            const splat: V = @splat(a);
            var i: usize = 0;
            while (i + vec_len <= dst.len) : (i += vec_len) {
                const d: V = dst[i..][0..vec_len].*;
                const s: V = src[i..][0..vec_len].*;
                dst[i..][0..vec_len].* = d + s * splat;
            }
            while (i < dst.len) : (i += 1) {
                dst[i] += a * src[i];
            }
        }

        // ----------------------------------------------------------------
        // Reductions
        // ----------------------------------------------------------------

        /// Sum all elements of v.
        ///
        /// Note: For small integer types (u8, i8), the result type is the same as the input type,
        /// so overflow may occur for large slices. Use a wider type (e.g., VectorOps(u32)) for safe accumulation.
        pub fn sum(v: []const T) T {
            var acc: V = @splat(0);
            var i: usize = 0;
            while (i + vec_len <= v.len) : (i += vec_len) {
                const chunk: V = v[i..][0..vec_len].*;
                acc += chunk;
            }
            var result = @reduce(.Add, acc);
            while (i < v.len) : (i += 1) {
                result += v[i];
            }
            return result;
        }

        /// Dot product of two equal-length slices.
        pub fn dot(a: []const T, b: []const T) T {
            std.debug.assert(a.len == b.len);
            var acc: V = @splat(0);
            var i: usize = 0;
            while (i + vec_len <= a.len) : (i += vec_len) {
                const va: V = a[i..][0..vec_len].*;
                const vb: V = b[i..][0..vec_len].*;
                acc += va * vb;
            }
            var result = @reduce(.Add, acc);
            while (i < a.len) : (i += 1) {
                result += a[i] * b[i];
            }
            return result;
        }

        /// Return the maximum element. Panics on empty slice.
        pub fn max(v: []const T) T {
            std.debug.assert(v.len > 0);
            var acc: V = @splat(v[0]);
            var i: usize = 0;
            while (i + vec_len <= v.len) : (i += vec_len) {
                const chunk: V = v[i..][0..vec_len].*;
                acc = @max(acc, chunk);
            }
            var result = @reduce(.Max, acc);
            while (i < v.len) : (i += 1) {
                if (v[i] > result) result = v[i];
            }
            return result;
        }

        /// Return the minimum element. Panics on empty slice.
        pub fn min(v: []const T) T {
            std.debug.assert(v.len > 0);
            var acc: V = @splat(v[0]);
            var i: usize = 0;
            while (i + vec_len <= v.len) : (i += vec_len) {
                const chunk: V = v[i..][0..vec_len].*;
                acc = @min(acc, chunk);
            }
            var result = @reduce(.Min, acc);
            while (i < v.len) : (i += 1) {
                if (v[i] < result) result = v[i];
            }
            return result;
        }

        /// Return the index of the first maximum element. Panics on empty slice.
        pub fn argmax(v: []const T) usize {
            std.debug.assert(v.len > 0);
            var best_val = v[0];
            var best_idx: usize = 0;
            for (v[1..], 1..) |x, i| {
                if (x > best_val) {
                    best_val = x;
                    best_idx = i;
                }
            }
            return best_idx;
        }

        /// Return the index of the first minimum element. Panics on empty slice.
        pub fn argmin(v: []const T) usize {
            std.debug.assert(v.len > 0);
            var best_val = v[0];
            var best_idx: usize = 0;
            for (v[1..], 1..) |x, i| {
                if (x < best_val) {
                    best_val = x;
                    best_idx = i;
                }
            }
            return best_idx;
        }

        // ----------------------------------------------------------------
        // Utilities
        // ----------------------------------------------------------------

        /// Copy src into dst. Slices must have equal length.
        pub fn copy(dst: []T, src: []const T) void {
            std.debug.assert(dst.len == src.len);
            @memcpy(dst, src);
        }

        /// Reverse v in-place.
        pub fn reverse(v: []T) void {
            var lo: usize = 0;
            var hi: usize = v.len;
            while (lo < hi) {
                hi -= 1;
                const tmp = v[lo];
                v[lo] = v[hi];
                v[hi] = tmp;
                lo += 1;
            }
        }

        // ----------------------------------------------------------------
        // Float-only operations
        // ----------------------------------------------------------------

        /// Normalize v so elements sum to 1 (L1 norm).
        /// When sum is zero, sets all elements to 1/n (uniform distribution),
        /// matching Easel's esl_vec_DNorm behavior.
        /// Only available for float types.
        pub fn normalize(v: []T) void {
            comptime if (!is_float) @compileError("normalize requires a float type");
            const s = sum(v);
            if (s == 0) {
                if (v.len > 0) {
                    const uniform: T = 1.0 / @as(T, @floatFromInt(v.len));
                    for (v) |*x| x.* = uniform;
                }
                return;
            }
            scale(v, 1 / s);
        }

        /// Compute Shannon entropy: -sum(p[i] * log2(p[i])) for p[i] > 0.
        /// Only available for float types.
        pub fn entropy(p: []const T) T {
            comptime if (!is_float) @compileError("entropy requires a float type");
            var h: T = 0;
            for (p) |x| {
                if (x > 0) h -= x * @log2(x);
            }
            return h;
        }

        /// Numerically stable log(sum(exp(v))): find max, then max + log(sum(exp(v[i]-max))).
        /// Guards against +inf max and skips extreme underflow terms.
        /// Reference: Easel esl_vec_DLogSum.
        /// Only available for float types.
        pub fn logSum(v: []const T) T {
            comptime if (!is_float) @compileError("logSum requires a float type");
            if (v.len == 0) return -math.inf(T);
            const m = max(v);
            if (math.isPositiveInf(m)) return math.inf(T);
            var s: T = 0;
            for (v) |x| {
                // Skip terms that would cause extreme underflow (or NaN from -inf - -inf).
                if (x >= m - 500) {
                    s += @exp(x - m);
                }
            }
            return m + @log(s);
        }

        /// Convert a log-probability vector to a normalized probability vector.
        /// Subtracts logSum, exponentiates, then re-normalizes.
        /// Reference: Easel esl_vec_DLogNorm.
        /// Only available for float types.
        pub fn logNorm(v: []T) void {
            comptime if (!is_float) @compileError("logNorm requires a float type");
            if (v.len == 0) return;
            const denom = logSum(v);
            for (v) |*x| {
                x.* = @exp(x.* - denom);
            }
            normalize(v);
        }

        /// Element-wise natural log: dest[i] = ln(src[i]).
        /// For src[i] <= 0, dest[i] = -inf.
        /// Only available for float types. Slices must have equal length.
        pub fn log(dest: []T, src: []const T) void {
            comptime if (!is_float) @compileError("log requires a float type");
            std.debug.assert(dest.len == src.len);
            for (dest, src) |*d, s| {
                d.* = if (s > 0) @log(s) else -math.inf(T);
            }
        }

        /// Element-wise exp: dest[i] = exp(src[i]).
        /// Only available for float types. Slices must have equal length.
        pub fn exp(dest: []T, src: []const T) void {
            comptime if (!is_float) @compileError("exp requires a float type");
            std.debug.assert(dest.len == src.len);
            for (dest, src) |*d, s| {
                d.* = @exp(s);
            }
        }

        /// Kullback-Leibler divergence D_KL(p||q) = sum(p[i] * log2(p[i]/q[i])) in bits.
        /// Terms where p[i] == 0 are skipped. Returns inf if any p[i] > 0 and q[i] == 0.
        /// Reference: Easel esl_vec_DRelEntropy (returns bits, not nats).
        /// Only available for float types. Slices must have equal length.
        pub fn relativeEntropy(p: []const T, q: []const T) T {
            comptime if (!is_float) @compileError("relativeEntropy requires a float type");
            std.debug.assert(p.len == q.len);
            var kl: T = 0;
            for (p, q) |pi, qi| {
                if (pi > 0) {
                    if (qi == 0) return math.inf(T);
                    kl += pi * @log2(pi / qi);
                }
            }
            return kl;
        }

        /// Sort elements in ascending order (in-place).
        pub fn sortIncreasing(v: []T) void {
            std.mem.sort(T, v, {}, std.sort.asc(T));
        }

        /// Sort elements in descending order (in-place).
        pub fn sortDecreasing(v: []T) void {
            std.mem.sort(T, v, {}, std.sort.desc(T));
        }

        /// Return true if all elements are finite (no NaN, no Inf).
        /// Only available for float types.
        pub fn validate(v: []const T) bool {
            comptime if (!is_float) @compileError("validate requires a float type");
            for (v) |x| {
                if (math.isNan(x) or math.isInf(x)) return false;
            }
            return true;
        }

        /// Return true if all elements are valid log-probabilities (<= 0, not NaN, not +inf).
        /// -inf is accepted (represents log(0) = probability 0).
        /// Only available for float types.
        pub fn logValidate(v: []const T) bool {
            comptime if (!is_float) @compileError("logValidate requires a float type");
            for (v) |x| {
                if (math.isNan(x) or math.isPositiveInf(x) or x > 0) return false;
            }
            return true;
        }
    };
}

/// Pre-instantiated convenience aliases.
pub const VecF64 = VectorOps(f64);
pub const VecF32 = VectorOps(f32);
pub const VecI32 = VectorOps(i32);
pub const VecU8 = VectorOps(u8);

// ============================================================================
// Tests
// ============================================================================

const testing = std.testing;

// ---- f64 tests ----

test "VecF64.set" {
    var buf: [8]f64 = undefined;
    VecF64.set(&buf, 3.14);
    for (buf) |x| try testing.expectApproxEqAbs(3.14, x, 1e-12);
}

test "VecF64.scale" {
    var buf = [_]f64{ 1, 2, 3, 4 };
    VecF64.scale(&buf, 2.0);
    try testing.expectApproxEqAbs(2.0, buf[0], 1e-12);
    try testing.expectApproxEqAbs(4.0, buf[1], 1e-12);
    try testing.expectApproxEqAbs(6.0, buf[2], 1e-12);
    try testing.expectApproxEqAbs(8.0, buf[3], 1e-12);
}

test "VecF64.add" {
    var dst = [_]f64{ 1, 2, 3, 4 };
    const src = [_]f64{ 10, 20, 30, 40 };
    VecF64.add(&dst, &src);
    try testing.expectApproxEqAbs(11.0, dst[0], 1e-12);
    try testing.expectApproxEqAbs(22.0, dst[1], 1e-12);
    try testing.expectApproxEqAbs(33.0, dst[2], 1e-12);
    try testing.expectApproxEqAbs(44.0, dst[3], 1e-12);
}

test "VecF64.addScaled" {
    var dst = [_]f64{ 0, 0, 0, 0 };
    const src = [_]f64{ 1, 2, 3, 4 };
    VecF64.addScaled(&dst, &src, 2.0);
    try testing.expectApproxEqAbs(2.0, dst[0], 1e-12);
    try testing.expectApproxEqAbs(4.0, dst[1], 1e-12);
    try testing.expectApproxEqAbs(6.0, dst[2], 1e-12);
    try testing.expectApproxEqAbs(8.0, dst[3], 1e-12);
}

test "VecF64.sum small" {
    const v = [_]f64{ 1, 2, 3, 4 };
    try testing.expectApproxEqAbs(10.0, VecF64.sum(&v), 1e-12);
}

test "VecF64.sum large (SIMD path)" {
    var v: [100]f64 = undefined;
    for (&v, 0..) |*x, i| x.* = @floatFromInt(i + 1);
    // sum 1..100 = 5050
    try testing.expectApproxEqAbs(5050.0, VecF64.sum(&v), 1e-9);
}

test "VecF64.dot" {
    const a = [_]f64{ 1, 2, 3 };
    const b = [_]f64{ 4, 5, 6 };
    try testing.expectApproxEqAbs(32.0, VecF64.dot(&a, &b), 1e-12);
}

test "VecF64.max and min" {
    const v = [_]f64{ 3, 1, 4, 1, 5, 9, 2, 6 };
    try testing.expectApproxEqAbs(9.0, VecF64.max(&v), 1e-12);
    try testing.expectApproxEqAbs(1.0, VecF64.min(&v), 1e-12);
}

test "VecF64.argmax and argmin" {
    const v = [_]f64{ 3, 1, 4, 1, 5, 9, 2, 6 };
    try testing.expectEqual(@as(usize, 5), VecF64.argmax(&v));
    try testing.expectEqual(@as(usize, 1), VecF64.argmin(&v));
}

test "VecF64.normalize" {
    var v = [_]f64{ 1, 2, 3, 4 };
    VecF64.normalize(&v);
    try testing.expectApproxEqAbs(0.1, v[0], 1e-12);
    try testing.expectApproxEqAbs(0.2, v[1], 1e-12);
    try testing.expectApproxEqAbs(0.3, v[2], 1e-12);
    try testing.expectApproxEqAbs(0.4, v[3], 1e-12);
}

test "VecF64.entropy uniform 4" {
    // Uniform distribution over 4 states has entropy 2.0 bits.
    var v = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    try testing.expectApproxEqAbs(2.0, VecF64.entropy(&v), 1e-12);
}

test "VecF64.logSum" {
    // log(exp(1) + exp(2)) = log(e + e^2)
    const v = [_]f64{ 1.0, 2.0 };
    const expected = @log(@exp(1.0) + @exp(2.0));
    try testing.expectApproxEqAbs(expected, VecF64.logSum(&v), 1e-12);
}

test "VecF64.copy independence" {
    const src = [_]f64{ 1, 2, 3, 4 };
    var dst: [4]f64 = undefined;
    VecF64.copy(&dst, &src);
    // Mutate dst, src must be unchanged.
    dst[0] = 99;
    try testing.expectApproxEqAbs(1.0, src[0], 1e-12);
    try testing.expectApproxEqAbs(99.0, dst[0], 1e-12);
}

test "VecF64.reverse" {
    var v = [_]f64{ 1, 2, 3, 4 };
    VecF64.reverse(&v);
    try testing.expectApproxEqAbs(4.0, v[0], 1e-12);
    try testing.expectApproxEqAbs(3.0, v[1], 1e-12);
    try testing.expectApproxEqAbs(2.0, v[2], 1e-12);
    try testing.expectApproxEqAbs(1.0, v[3], 1e-12);
}

// ---- i32 tests ----

test "VecI32.sum" {
    const v = [_]i32{ 1, 2, 3, 4 };
    try testing.expectEqual(@as(i32, 10), VecI32.sum(&v));
}

test "VecI32.max and min" {
    const v = [_]i32{ -5, 3, 0, 7, -1 };
    try testing.expectEqual(@as(i32, 7), VecI32.max(&v));
    try testing.expectEqual(@as(i32, -5), VecI32.min(&v));
}

test "VecI32.argmax" {
    const v = [_]i32{ 1, 8, 3, 2 };
    try testing.expectEqual(@as(usize, 1), VecI32.argmax(&v));
}

// ---- u8 tests ----

test "VecU8.sum" {
    const v = [_]u8{ 10, 20, 30, 40 };
    try testing.expectEqual(@as(u8, 100), VecU8.sum(&v));
}

test "VecU8.max" {
    const v = [_]u8{ 5, 200, 13, 99 };
    try testing.expectEqual(@as(u8, 200), VecU8.max(&v));
}

// ---- Large array SIMD correctness ----

test "VecF64 large array sum matches scalar" {
    var buf: [1024]f64 = undefined;
    for (&buf, 0..) |*x, i| x.* = @floatFromInt(i);
    // Compare SIMD result with naive scalar sum.
    var scalar_sum: f64 = 0;
    for (buf) |x| scalar_sum += x;
    const simd_sum = VecF64.sum(&buf);
    try testing.expectApproxEqAbs(scalar_sum, simd_sum, 1e-6);
}

test "VecI32 large array sum matches scalar" {
    var buf: [1000]i32 = undefined;
    for (&buf, 0..) |*x, i| x.* = @intCast(i % 100);
    var scalar_sum: i32 = 0;
    for (buf) |x| scalar_sum += x;
    try testing.expectEqual(scalar_sum, VecI32.sum(&buf));
}

// ---- New f64 tests: logNorm, log, exp, relativeEntropy, sort, validate ----

test "VecF64.logNorm" {
    // Start with log-probabilities that don't sum to 1 in probability space.
    // After logNorm, sum(exp(v[i])) should equal 1.
    var v = [_]f64{ -1.0, -2.0, -3.0 };
    VecF64.logNorm(&v);
    // Verify sum of exp(v[i]) == 1.
    var prob_sum: f64 = 0;
    for (v) |x| prob_sum += @exp(x);
    try testing.expectApproxEqAbs(1.0, prob_sum, 1e-12);
}

test "VecF64.logNorm preserves relative order" {
    var v = [_]f64{ 0.0, 1.0, 2.0, 3.0 };
    VecF64.logNorm(&v);
    // Values should still be increasing.
    try testing.expect(v[0] < v[1]);
    try testing.expect(v[1] < v[2]);
    try testing.expect(v[2] < v[3]);
}

test "VecF64.logNorm empty" {
    var v = [_]f64{};
    VecF64.logNorm(&v); // should not panic
}

test "VecF64.log" {
    const src = [_]f64{ 1.0, @exp(1.0), @exp(2.0) };
    var dest: [3]f64 = undefined;
    VecF64.log(&dest, &src);
    try testing.expectApproxEqAbs(0.0, dest[0], 1e-12);
    try testing.expectApproxEqAbs(1.0, dest[1], 1e-12);
    try testing.expectApproxEqAbs(2.0, dest[2], 1e-12);
}

test "VecF64.log handles zero" {
    const src = [_]f64{ 0.0, 1.0 };
    var dest: [2]f64 = undefined;
    VecF64.log(&dest, &src);
    try testing.expect(dest[0] == -std.math.inf(f64));
    try testing.expectApproxEqAbs(0.0, dest[1], 1e-12);
}

test "VecF64.exp" {
    const src = [_]f64{ 0.0, 1.0, 2.0 };
    var dest: [3]f64 = undefined;
    VecF64.exp(&dest, &src);
    try testing.expectApproxEqAbs(1.0, dest[0], 1e-12);
    try testing.expectApproxEqAbs(@exp(1.0), dest[1], 1e-12);
    try testing.expectApproxEqAbs(@exp(2.0), dest[2], 1e-12);
}

test "VecF64.log and exp roundtrip" {
    const original = [_]f64{ 0.25, 0.5, 1.0, 2.0, 10.0 };
    var logged: [5]f64 = undefined;
    var restored: [5]f64 = undefined;
    VecF64.log(&logged, &original);
    VecF64.exp(&restored, &logged);
    for (original, restored) |o, r| {
        try testing.expectApproxEqAbs(o, r, 1e-12);
    }
}

test "VecF64.relativeEntropy identical distributions" {
    const p = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    const q = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    try testing.expectApproxEqAbs(0.0, VecF64.relativeEntropy(&p, &q), 1e-12);
}

test "VecF64.relativeEntropy known value" {
    // D_KL([0.5, 0.5, 0, 0] || [0.25, 0.25, 0.25, 0.25])
    // = 0.5*ln(0.5/0.25) + 0.5*ln(0.5/0.25) = ln(2)
    const p = [_]f64{ 0.5, 0.5, 0, 0 };
    const q = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    try testing.expectApproxEqAbs(@log(2.0), VecF64.relativeEntropy(&p, &q), 1e-12);
}

test "VecF64.relativeEntropy q_i zero with p_i positive returns inf" {
    const p = [_]f64{ 0.5, 0.5 };
    const q = [_]f64{ 1.0, 0.0 };
    try testing.expect(VecF64.relativeEntropy(&p, &q) == std.math.inf(f64));
}

test "VecF64.sortIncreasing" {
    var v = [_]f64{ 3, 1, 4, 1, 5, 9, 2, 6 };
    VecF64.sortIncreasing(&v);
    try testing.expectApproxEqAbs(1.0, v[0], 1e-12);
    try testing.expectApproxEqAbs(1.0, v[1], 1e-12);
    try testing.expectApproxEqAbs(2.0, v[2], 1e-12);
    try testing.expectApproxEqAbs(3.0, v[3], 1e-12);
    try testing.expectApproxEqAbs(4.0, v[4], 1e-12);
    try testing.expectApproxEqAbs(5.0, v[5], 1e-12);
    try testing.expectApproxEqAbs(6.0, v[6], 1e-12);
    try testing.expectApproxEqAbs(9.0, v[7], 1e-12);
}

test "VecF64.sortDecreasing" {
    var v = [_]f64{ 3, 1, 4, 1, 5, 9, 2, 6 };
    VecF64.sortDecreasing(&v);
    try testing.expectApproxEqAbs(9.0, v[0], 1e-12);
    try testing.expectApproxEqAbs(6.0, v[1], 1e-12);
    try testing.expectApproxEqAbs(5.0, v[2], 1e-12);
    try testing.expectApproxEqAbs(4.0, v[3], 1e-12);
    try testing.expectApproxEqAbs(3.0, v[4], 1e-12);
    try testing.expectApproxEqAbs(2.0, v[5], 1e-12);
    try testing.expectApproxEqAbs(1.0, v[6], 1e-12);
    try testing.expectApproxEqAbs(1.0, v[7], 1e-12);
}

test "VecI32.sortIncreasing" {
    var v = [_]i32{ 5, 3, 8, 1, 4 };
    VecI32.sortIncreasing(&v);
    try testing.expectEqual(@as(i32, 1), v[0]);
    try testing.expectEqual(@as(i32, 3), v[1]);
    try testing.expectEqual(@as(i32, 4), v[2]);
    try testing.expectEqual(@as(i32, 5), v[3]);
    try testing.expectEqual(@as(i32, 8), v[4]);
}

test "VecI32.sortDecreasing" {
    var v = [_]i32{ 5, 3, 8, 1, 4 };
    VecI32.sortDecreasing(&v);
    try testing.expectEqual(@as(i32, 8), v[0]);
    try testing.expectEqual(@as(i32, 5), v[1]);
    try testing.expectEqual(@as(i32, 4), v[2]);
    try testing.expectEqual(@as(i32, 3), v[3]);
    try testing.expectEqual(@as(i32, 1), v[4]);
}

test "VecF64.validate all finite" {
    const v = [_]f64{ 1.0, -2.0, 0.0, 3.14 };
    try testing.expect(VecF64.validate(&v));
}

test "VecF64.validate with NaN" {
    const v = [_]f64{ 1.0, std.math.nan(f64), 3.0 };
    try testing.expect(!VecF64.validate(&v));
}

test "VecF64.validate with Inf" {
    const v = [_]f64{ 1.0, std.math.inf(f64), 3.0 };
    try testing.expect(!VecF64.validate(&v));
}

test "VecF64.validate with negative Inf" {
    const v = [_]f64{ 1.0, -std.math.inf(f64), 3.0 };
    try testing.expect(!VecF64.validate(&v));
}

test "VecF64.validate empty" {
    const v = [_]f64{};
    try testing.expect(VecF64.validate(&v));
}

test "VecF64.logValidate valid log-probabilities" {
    const v = [_]f64{ -1.0, -2.0, -0.5, 0.0 };
    try testing.expect(VecF64.logValidate(&v));
}

test "VecF64.logValidate rejects positive values" {
    const v = [_]f64{ -1.0, 0.5, -2.0 };
    try testing.expect(!VecF64.logValidate(&v));
}

test "VecF64.logValidate rejects NaN" {
    const v = [_]f64{ -1.0, std.math.nan(f64) };
    try testing.expect(!VecF64.logValidate(&v));
}

test "VecF64.logValidate rejects negative Inf" {
    const v = [_]f64{ -1.0, -std.math.inf(f64) };
    try testing.expect(!VecF64.logValidate(&v));
}

test "VecF64.logValidate empty" {
    const v = [_]f64{};
    try testing.expect(VecF64.logValidate(&v));
}

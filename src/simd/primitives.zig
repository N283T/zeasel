const std = @import("std");

// ============================================================================
// Lane shift operations (for Farrar's striped DP)
// ============================================================================

/// Shift vector lanes right by 1, filling the vacated leftmost lane with `fill`.
/// Example: rightShift({A,B,C,D}, X) = {X,A,B,C}
pub fn rightShift(comptime N: comptime_int, comptime T: type, v: @Vector(N, T), fill: T) @Vector(N, T) {
    const fill_vec: @Vector(N, T) = @splat(fill);
    comptime var indices: [N]i32 = undefined;
    indices[0] = ~@as(i32, 0); // Take from fill_vec index 0
    inline for (1..N) |i| {
        indices[i] = @intCast(i - 1); // Take from v
    }
    return @shuffle(T, v, fill_vec, indices);
}

/// Shift vector lanes left by 1, filling the vacated rightmost lane with `fill`.
/// Example: leftShift({A,B,C,D}, X) = {B,C,D,X}
pub fn leftShift(comptime N: comptime_int, comptime T: type, v: @Vector(N, T), fill: T) @Vector(N, T) {
    const fill_vec: @Vector(N, T) = @splat(fill);
    comptime var indices: [N]i32 = undefined;
    inline for (0..N - 1) |i| {
        indices[i] = @intCast(i + 1); // Take from v
    }
    indices[N - 1] = ~@as(i32, 0); // Take from fill_vec index 0
    return @shuffle(T, v, fill_vec, indices);
}

// ============================================================================
// Horizontal reductions
// ============================================================================

/// Horizontal max across all lanes.
pub fn hmax(comptime N: comptime_int, comptime T: type, v: @Vector(N, T)) T {
    return @reduce(.Max, v);
}

/// Horizontal min across all lanes.
pub fn hmin(comptime N: comptime_int, comptime T: type, v: @Vector(N, T)) T {
    return @reduce(.Min, v);
}

/// Horizontal sum across all lanes.
pub fn hsum(comptime N: comptime_int, comptime T: type, v: @Vector(N, T)) T {
    return @reduce(.Add, v);
}

// ============================================================================
// Comparison predicates
// ============================================================================

/// Returns true if any lane in `a` is greater than the corresponding lane in `b`.
pub fn anyGt(comptime N: comptime_int, comptime T: type, a: @Vector(N, T), b: @Vector(N, T)) bool {
    const cmp = a > b;
    return @reduce(.Or, cmp);
}

/// Returns true if any lane in `a` is greater than or equal to the corresponding lane in `b`.
pub fn anyGte(comptime N: comptime_int, comptime T: type, a: @Vector(N, T), b: @Vector(N, T)) bool {
    const cmp = a >= b;
    return @reduce(.Or, cmp);
}

// ============================================================================
// Select / blend (branchless conditional)
// ============================================================================

/// Branchless select: for each lane, result[i] = if mask[i] then b[i] else a[i].
pub fn select(comptime N: comptime_int, comptime T: type, mask: @Vector(N, bool), a: @Vector(N, T), b: @Vector(N, T)) @Vector(N, T) {
    return @select(T, mask, b, a);
}

// ============================================================================
// Saturating arithmetic
// ============================================================================

/// Saturating add for integer vectors.
pub fn addSat(comptime N: comptime_int, comptime T: type, a: @Vector(N, T), b: @Vector(N, T)) @Vector(N, T) {
    return a +| b;
}

/// Saturating subtract for integer vectors.
pub fn subSat(comptime N: comptime_int, comptime T: type, a: @Vector(N, T), b: @Vector(N, T)) @Vector(N, T) {
    return a -| b;
}

// ============================================================================
// Tests
// ============================================================================

fn expectVecEqual(comptime N: comptime_int, comptime T: type, expected: @Vector(N, T), actual: @Vector(N, T)) !void {
    const exp_arr: [N]T = expected;
    const act_arr: [N]T = actual;
    for (0..N) |i| {
        if (comptime @typeInfo(T) == .float) {
            if (std.math.isNan(exp_arr[i]) and std.math.isNan(act_arr[i])) continue;
            if (std.math.isInf(exp_arr[i]) and std.math.isInf(act_arr[i])) {
                if ((exp_arr[i] > 0) == (act_arr[i] > 0)) continue;
            }
        }
        try std.testing.expectEqual(exp_arr[i], act_arr[i]);
    }
}

// ---------------------------------------------------------------------------
// rightShift tests
// ---------------------------------------------------------------------------

test "rightShift: 4-wide f32" {
    const v: @Vector(4, f32) = .{ 1.0, 2.0, 3.0, 4.0 };
    const result = rightShift(4, f32, v, -std.math.inf(f32));
    try expectVecEqual(4, f32, .{ -std.math.inf(f32), 1.0, 2.0, 3.0 }, result);
}

test "rightShift: 4-wide i32" {
    const v: @Vector(4, i32) = .{ 10, 20, 30, 40 };
    const result = rightShift(4, i32, v, -999);
    try expectVecEqual(4, i32, .{ -999, 10, 20, 30 }, result);
}

test "rightShift: 16-wide u8 (MSV filter width)" {
    const v: @Vector(16, u8) = .{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    const result = rightShift(16, u8, v, 0);
    try expectVecEqual(16, u8, .{ 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 }, result);
}

test "rightShift: 8-wide i16" {
    const v: @Vector(8, i16) = .{ 100, 200, 300, 400, 500, 600, 700, 800 };
    const result = rightShift(8, i16, v, 0);
    try expectVecEqual(8, i16, .{ 0, 100, 200, 300, 400, 500, 600, 700 }, result);
}

// ---------------------------------------------------------------------------
// leftShift tests
// ---------------------------------------------------------------------------

test "leftShift: 4-wide i16" {
    const v: @Vector(4, i16) = .{ 10, 20, 30, 40 };
    const result = leftShift(4, i16, v, 0);
    try expectVecEqual(4, i16, .{ 20, 30, 40, 0 }, result);
}

test "leftShift: 4-wide f32" {
    const v: @Vector(4, f32) = .{ 1.0, 2.0, 3.0, 4.0 };
    const result = leftShift(4, f32, v, 0.0);
    try expectVecEqual(4, f32, .{ 2.0, 3.0, 4.0, 0.0 }, result);
}

test "leftShift: 16-wide u8" {
    const v: @Vector(16, u8) = .{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    const result = leftShift(16, u8, v, 255);
    try expectVecEqual(16, u8, .{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 255 }, result);
}

// ---------------------------------------------------------------------------
// leftShift / rightShift roundtrip
// ---------------------------------------------------------------------------

test "rightShift then leftShift roundtrip" {
    const v: @Vector(4, i32) = .{ 10, 20, 30, 40 };
    const shifted = rightShift(4, i32, v, 0);
    const back = leftShift(4, i32, shifted, 0);
    // After rightShift: {0, 10, 20, 30}, leftShift: {10, 20, 30, 0}
    // Original last element (40) is lost, last lane becomes fill (0)
    try expectVecEqual(4, i32, .{ 10, 20, 30, 0 }, back);
}

// ---------------------------------------------------------------------------
// hmax tests
// ---------------------------------------------------------------------------

test "hmax: u8 16-wide" {
    const v: @Vector(16, u8) = .{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 255 };
    try std.testing.expectEqual(@as(u8, 255), hmax(16, u8, v));
}

test "hmax: i16 8-wide" {
    const v: @Vector(8, i16) = .{ -100, 50, 30, -20, 10, 0, -5, 100 };
    try std.testing.expectEqual(@as(i16, 100), hmax(8, i16, v));
}

test "hmax: f32 4-wide" {
    const v: @Vector(4, f32) = .{ 1.5, -3.0, 7.25, 2.0 };
    try std.testing.expectEqual(@as(f32, 7.25), hmax(4, f32, v));
}

test "hmax: i8 16-wide with negatives" {
    const v: @Vector(16, i8) = .{ -128, -1, -50, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14 };
    try std.testing.expectEqual(@as(i8, -1), hmax(16, i8, v));
}

// ---------------------------------------------------------------------------
// hmin tests
// ---------------------------------------------------------------------------

test "hmin: u8 16-wide" {
    const v: @Vector(16, u8) = .{ 10, 20, 3, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 255 };
    try std.testing.expectEqual(@as(u8, 3), hmin(16, u8, v));
}

test "hmin: i16 8-wide" {
    const v: @Vector(8, i16) = .{ 100, -200, 300, 400, 500, 600, 700, 800 };
    try std.testing.expectEqual(@as(i16, -200), hmin(8, i16, v));
}

test "hmin: f32 4-wide" {
    const v: @Vector(4, f32) = .{ 1.5, -3.0, 7.25, 2.0 };
    try std.testing.expectEqual(@as(f32, -3.0), hmin(4, f32, v));
}

// ---------------------------------------------------------------------------
// hsum tests
// ---------------------------------------------------------------------------

test "hsum: i32 4-wide" {
    const v: @Vector(4, i32) = .{ 1, 2, 3, 4 };
    try std.testing.expectEqual(@as(i32, 10), hsum(4, i32, v));
}

test "hsum: f32 4-wide" {
    const v: @Vector(4, f32) = .{ 1.0, 2.0, 3.0, 4.0 };
    try std.testing.expectEqual(@as(f32, 10.0), hsum(4, f32, v));
}

test "hsum: u8 16-wide" {
    const v: @Vector(16, u8) = .{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    try std.testing.expectEqual(@as(u8, 16), hsum(16, u8, v));
}

// ---------------------------------------------------------------------------
// anyGt tests
// ---------------------------------------------------------------------------

test "anyGt: detects threshold exceeded (f32)" {
    const scores: @Vector(4, f32) = .{ 1.0, 2.0, 3.0, 15.0 };
    const threshold: @Vector(4, f32) = @splat(10.0);
    try std.testing.expect(anyGt(4, f32, scores, threshold));
}

test "anyGt: below threshold (f32)" {
    const scores: @Vector(4, f32) = .{ 1.0, 2.0, 3.0, 4.0 };
    const threshold: @Vector(4, f32) = @splat(10.0);
    try std.testing.expect(!anyGt(4, f32, scores, threshold));
}

test "anyGt: equal values not greater" {
    const a: @Vector(4, i32) = .{ 5, 5, 5, 5 };
    const b: @Vector(4, i32) = .{ 5, 5, 5, 5 };
    try std.testing.expect(!anyGt(4, i32, a, b));
}

test "anyGt: u8 16-wide" {
    const a: @Vector(16, u8) = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
    const b: @Vector(16, u8) = @splat(0);
    try std.testing.expect(anyGt(16, u8, a, b));
}

// ---------------------------------------------------------------------------
// anyGte tests
// ---------------------------------------------------------------------------

test "anyGte: equal values detected" {
    const a: @Vector(4, i32) = .{ 1, 2, 3, 5 };
    const b: @Vector(4, i32) = .{ 10, 10, 10, 5 };
    try std.testing.expect(anyGte(4, i32, a, b));
}

test "anyGte: strictly less" {
    const a: @Vector(4, i32) = .{ 1, 2, 3, 4 };
    const b: @Vector(4, i32) = .{ 10, 10, 10, 10 };
    try std.testing.expect(!anyGte(4, i32, a, b));
}

// ---------------------------------------------------------------------------
// select tests
// ---------------------------------------------------------------------------

test "select: branchless max (f32)" {
    const a: @Vector(4, f32) = .{ 1.0, 5.0, 3.0, 7.0 };
    const b: @Vector(4, f32) = .{ 4.0, 2.0, 6.0, 1.0 };
    const mask = a > b;
    const result = select(4, f32, mask, b, a); // where a > b pick a, else b
    try expectVecEqual(4, f32, .{ 4.0, 5.0, 6.0, 7.0 }, result);
}

test "select: all true" {
    const a: @Vector(4, i32) = .{ 1, 2, 3, 4 };
    const b: @Vector(4, i32) = .{ 10, 20, 30, 40 };
    const mask: @Vector(4, bool) = .{ true, true, true, true };
    const result = select(4, i32, mask, a, b);
    try expectVecEqual(4, i32, .{ 10, 20, 30, 40 }, result);
}

test "select: all false" {
    const a: @Vector(4, i32) = .{ 1, 2, 3, 4 };
    const b: @Vector(4, i32) = .{ 10, 20, 30, 40 };
    const mask: @Vector(4, bool) = .{ false, false, false, false };
    const result = select(4, i32, mask, a, b);
    try expectVecEqual(4, i32, .{ 1, 2, 3, 4 }, result);
}

test "select: mixed mask i16" {
    const a: @Vector(8, i16) = .{ 1, 2, 3, 4, 5, 6, 7, 8 };
    const b: @Vector(8, i16) = .{ 10, 20, 30, 40, 50, 60, 70, 80 };
    const mask: @Vector(8, bool) = .{ true, false, true, false, true, false, true, false };
    const result = select(8, i16, mask, a, b);
    try expectVecEqual(8, i16, .{ 10, 2, 30, 4, 50, 6, 70, 8 }, result);
}

// ---------------------------------------------------------------------------
// addSat tests
// ---------------------------------------------------------------------------

test "addSat: u8 saturating" {
    const a: @Vector(4, u8) = .{ 200, 100, 50, 255 };
    const b: @Vector(4, u8) = .{ 100, 100, 100, 1 };
    const result = addSat(4, u8, a, b);
    try expectVecEqual(4, u8, .{ 255, 200, 150, 255 }, result);
}

test "addSat: i8 saturating" {
    const a: @Vector(4, i8) = .{ 100, -100, 50, 127 };
    const b: @Vector(4, i8) = .{ 50, -50, -100, 1 };
    const result = addSat(4, i8, a, b);
    try expectVecEqual(4, i8, .{ 127, -128, -50, 127 }, result);
}

test "addSat: i16 saturating" {
    const a: @Vector(8, i16) = .{ 32000, -32000, 100, 0, 0, 0, 0, 0 };
    const b: @Vector(8, i16) = .{ 1000, -1000, 200, 0, 0, 0, 0, 0 };
    const result = addSat(8, i16, a, b);
    try expectVecEqual(8, i16, .{ 32767, -32768, 300, 0, 0, 0, 0, 0 }, result);
}

test "addSat: u8 16-wide MSV style" {
    var a: @Vector(16, u8) = @splat(200);
    var b: @Vector(16, u8) = @splat(100);
    const result = addSat(16, u8, a, b);
    const expected: @Vector(16, u8) = @splat(255);
    try expectVecEqual(16, u8, expected, result);

    // Non-saturating case
    a = @splat(50);
    b = @splat(30);
    const result2 = addSat(16, u8, a, b);
    const expected2: @Vector(16, u8) = @splat(80);
    try expectVecEqual(16, u8, expected2, result2);
}

// ---------------------------------------------------------------------------
// subSat tests
// ---------------------------------------------------------------------------

test "subSat: u8 saturating" {
    const a: @Vector(4, u8) = .{ 200, 100, 50, 0 };
    const b: @Vector(4, u8) = .{ 100, 100, 100, 1 };
    const result = subSat(4, u8, a, b);
    try expectVecEqual(4, u8, .{ 100, 0, 0, 0 }, result);
}

test "subSat: i8 saturating" {
    const a: @Vector(4, i8) = .{ -100, 100, 0, -128 };
    const b: @Vector(4, i8) = .{ 50, -50, 0, 1 };
    const result = subSat(4, i8, a, b);
    try expectVecEqual(4, i8, .{ -128, 127, 0, -128 }, result);
}

test "subSat: u8 16-wide MSV style" {
    const a: @Vector(16, u8) = @splat(10);
    const b: @Vector(16, u8) = @splat(20);
    const result = subSat(16, u8, a, b);
    const expected: @Vector(16, u8) = @splat(0);
    try expectVecEqual(16, u8, expected, result);
}

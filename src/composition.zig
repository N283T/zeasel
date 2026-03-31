// Background amino acid residue frequency vectors.
//
// These are standard background compositions used by substitution matrices
// and HMMER profile HMM parameterization. All vectors are in alphabetic order
// matching the standard zeasel amino acid alphabet (A=0, C=1, ..., Y=19).

const std = @import("std");

/// BLOSUM62 background frequencies (Henikoff & Henikoff, 1992).
pub const bl62: [20]f64 = .{
    0.074, // A
    0.025, // C
    0.054, // D
    0.054, // E
    0.047, // F
    0.074, // G
    0.026, // H
    0.068, // I
    0.058, // K
    0.099, // L
    0.025, // M
    0.045, // N
    0.039, // P
    0.034, // Q
    0.052, // R
    0.057, // S
    0.051, // T
    0.073, // V
    0.013, // W
    0.032, // Y
};

/// WAG model background frequencies (Whelan & Goldman, 2001).
pub const wag: [20]f64 = .{
    0.086628, // A
    0.019308, // C
    0.057045, // D
    0.058059, // E
    0.038432, // F
    0.083252, // G
    0.024431, // H
    0.048466, // I
    0.062029, // K
    0.086209, // L
    0.019503, // M
    0.039089, // N
    0.045763, // P
    0.036728, // Q
    0.043972, // R
    0.069518, // S
    0.061013, // T
    0.070896, // V
    0.014386, // W
    0.035274, // Y
};

/// Swiss-Prot release 34 background frequencies (21.2M residues).
pub const sw34: [20]f64 = .{
    0.075520, // A
    0.016973, // C
    0.053029, // D
    0.063204, // E
    0.040762, // F
    0.068448, // G
    0.022406, // H
    0.057284, // I
    0.059398, // K
    0.093399, // L
    0.023569, // M
    0.045293, // N
    0.049262, // P
    0.040231, // Q
    0.051573, // R
    0.072214, // S
    0.057454, // T
    0.065252, // V
    0.012513, // W
    0.031985, // Y
};

/// Swiss-Prot release 50.8 background frequencies (86.0M residues, Oct 2006).
pub const sw50: [20]f64 = .{
    0.0787945, // A
    0.0151600, // C
    0.0535222, // D
    0.0668298, // E
    0.0397062, // F
    0.0695071, // G
    0.0229198, // H
    0.0590092, // I
    0.0594422, // K
    0.0963728, // L
    0.0237718, // M
    0.0414386, // N
    0.0482904, // P
    0.0395639, // Q
    0.0540978, // R
    0.0683364, // S
    0.0540687, // T
    0.0673417, // V
    0.0114135, // W
    0.0304133, // Y
};

/// Uniform background frequencies (1/20 each).
pub const uniform: [20]f64 = .{0.05} ** 20;

// --- Tests ---

test "bl62: sums to ~1.0" {
    var sum: f64 = 0;
    for (bl62) |f| sum += f;
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), sum, 1e-6);
}

test "wag: sums to ~1.0" {
    var sum: f64 = 0;
    for (wag) |f| sum += f;
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), sum, 1e-4);
}

test "sw34: sums to ~1.0" {
    var sum: f64 = 0;
    for (sw34) |f| sum += f;
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), sum, 1e-3);
}

test "sw50: sums to ~1.0" {
    var sum: f64 = 0;
    for (sw50) |f| sum += f;
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), sum, 1e-4);
}

test "uniform: sums to 1.0" {
    var sum: f64 = 0;
    for (uniform) |f| sum += f;
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), sum, 1e-9);
}

test "bl62: all positive" {
    for (bl62) |f| {
        try std.testing.expect(f > 0);
    }
}

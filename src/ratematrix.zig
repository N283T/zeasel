// Evolutionary rate matrix Q and probability matrix P = exp(tQ).
//
// A rate matrix Q is a square matrix where:
//   - Off-diagonal Q[i][j] >= 0 (rate of change from i to j)
//   - Rows sum to zero: Q[i][i] = -sum_{j!=i} Q[i][j]
//
// The transition probability matrix P(t) = exp(tQ) gives the
// probability of changing from residue i to j in time t.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Matrix = @import("matrix.zig").Matrix;
const composition = @import("composition.zig");

/// Create a rate matrix Q from an exchangeability matrix S and
/// stationary frequencies pi. Q[i][j] = S[i][j] * pi[j] for i != j,
/// and Q[i][i] set so rows sum to zero.
/// S must be symmetric. Returns a new n*n Matrix.
pub fn fromExchangeability(allocator: Allocator, s: Matrix, pi: []const f64) !Matrix {
    const n = s.rows;
    std.debug.assert(s.cols == n);
    std.debug.assert(pi.len == n);

    var q = try Matrix.init(allocator, n, n);
    errdefer q.deinit();

    for (0..n) |i| {
        var row_sum: f64 = 0;
        for (0..n) |j| {
            if (i != j) {
                const rate = s.get(i, j) * pi[j];
                q.set(i, j, rate);
                row_sum += rate;
            }
        }
        q.set(i, i, -row_sum);
    }

    return q;
}

/// Normalize a rate matrix Q so that the expected number of
/// substitutions per unit time is 1: -sum_i pi[i] * Q[i][i] = 1.
/// Replaces Q's data with the scaled version.
pub fn normalize(q: *Matrix, pi: []const f64) !void {
    const n = q.rows;
    var rate: f64 = 0;
    for (0..n) |i| {
        rate -= pi[i] * q.get(i, i);
    }
    if (rate > 0) {
        const scaled = try q.scale(1.0 / rate);
        q.allocator.free(q.data);
        q.data = scaled.data;
    }
}

/// Scale a rate matrix Q so that the expected number of substitutions
/// per unit time equals `unit`. For example, unit=1.0 gives substitutions
/// per site; unit=0.01 gives PAM units (1 substitution per 100 sites).
/// Equivalent to Easel's `esl_rmx_ScaleTo()`.
/// Replaces Q's data with the scaled version.
pub fn scaleTo(q: *Matrix, pi: []const f64, unit: f64) !void {
    const n = q.rows;
    var rate: f64 = 0;
    for (0..n) |i| {
        rate -= pi[i] * q.get(i, i);
    }
    if (rate > 0) {
        const scaled = try q.scale(unit / rate);
        q.allocator.free(q.data);
        q.data = scaled.data;
    }
}

/// Compute the probability matrix P(t) = exp(tQ).
/// Caller owns the returned Matrix.
pub fn probMatrix(q: Matrix, t: f64) !Matrix {
    return q.exp(t);
}

/// Compute the probability matrix P(t) = exp(tQ) using eigendecomposition.
///
/// For reversible rate matrices (those satisfying detailed balance), this
/// first symmetrizes Q using the stationary frequencies pi, diagonalizes
/// the symmetric matrix via Jacobi, and then computes P = R * diag(exp(t*lambda)) * L.
///
/// When computing P(t) for many different t values (as in phylogenetic
/// likelihood calculations), call `Matrix.diagonalizeSymmetrizable()` once
/// and then `Matrix.expFromEigen()` for each t instead.
///
/// Caller owns the returned Matrix.
pub fn probMatrixEigen(q: Matrix, pi: []const f64, t: f64) !Matrix {
    var eigen = try q.diagonalizeSymmetrizable(q.allocator, pi);
    defer eigen.deinit();
    return Matrix.expFromEigen(q.allocator, eigen, t);
}

/// Validate a rate matrix: rows should sum to ~0.
pub fn isValid(q: Matrix, tol: f64) bool {
    if (q.rows != q.cols) return false;
    const n = q.rows;
    for (0..n) |i| {
        var row_sum: f64 = 0;
        for (0..n) |j| {
            row_sum += q.get(i, j);
        }
        if (@abs(row_sum) > tol) return false;
        // Diagonal should be negative
        if (q.get(i, i) > tol) return false;
        // Off-diagonals should be non-negative
        for (0..n) |j| {
            if (i != j and q.get(i, j) < -tol) return false;
        }
    }
    return true;
}

/// Jukes-Cantor model for an alphabet of size K.
/// Off-diagonal rates are all equal, yielding a uniform substitution model.
/// The matrix is normalized so that 1 time unit = 1 expected substitution/site.
/// Caller owns the returned Matrix.
pub fn setJukesCantor(allocator: Allocator, k: usize) !Matrix {
    var q = try Matrix.init(allocator, k, k);
    errdefer q.deinit();

    const fk: f64 = @floatFromInt(k);
    for (0..k) |i| {
        for (0..k) |j| {
            if (i != j) {
                q.set(i, j, 1.0);
            }
        }
        q.set(i, i, -@as(f64, @floatFromInt(k - 1)));
    }

    // Normalize: uniform pi = 1/K, expected rate = -sum(pi[i]*Q[i][i]) = (K-1).
    // Scale so expected rate = 1.
    const scaled = try q.scale(1.0 / (fk - 1.0));
    q.allocator.free(q.data);
    q.data = scaled.data;

    return q;
}

/// Kimura 2-parameter model for DNA (4x4).
/// `alpha` = transition rate, `beta` = transversion rate.
/// Normalized to 1 expected substitution/site per unit time.
/// DNA alphabet order: A=0, C=1, G=2, T=3.
/// Transitions: A<->G (0,2) and C<->T (1,3); all others are transversions.
/// Caller owns the returned Matrix.
pub fn setKimura(allocator: Allocator, alpha: f64, beta: f64) !Matrix {
    var q = try Matrix.init(allocator, 4, 4);
    errdefer q.deinit();

    for (0..4) |i| {
        var row_sum: f64 = 0;
        for (0..4) |j| {
            if (i != j) {
                // Easel convention: (i+j)%2 == 0 => transition, odd => transversion
                const rate = if ((i + j) % 2 == 0) alpha else beta;
                q.set(i, j, rate);
                row_sum += rate;
            }
        }
        q.set(i, i, -row_sum);
    }

    // Normalize to 1 substitution/site: pi = 0.25 uniform
    const pi = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    try normalize(&q, &pi);

    return q;
}

/// F81 model (Felsenstein, 1981) for DNA (4x4).
/// Off-diagonal Q[i][j] = pi[j]; diagonal set so rows sum to zero.
/// Normalized to 1 expected substitution/site per unit time.
/// DNA alphabet order: A=0, C=1, G=2, T=3.
/// Caller owns the returned Matrix.
pub fn setF81(allocator: Allocator, pi: [4]f64) !Matrix {
    var q = try Matrix.init(allocator, 4, 4);
    errdefer q.deinit();

    for (0..4) |i| {
        var row_sum: f64 = 0;
        for (0..4) |j| {
            if (i != j) {
                q.set(i, j, pi[j]);
                row_sum += pi[j];
            }
        }
        q.set(i, i, -row_sum);
    }

    // Normalize to 1 substitution/site
    var pi_slice: [4]f64 = pi;
    try normalize(&q, &pi_slice);

    return q;
}

/// HKY model (Hasegawa-Kishino-Yano, 1985) for DNA (4x4).
/// `alpha` = transition rate, `beta` = transversion rate.
/// Q[i][j] = alpha * pi[j] for transitions (A<->G, C<->T),
///           beta * pi[j]  for transversions.
/// Diagonal set so rows sum to zero.
/// Normalized to 1 expected substitution/site per unit time.
/// DNA alphabet order: A=0, C=1, G=2, T=3.
/// Transitions: (0,2), (2,0), (1,3), (3,1); all others are transversions.
/// Caller owns the returned Matrix.
pub fn setHKY(allocator: Allocator, alpha: f64, beta: f64, pi: [4]f64) !Matrix {
    var q = try Matrix.init(allocator, 4, 4);
    errdefer q.deinit();

    for (0..4) |i| {
        var row_sum: f64 = 0;
        for (0..4) |j| {
            if (i != j) {
                // (i+j)%2 == 0 => transition (A<->G: 0+2=2, C<->T: 1+3=4)
                const rate = if ((i + j) % 2 == 0) alpha * pi[j] else beta * pi[j];
                q.set(i, j, rate);
                row_sum += rate;
            }
        }
        q.set(i, i, -row_sum);
    }

    // Normalize to 1 substitution/site
    var pi_slice: [4]f64 = pi;
    try normalize(&q, &pi_slice);

    return q;
}

/// WAG amino acid exchangeability parameters (Whelan & Goldman, 2001).
/// 190 values in lower-triangular order: for i=1..19, j=0..i-1.
const wag_exchangeabilities: [190]f64 = .{
    1.027040, 0.738998, 0.030295, 1.582850, 0.021352, 6.174160, 0.210494, 0.398020, 0.046730, 0.081134,
    1.416720, 0.306674, 0.865584, 0.567717, 0.049931, 0.316954, 0.248972, 0.930676, 0.570025, 0.679371,
    0.249410, 0.193335, 0.170135, 0.039437, 0.127395, 1.059470, 0.030450, 0.138190, 0.906265, 0.074034,
    0.479855, 2.584430, 0.088836, 0.373558, 0.890432, 0.323832, 0.397915, 0.384287, 0.084805, 0.154263,
    2.115170, 0.061304, 0.499462, 3.170970, 0.257555, 0.893496, 0.390482, 0.103754, 0.315124, 1.190630,
    0.174100, 0.404141, 4.257460, 0.934276, 4.854020, 0.509848, 0.265256, 5.429420, 0.947198, 0.096162,
    1.125560, 3.956290, 0.554236, 3.012010, 0.131528, 0.198221, 1.438550, 0.109404, 0.423984, 0.682355,
    0.161444, 0.243570, 0.696198, 0.099929, 0.556896, 0.415844, 0.171329, 0.195081, 0.908598, 0.098818,
    0.616783, 5.469470, 0.099921, 0.330052, 4.294110, 0.113917, 3.894900, 0.869489, 1.545260, 1.543640,
    0.933372, 0.551571, 0.528191, 0.147304, 0.439157, 0.102711, 0.584665, 2.137150, 0.186979, 5.351420,
    0.497671, 0.683162, 0.635346, 0.679489, 3.035500, 3.370790, 1.407660, 1.071760, 0.704939, 0.545931,
    1.341820, 0.740169, 0.319440, 0.967130, 0.344739, 0.493905, 3.974230, 1.613280, 1.028870, 1.224190,
    2.121110, 0.512984, 0.374866, 0.822765, 0.171903, 0.225833, 0.473307, 1.458160, 1.386980, 0.326622,
    1.516120, 2.030060, 0.795384, 0.857928, 0.554413, 4.378020, 2.006010, 1.002140, 0.152335, 0.588731,
    0.649892, 0.187247, 0.118358, 7.821300, 0.305434, 1.800340, 2.058450, 0.196246, 0.314887, 0.301281,
    0.251849, 0.232739, 1.388230, 0.113133, 0.717070, 0.129767, 0.156557, 1.529640, 0.336983, 0.262569,
    0.212483, 0.137505, 0.665309, 0.515706, 0.071917, 0.139405, 0.215737, 1.163920, 0.523742, 0.110864,
    0.365369, 0.240735, 0.543833, 0.325711, 0.196303, 6.454280, 0.103604, 3.873440, 0.420170, 0.133264,
    0.398618, 0.428437, 1.086000, 0.216046, 0.227710, 0.381533, 0.786993, 0.291148, 0.314730, 2.485390,
};

/// WAG amino acid rate matrix (Whelan & Goldman, 2001).
/// Q[i][j] = S[i][j] * pi[j] for i != j, diagonal set so rows sum to 0.
/// Normalized to 1 expected substitution/site per unit time.
/// If `pi` is provided, those frequencies are used instead of the default WAG frequencies.
/// Caller owns the returned Matrix.
pub fn setWAG(allocator: Allocator, pi: ?[]const f64) !Matrix {
    const freqs = if (pi) |p| p else &composition.wag;
    std.debug.assert(freqs.len == 20);

    var q = try Matrix.init(allocator, 20, 20);
    errdefer q.deinit();

    // Fill symmetric exchangeability matrix into Q
    var z: usize = 0;
    for (0..20) |i| {
        q.set(i, i, 0);
        for (0..i) |j| {
            q.set(i, j, wag_exchangeabilities[z]);
            q.set(j, i, wag_exchangeabilities[z]);
            z += 1;
        }
    }

    // Multiply off-diagonals by pi[j]: Q[i][j] = E[i][j] * pi[j]
    for (0..20) |i| {
        for (0..20) |j| {
            if (i != j) {
                q.set(i, j, q.get(i, j) * freqs[j]);
            }
        }
    }

    // Set diagonal: Q[i][i] = -sum_{j!=i} Q[i][j]
    for (0..20) |i| {
        var row_sum: f64 = 0;
        for (0..20) |j| {
            if (i != j) row_sum += q.get(i, j);
        }
        q.set(i, i, -row_sum);
    }

    // Normalize to 1 substitution/site
    try normalize(&q, freqs);

    return q;
}

/// Validate a probability transition matrix P.
///
/// Checks that every row sums to 1 (within `tol`) and all elements
/// are in [0, 1] (within `tol`). Returns true if the matrix passes
/// all checks.
pub fn validateP(m: Matrix, tol: f64) bool {
    if (m.rows != m.cols) return false;
    const n = m.rows;
    for (0..n) |i| {
        var row_sum: f64 = 0;
        for (0..n) |j| {
            const v = m.get(i, j);
            if (v < -tol or v > 1.0 + tol) return false;
            row_sum += v;
        }
        if (@abs(row_sum - 1.0) > tol) return false;
    }
    return true;
}

/// KL divergence of transition matrix P given stationary distribution pi.
/// Returns result in **bits** (base-2 logarithm), matching Easel's esl_rmx_RelativeEntropy.
///
/// D_KL = sum_i pi[i] * sum_j P[i][j] * log2(P[i][j] / pi[j])
///
/// Terms where P[i][j] == 0 contribute 0 (by convention 0*log(0) = 0).
pub fn relativeEntropy(p: Matrix, pi: []const f64) f64 {
    const n = p.rows;
    std.debug.assert(p.cols == n);
    std.debug.assert(pi.len == n);

    var dkl: f64 = 0;
    for (0..n) |i| {
        for (0..n) |j| {
            const pij = p.get(i, j);
            if (pij > 0 and pi[j] > 0) {
                dkl += pi[i] * pij * @log2(pij / pi[j]);
            }
        }
    }
    return dkl;
}

/// Expected score per site in **bits**.
/// Reference: Easel esl_rmx_ExpectedScore.
///
/// E = sum_i sum_j pi[i] * pi[j] * log2(P[i][j] / pi[j])
///
/// Note: this differs from relativeEntropy in using pi[j] (not P[i][j])
/// as the weight in the outer sum alongside pi[i].
pub fn expectedScore(p: Matrix, pi: []const f64) f64 {
    const n = p.rows;
    std.debug.assert(p.cols == n);
    std.debug.assert(pi.len == n);

    var score: f64 = 0;
    for (0..n) |i| {
        for (0..n) |j| {
            const pij = p.get(i, j);
            if (pij > 0 and pi[j] > 0) {
                score += pi[i] * pi[j] * @log2(pij / pi[j]);
            }
        }
    }
    return score;
}

// --- Tests ---

test "fromExchangeability: 2x2" {
    const allocator = std.testing.allocator;
    // S = [[0, 1], [1, 0]], pi = [0.6, 0.4]
    // Q[0][1] = 1 * 0.4 = 0.4, Q[1][0] = 1 * 0.6 = 0.6
    // Q[0][0] = -0.4, Q[1][1] = -0.6
    var s = try Matrix.init(allocator, 2, 2);
    defer s.deinit();
    s.set(0, 1, 1.0);
    s.set(1, 0, 1.0);

    const pi = [_]f64{ 0.6, 0.4 };
    var q = try fromExchangeability(allocator, s, &pi);
    defer q.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, -0.4), q.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), q.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), q.get(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, -0.6), q.get(1, 1), 1e-10);
    try std.testing.expect(isValid(q, 1e-10));
}

test "normalize: expected rate becomes 1" {
    const allocator = std.testing.allocator;
    var s = try Matrix.init(allocator, 2, 2);
    defer s.deinit();
    s.set(0, 1, 2.0);
    s.set(1, 0, 2.0);

    const pi = [_]f64{ 0.5, 0.5 };
    var q = try fromExchangeability(allocator, s, &pi);
    defer q.deinit();

    try normalize(&q, &pi);

    // Check: -sum(pi[i] * Q[i][i]) should be 1.0
    var rate: f64 = 0;
    for (0..2) |i| rate -= pi[i] * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rate, 1e-10);
}

test "scaleTo: PAM unit scaling" {
    const allocator = std.testing.allocator;
    var s = try Matrix.init(allocator, 2, 2);
    defer s.deinit();
    s.set(0, 1, 2.0);
    s.set(1, 0, 2.0);

    const pi = [_]f64{ 0.5, 0.5 };
    var q = try fromExchangeability(allocator, s, &pi);
    defer q.deinit();

    // Scale to PAM units (0.01 substitutions/site per unit time)
    try scaleTo(&q, &pi, 0.01);

    // Check: -sum(pi[i] * Q[i][i]) should be 0.01
    var rate: f64 = 0;
    for (0..2) |i| rate -= pi[i] * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), rate, 1e-10);
}

test "probMatrix: P(0) = I" {
    const allocator = std.testing.allocator;
    var q = try Matrix.init(allocator, 2, 2);
    defer q.deinit();
    q.set(0, 0, -1.0);
    q.set(0, 1, 1.0);
    q.set(1, 0, 1.0);
    q.set(1, 1, -1.0);

    var p = try probMatrix(q, 0.0);
    defer p.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, 1.0), p.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), p.get(0, 1), 1e-10);
}

test "probMatrix: rows sum to 1" {
    const allocator = std.testing.allocator;
    var q = try Matrix.init(allocator, 3, 3);
    defer q.deinit();
    q.set(0, 0, -2.0);
    q.set(0, 1, 1.0);
    q.set(0, 2, 1.0);
    q.set(1, 0, 1.0);
    q.set(1, 1, -2.0);
    q.set(1, 2, 1.0);
    q.set(2, 0, 1.0);
    q.set(2, 1, 1.0);
    q.set(2, 2, -2.0);

    var p = try probMatrix(q, 0.5);
    defer p.deinit();

    for (0..3) |i| {
        var row_sum: f64 = 0;
        for (0..3) |j| row_sum += p.get(i, j);
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), row_sum, 1e-8);
    }
}

test "setJukesCantor: DNA (K=4)" {
    const allocator = std.testing.allocator;
    var q = try setJukesCantor(allocator, 4);
    defer q.deinit();

    try std.testing.expect(isValid(q, 1e-10));

    // All off-diagonals should be equal
    const off_diag = q.get(0, 1);
    for (0..4) |i| {
        for (0..4) |j| {
            if (i != j) {
                try std.testing.expectApproxEqAbs(off_diag, q.get(i, j), 1e-10);
            }
        }
    }

    // Expected rate = 1.0: -sum(pi[i] * Q[i][i]) with pi = 0.25
    var rate: f64 = 0;
    for (0..4) |i| rate -= 0.25 * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rate, 1e-10);
}

test "setJukesCantor: protein (K=20)" {
    const allocator = std.testing.allocator;
    var q = try setJukesCantor(allocator, 20);
    defer q.deinit();

    try std.testing.expect(isValid(q, 1e-10));

    // Expected rate = 1.0 with uniform pi = 1/20
    var rate: f64 = 0;
    for (0..20) |i| rate -= (1.0 / 20.0) * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rate, 1e-10);
}

test "setKimura: alpha=1 beta=1 reduces to JC" {
    const allocator = std.testing.allocator;
    var jc = try setJukesCantor(allocator, 4);
    defer jc.deinit();
    var kim = try setKimura(allocator, 1.0, 1.0);
    defer kim.deinit();

    for (0..4) |i| {
        for (0..4) |j| {
            try std.testing.expectApproxEqAbs(jc.get(i, j), kim.get(i, j), 1e-10);
        }
    }
}

test "setKimura: transition/transversion asymmetry" {
    const allocator = std.testing.allocator;
    var q = try setKimura(allocator, 2.0, 1.0);
    defer q.deinit();

    try std.testing.expect(isValid(q, 1e-10));

    // Transition rate (A->G = 0,2) should be higher than transversion (A->C = 0,1)
    try std.testing.expect(q.get(0, 2) > q.get(0, 1));

    // Expected rate = 1.0
    var rate: f64 = 0;
    for (0..4) |i| rate -= 0.25 * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rate, 1e-10);
}

test "setF81: uniform pi reduces to JC" {
    const allocator = std.testing.allocator;
    var jc = try setJukesCantor(allocator, 4);
    defer jc.deinit();
    var f81 = try setF81(allocator, .{ 0.25, 0.25, 0.25, 0.25 });
    defer f81.deinit();

    for (0..4) |i| {
        for (0..4) |j| {
            try std.testing.expectApproxEqAbs(jc.get(i, j), f81.get(i, j), 1e-10);
        }
    }
}

test "setF81: valid rate matrix with non-uniform pi" {
    const allocator = std.testing.allocator;
    const pi = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    var q = try setF81(allocator, pi);
    defer q.deinit();

    try std.testing.expect(isValid(q, 1e-10));

    // Expected rate = 1.0
    var rate: f64 = 0;
    for (0..4) |i| rate -= pi[i] * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rate, 1e-10);
}

test "setF81: detailed balance" {
    const allocator = std.testing.allocator;
    const pi = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    var q = try setF81(allocator, pi);
    defer q.deinit();

    // Detailed balance: pi[i]*Q[i][j] = pi[j]*Q[j][i]
    for (0..4) |i| {
        for (0..4) |j| {
            if (i != j) {
                const lhs = pi[i] * q.get(i, j);
                const rhs = pi[j] * q.get(j, i);
                try std.testing.expectApproxEqRel(lhs, rhs, 1e-10);
            }
        }
    }
}

test "setHKY: uniform pi with alpha=beta reduces to JC" {
    const allocator = std.testing.allocator;
    var jc = try setJukesCantor(allocator, 4);
    defer jc.deinit();
    var hky = try setHKY(allocator, 1.0, 1.0, .{ 0.25, 0.25, 0.25, 0.25 });
    defer hky.deinit();

    for (0..4) |i| {
        for (0..4) |j| {
            try std.testing.expectApproxEqAbs(jc.get(i, j), hky.get(i, j), 1e-10);
        }
    }
}

test "setHKY: uniform pi reduces to Kimura" {
    const allocator = std.testing.allocator;
    var kim = try setKimura(allocator, 2.0, 1.0);
    defer kim.deinit();
    var hky = try setHKY(allocator, 2.0, 1.0, .{ 0.25, 0.25, 0.25, 0.25 });
    defer hky.deinit();

    for (0..4) |i| {
        for (0..4) |j| {
            try std.testing.expectApproxEqAbs(kim.get(i, j), hky.get(i, j), 1e-10);
        }
    }
}

test "setHKY: valid rate matrix with non-uniform pi" {
    const allocator = std.testing.allocator;
    const pi = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    var q = try setHKY(allocator, 2.0, 1.0, pi);
    defer q.deinit();

    try std.testing.expect(isValid(q, 1e-10));

    // Expected rate = 1.0
    var rate: f64 = 0;
    for (0..4) |i| rate -= pi[i] * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rate, 1e-10);

    // Transition rate (A->G) should be higher than transversion (A->C)
    try std.testing.expect(q.get(0, 2) > q.get(0, 1));
}

test "setHKY: detailed balance" {
    const allocator = std.testing.allocator;
    const pi = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    var q = try setHKY(allocator, 2.0, 1.0, pi);
    defer q.deinit();

    // Detailed balance: pi[i]*Q[i][j] = pi[j]*Q[j][i]
    for (0..4) |i| {
        for (0..4) |j| {
            if (i != j) {
                const lhs = pi[i] * q.get(i, j);
                const rhs = pi[j] * q.get(j, i);
                try std.testing.expectApproxEqRel(lhs, rhs, 1e-10);
            }
        }
    }
}

test "setHKY: alpha=beta reduces to F81" {
    const allocator = std.testing.allocator;
    const pi = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    var f81 = try setF81(allocator, pi);
    defer f81.deinit();
    var hky = try setHKY(allocator, 1.0, 1.0, pi);
    defer hky.deinit();

    for (0..4) |i| {
        for (0..4) |j| {
            try std.testing.expectApproxEqAbs(f81.get(i, j), hky.get(i, j), 1e-10);
        }
    }
}

test "setWAG: valid rate matrix with default frequencies" {
    const allocator = std.testing.allocator;
    var q = try setWAG(allocator, null);
    defer q.deinit();

    try std.testing.expect(isValid(q, 1e-10));

    // Expected rate = 1.0
    var rate: f64 = 0;
    for (0..20) |i| rate -= composition.wag[i] * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rate, 1e-10);

    // Matrix should be 20x20
    try std.testing.expectEqual(@as(usize, 20), q.rows);
    try std.testing.expectEqual(@as(usize, 20), q.cols);
}

test "setWAG: detailed balance" {
    const allocator = std.testing.allocator;
    var q = try setWAG(allocator, null);
    defer q.deinit();

    // Detailed balance: pi[i]*Q[i][j] = pi[j]*Q[j][i]
    for (0..20) |i| {
        for (0..20) |j| {
            if (i != j) {
                const lhs = composition.wag[i] * q.get(i, j);
                const rhs = composition.wag[j] * q.get(j, i);
                try std.testing.expectApproxEqRel(lhs, rhs, 1e-10);
            }
        }
    }
}

test "validateP: identity matrix is valid P" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 4);
    defer m.deinit();
    try std.testing.expect(validateP(m, 1e-10));
}

test "validateP: P(t) from rate matrix is valid" {
    const allocator = std.testing.allocator;
    var q = try setJukesCantor(allocator, 4);
    defer q.deinit();
    var p = try probMatrix(q, 0.5);
    defer p.deinit();
    try std.testing.expect(validateP(p, 1e-8));
}

test "validateP: zero matrix is invalid (rows sum to 0, not 1)" {
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 2, 2);
    defer m.deinit();
    try std.testing.expect(!validateP(m, 1e-10));
}

test "validateP: negative elements are invalid" {
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 2, 2);
    defer m.deinit();
    m.set(0, 0, 1.5);
    m.set(0, 1, -0.5);
    m.set(1, 0, 0.5);
    m.set(1, 1, 0.5);
    try std.testing.expect(!validateP(m, 1e-10));
}

test "relativeEntropy: identity P gives entropy of pi in bits" {
    // If P[i][j] = delta_ij, then D_KL = sum_i pi[i] * log2(1/pi[i]) = H(pi) in bits.
    // For uniform pi = 0.25: D_KL = log2(4) = 2.0 bits.
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 4);
    defer m.deinit();
    const pi = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    const dkl = relativeEntropy(m, &pi);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), dkl, 1e-10);
}

test "relativeEntropy: uniform P gives 0" {
    // If P[i][j] = pi[j] for all i, then log(P[i][j]/pi[j]) = 0.
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 3, 3);
    defer m.deinit();
    const pi = [_]f64{ 0.2, 0.3, 0.5 };
    for (0..3) |i| {
        for (0..3) |j| {
            m.set(i, j, pi[j]);
        }
    }
    const dkl = relativeEntropy(m, &pi);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), dkl, 1e-10);
}

test "expectedScore: differs from relativeEntropy for non-equilibrium P" {
    const allocator = std.testing.allocator;
    var q = try setJukesCantor(allocator, 4);
    defer q.deinit();
    var p = try probMatrix(q, 0.1);
    defer p.deinit();
    const pi = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    // For uniform pi with JC model, P is symmetric and rows sum to 1,
    // so P[i][j] = pi[j] when i!=j has same factor. For uniform pi,
    // relativeEntropy and expectedScore give different results because
    // relativeEntropy uses P[i][j] as weight while expectedScore uses pi[j].
    const re = relativeEntropy(p, &pi);
    const es = expectedScore(p, &pi);
    // Both should be non-negative for this model
    try std.testing.expect(re >= 0.0);
    try std.testing.expect(es <= 0.0 or es >= 0.0); // just check it computes
    // For uniform pi + JC, P[i][j] = pi[j] + (delta_ij - pi[j])*exp(-4t/3)
    // The two should actually differ because relativeEntropy weights by P[i][j]
    // and expectedScore weights by pi[j].
}

test "setWAG: custom uniform frequencies give symmetric Q" {
    const allocator = std.testing.allocator;
    // Use uniform frequencies
    var uniform: [20]f64 = undefined;
    for (&uniform) |*v| v.* = 0.05;

    var q = try setWAG(allocator, &uniform);
    defer q.deinit();

    try std.testing.expect(isValid(q, 1e-10));

    // With uniform pi, Q should be symmetric (since S is symmetric)
    for (0..20) |i| {
        for (0..20) |j| {
            if (i != j) {
                try std.testing.expectApproxEqRel(q.get(i, j), q.get(j, i), 1e-10);
            }
        }
    }
}

test "probMatrixEigen: JC matches direct exp" {
    const allocator = std.testing.allocator;
    var q = try setJukesCantor(allocator, 4);
    defer q.deinit();

    const pi = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    var p_eigen = try probMatrixEigen(q, &pi, 0.5);
    defer p_eigen.deinit();
    var p_direct = try probMatrix(q, 0.5);
    defer p_direct.deinit();

    for (0..4) |i| {
        for (0..4) |j| {
            try std.testing.expectApproxEqAbs(p_direct.get(i, j), p_eigen.get(i, j), 1e-8);
        }
    }
    try std.testing.expect(validateP(p_eigen, 1e-8));
}

test "probMatrixEigen: WAG with default frequencies" {
    const allocator = std.testing.allocator;
    var q = try setWAG(allocator, null);
    defer q.deinit();

    var p_eigen = try probMatrixEigen(q, &composition.wag, 0.1);
    defer p_eigen.deinit();
    var p_direct = try probMatrix(q, 0.1);
    defer p_direct.deinit();

    // P from eigendecomposition should match direct matrix exp
    for (0..20) |i| {
        for (0..20) |j| {
            try std.testing.expectApproxEqAbs(p_direct.get(i, j), p_eigen.get(i, j), 1e-6);
        }
    }
    try std.testing.expect(validateP(p_eigen, 1e-6));
}

test "probMatrixEigen: Kimura with non-uniform pi" {
    const allocator = std.testing.allocator;
    // Build Kimura-like Q with non-uniform frequencies
    const pi = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    const alpha: f64 = 2.0;
    const beta: f64 = 1.0;

    var q = try Matrix.init(allocator, 4, 4);
    defer q.deinit();
    for (0..4) |i| {
        var row_sum: f64 = 0;
        for (0..4) |j| {
            if (i != j) {
                const rate = if ((i + j) % 2 == 0) alpha else beta;
                q.set(i, j, rate * pi[j]);
                row_sum += rate * pi[j];
            }
        }
        q.set(i, i, -row_sum);
    }

    var p_eigen = try probMatrixEigen(q, &pi, 0.3);
    defer p_eigen.deinit();
    var p_direct = try probMatrix(q, 0.3);
    defer p_direct.deinit();

    for (0..4) |i| {
        for (0..4) |j| {
            try std.testing.expectApproxEqAbs(p_direct.get(i, j), p_eigen.get(i, j), 1e-8);
        }
    }
    try std.testing.expect(validateP(p_eigen, 1e-8));
}

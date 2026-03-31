// Dense matrix of f64 values with row-major storage.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const Matrix = struct {
    data: []f64,
    rows: usize,
    cols: usize,
    allocator: Allocator,

    /// Allocate a rows×cols matrix initialised to zero.
    pub fn init(allocator: Allocator, rows: usize, cols: usize) !Matrix {
        const data = try allocator.alloc(f64, rows * cols);
        @memset(data, 0.0);
        return Matrix{ .data = data, .rows = rows, .cols = cols, .allocator = allocator };
    }

    /// Allocate an n×n identity matrix.
    pub fn initIdentity(allocator: Allocator, n: usize) !Matrix {
        var m = try Matrix.init(allocator, n, n);
        for (0..n) |i| m.data[i * n + i] = 1.0;
        return m;
    }

    pub fn deinit(self: *Matrix) void {
        self.allocator.free(self.data);
    }

    pub fn get(self: Matrix, i: usize, j: usize) f64 {
        return self.data[i * self.cols + j];
    }

    pub fn set(self: *Matrix, i: usize, j: usize, val: f64) void {
        self.data[i * self.cols + j] = val;
    }

    /// Return a new matrix equal to self + other. Dimensions must match.
    pub fn add(self: Matrix, other: Matrix) !Matrix {
        if (self.rows != other.rows or self.cols != other.cols) return error.DimensionMismatch;
        const result = try Matrix.init(self.allocator, self.rows, self.cols);
        for (result.data, self.data, other.data) |*r, a, b| r.* = a + b;
        return result;
    }

    /// Scale all elements in-place by factor.
    pub fn scale(self: *Matrix, factor: f64) void {
        for (self.data) |*v| v.* *= factor;
    }

    /// Return a new matrix C = A * B.
    pub fn multiply(allocator: Allocator, a: Matrix, b: Matrix) !Matrix {
        if (a.cols != b.rows) return error.DimensionMismatch;
        var c = try Matrix.init(allocator, a.rows, b.cols);
        for (0..a.rows) |i| {
            for (0..b.cols) |j| {
                var dot: f64 = 0.0;
                for (0..a.cols) |k| dot += a.get(i, k) * b.get(k, j);
                c.set(i, j, dot);
            }
        }
        return c;
    }

    /// Return a new matrix equal to self transposed.
    pub fn transpose(self: Matrix) !Matrix {
        var t = try Matrix.init(self.allocator, self.cols, self.rows);
        for (0..self.rows) |i| {
            for (0..self.cols) |j| {
                t.set(j, i, self.get(i, j));
            }
        }
        return t;
    }

    /// Return a deep copy of this matrix (same allocator).
    pub fn clone(self: Matrix) !Matrix {
        const data = try self.allocator.dupe(f64, self.data);
        return Matrix{ .data = data, .rows = self.rows, .cols = self.cols, .allocator = self.allocator };
    }

    /// Return true when self is square and |a[i,j] - a[j,i]| <= tol for all i,j.
    pub fn isSymmetric(self: Matrix, tol: f64) bool {
        if (self.rows != self.cols) return false;
        const n = self.rows;
        for (0..n) |i| {
            for (i + 1..n) |j| {
                if (@abs(self.get(i, j) - self.get(j, i)) > tol) return false;
            }
        }
        return true;
    }

    /// Create an upper triangular matrix of size n.
    pub fn initUpper(allocator: Allocator, n: usize) !Matrix {
        return Matrix.init(allocator, n, n);
    }

    /// Sum of all elements.
    pub fn sum(self: Matrix) f64 {
        var s: f64 = 0;
        for (self.data) |v| s += v;
        return s;
    }

    /// Minimum element value.
    pub fn min(self: Matrix) f64 {
        var m: f64 = self.data[0];
        for (self.data[1..]) |v| {
            if (v < m) m = v;
        }
        return m;
    }

    /// Maximum element value.
    pub fn max(self: Matrix) f64 {
        var m: f64 = self.data[0];
        for (self.data[1..]) |v| {
            if (v > m) m = v;
        }
        return m;
    }

    /// Frobenius norm: sqrt(sum of squared elements).
    pub fn frobeniusNorm(self: Matrix) f64 {
        var s: f64 = 0;
        for (self.data) |v| s += v * v;
        return @sqrt(s);
    }

    /// LU decomposition with partial pivoting (Doolittle).
    /// Returns L, U matrices and pivot permutation array.
    /// A = P * L * U where pivot[i] gives the row that was swapped to position i.
    pub fn luDecompose(self: Matrix) !struct { l: Matrix, u: Matrix, pivot: []usize } {
        if (self.rows != self.cols) return error.DimensionMismatch;
        const n = self.rows;

        // Work on a copy
        var a = try self.clone();
        errdefer a.deinit();

        const pivot = try self.allocator.alloc(usize, n);
        errdefer self.allocator.free(pivot);
        for (0..n) |i| pivot[i] = i;

        // Gaussian elimination with partial pivoting
        for (0..n) |k| {
            // Find pivot
            var max_val: f64 = @abs(a.get(k, k));
            var max_row: usize = k;
            for (k + 1..n) |i| {
                const v = @abs(a.get(i, k));
                if (v > max_val) {
                    max_val = v;
                    max_row = i;
                }
            }

            // Swap rows
            if (max_row != k) {
                const tmp = pivot[k];
                pivot[k] = pivot[max_row];
                pivot[max_row] = tmp;
                for (0..n) |j| {
                    const t = a.get(k, j);
                    a.set(k, j, a.get(max_row, j));
                    a.set(max_row, j, t);
                }
            }

            const akk = a.get(k, k);
            if (@abs(akk) < 1e-15) return error.SingularMatrix;

            // Eliminate below
            for (k + 1..n) |i| {
                const factor = a.get(i, k) / akk;
                a.set(i, k, factor); // Store L factor in lower triangle
                for (k + 1..n) |j| {
                    a.set(i, j, a.get(i, j) - factor * a.get(k, j));
                }
            }
        }

        // Extract L and U
        var l = try Matrix.initIdentity(self.allocator, n);
        errdefer l.deinit();
        var u = try Matrix.init(self.allocator, n, n);
        errdefer u.deinit();

        for (0..n) |i| {
            for (0..n) |j| {
                if (i > j) {
                    l.set(i, j, a.get(i, j));
                } else {
                    u.set(i, j, a.get(i, j));
                }
            }
        }

        a.deinit();
        return .{ .l = l, .u = u, .pivot = pivot };
    }

    /// Invert a square matrix using LU decomposition.
    /// Returns a new matrix. Caller owns the result.
    pub fn invert(self: Matrix) !Matrix {
        if (self.rows != self.cols) return error.DimensionMismatch;
        const n = self.rows;

        var lu = try self.luDecompose();
        defer {
            lu.l.deinit();
            lu.u.deinit();
            self.allocator.free(lu.pivot);
        }

        var inv = try Matrix.init(self.allocator, n, n);
        errdefer inv.deinit();

        // Solve L*U*x = e_j for each column j, using the pivot permutation
        var col = try self.allocator.alloc(f64, n);
        defer self.allocator.free(col);

        for (0..n) |j| {
            // Set up permuted identity column
            @memset(col, 0);
            col[lu.pivot[j]] = 1.0;

            // Apply pivot to get the right-hand side
            var b = try self.allocator.alloc(f64, n);
            defer self.allocator.free(b);
            for (0..n) |i| b[i] = 0;
            b[j] = 1.0;

            // Permute b
            var pb = try self.allocator.alloc(f64, n);
            defer self.allocator.free(pb);
            for (0..n) |i| pb[i] = b[lu.pivot[i]];

            // Forward substitution: L * y = Pb
            var y = try self.allocator.alloc(f64, n);
            defer self.allocator.free(y);
            for (0..n) |i| {
                var s: f64 = pb[i];
                for (0..i) |k| s -= lu.l.get(i, k) * y[k];
                y[i] = s;
            }

            // Back substitution: U * x = y
            var ii: usize = n;
            while (ii > 0) {
                ii -= 1;
                var s: f64 = y[ii];
                for (ii + 1..n) |k| s -= lu.u.get(ii, k) * col[k];
                col[ii] = s / lu.u.get(ii, ii);
            }

            for (0..n) |i| inv.set(i, j, col[i]);
        }

        return inv;
    }

    /// Matrix exponentiation: compute exp(t * M) using scaling and squaring
    /// with Padé approximation. M must be square.
    /// Returns a new matrix. Caller owns the result.
    pub fn exp(self: Matrix, t: f64) !Matrix {
        if (self.rows != self.cols) return error.DimensionMismatch;
        const n = self.rows;

        // Scale: find s such that ||t*M/2^s|| < 1
        var norm: f64 = 0;
        for (self.data) |v| {
            const av = @abs(t * v);
            if (av > norm) norm = av;
        }

        var s: u32 = 0;
        while (norm > 0.5) {
            norm /= 2.0;
            s += 1;
        }

        const sfactor = t / @as(f64, @floatFromInt(@as(u64, 1) << @intCast(s)));

        // Scaled matrix A = sfactor * M
        var a_mat = try self.clone();
        defer a_mat.deinit();
        a_mat.scale(sfactor);

        // Padé(6,6) approximation: exp(A) ≈ (D + N) * (D - N)^(-1)
        // where N = sum_{k=1}^{q} c_k * A^k, D = sum_{k=0}^{q} (-1)^k * c_k * A^k
        // Using Taylor series truncated at order 12 for simplicity:
        // exp(A) ≈ I + A + A^2/2 + A^3/6 + ... + A^12/12!

        var result = try Matrix.initIdentity(self.allocator, n);
        errdefer result.deinit();

        var term = try Matrix.initIdentity(self.allocator, n);
        defer term.deinit();

        // Compute Taylor series: exp(A) = I + A + A^2/2! + ... + A^12/12!
        for (1..13) |k| {
            var new_term = try Matrix.multiply(self.allocator, term, a_mat);
            new_term.scale(1.0 / @as(f64, @floatFromInt(k)));
            term.deinit();
            term = new_term;

            // Add term to result
            for (0..n * n) |idx| {
                result.data[idx] += term.data[idx];
            }
        }

        // Square s times: exp(t*M) = (exp(t*M/2^s))^(2^s)
        for (0..s) |_| {
            const squared = try Matrix.multiply(self.allocator, result, result);
            result.deinit();
            result = squared;
        }

        return result;
    }
};

// --- Tests ---

test "Matrix.init: zero-initialised" {
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 2, 3);
    defer m.deinit();
    for (m.data) |v| try std.testing.expectEqual(@as(f64, 0.0), v);
    try std.testing.expectEqual(@as(usize, 2), m.rows);
    try std.testing.expectEqual(@as(usize, 3), m.cols);
}

test "Matrix.initIdentity" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 3);
    defer m.deinit();
    for (0..3) |i| {
        for (0..3) |j| {
            const expected: f64 = if (i == j) 1.0 else 0.0;
            try std.testing.expectEqual(expected, m.get(i, j));
        }
    }
}

test "Matrix.get and set" {
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 2, 2);
    defer m.deinit();
    m.set(0, 1, 3.14);
    try std.testing.expectApproxEqAbs(@as(f64, 3.14), m.get(0, 1), 1e-10);
    try std.testing.expectEqual(@as(f64, 0.0), m.get(1, 0));
}

test "Matrix.add" {
    const allocator = std.testing.allocator;
    var a = try Matrix.init(allocator, 2, 2);
    defer a.deinit();
    var b = try Matrix.init(allocator, 2, 2);
    defer b.deinit();
    a.set(0, 0, 1.0);
    a.set(1, 1, 2.0);
    b.set(0, 0, 3.0);
    b.set(0, 1, 4.0);
    var c = try a.add(b);
    defer c.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), c.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), c.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), c.get(1, 1), 1e-10);
}

test "Matrix.scale in-place" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 2);
    defer m.deinit();
    m.scale(3.0);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), m.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), m.get(0, 1), 1e-10);
}

test "Matrix.multiply: 2x2 * 2x2" {
    const allocator = std.testing.allocator;
    // A = [[1,2],[3,4]]  B = [[5,6],[7,8]]
    // C = [[1*5+2*7, 1*6+2*8],[3*5+4*7, 3*6+4*8]] = [[19,22],[43,50]]
    var a = try Matrix.init(allocator, 2, 2);
    defer a.deinit();
    a.set(0, 0, 1);
    a.set(0, 1, 2);
    a.set(1, 0, 3);
    a.set(1, 1, 4);
    var b = try Matrix.init(allocator, 2, 2);
    defer b.deinit();
    b.set(0, 0, 5);
    b.set(0, 1, 6);
    b.set(1, 0, 7);
    b.set(1, 1, 8);
    var c = try Matrix.multiply(allocator, a, b);
    defer c.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 19.0), c.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 22.0), c.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 43.0), c.get(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 50.0), c.get(1, 1), 1e-10);
}

test "Matrix.multiply: dimension mismatch" {
    const allocator = std.testing.allocator;
    var a = try Matrix.init(allocator, 2, 3);
    defer a.deinit();
    var b = try Matrix.init(allocator, 2, 2);
    defer b.deinit();
    try std.testing.expectError(error.DimensionMismatch, Matrix.multiply(allocator, a, b));
}

test "Matrix.transpose" {
    const allocator = std.testing.allocator;
    // [[1,2,3],[4,5,6]] -> [[1,4],[2,5],[3,6]]
    var m = try Matrix.init(allocator, 2, 3);
    defer m.deinit();
    m.set(0, 0, 1);
    m.set(0, 1, 2);
    m.set(0, 2, 3);
    m.set(1, 0, 4);
    m.set(1, 1, 5);
    m.set(1, 2, 6);
    var t = try m.transpose();
    defer t.deinit();
    try std.testing.expectEqual(@as(usize, 3), t.rows);
    try std.testing.expectEqual(@as(usize, 2), t.cols);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), t.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), t.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), t.get(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), t.get(2, 1), 1e-10);
}

test "Matrix.clone is independent" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 2);
    defer m.deinit();
    var c = try m.clone();
    defer c.deinit();
    c.set(0, 0, 99.0);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), m.get(0, 0), 1e-10);
}

test "Matrix.sum" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 3);
    defer m.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), m.sum(), 1e-10);
}

test "Matrix.min and max" {
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 2, 2);
    defer m.deinit();
    m.set(0, 0, -5.0);
    m.set(0, 1, 3.0);
    m.set(1, 0, 0.0);
    m.set(1, 1, 7.0);
    try std.testing.expectApproxEqAbs(@as(f64, -5.0), m.min(), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 7.0), m.max(), 1e-10);
}

test "Matrix.frobeniusNorm: identity" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 3);
    defer m.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, @sqrt(3.0)), m.frobeniusNorm(), 1e-10);
}

test "Matrix.luDecompose: 2x2" {
    const allocator = std.testing.allocator;
    var a = try Matrix.init(allocator, 2, 2);
    defer a.deinit();
    a.set(0, 0, 4.0);
    a.set(0, 1, 3.0);
    a.set(1, 0, 6.0);
    a.set(1, 1, 3.0);

    var lu = try a.luDecompose();
    defer {
        lu.l.deinit();
        lu.u.deinit();
        allocator.free(lu.pivot);
    }

    // Verify L*U reconstructs A (with permutation)
    var prod = try Matrix.multiply(allocator, lu.l, lu.u);
    defer prod.deinit();

    // The product should match the permuted A
    for (0..2) |i| {
        for (0..2) |j| {
            try std.testing.expectApproxEqAbs(a.get(lu.pivot[i], j), prod.get(i, j), 1e-10);
        }
    }
}

test "Matrix.invert: 2x2" {
    const allocator = std.testing.allocator;
    var a = try Matrix.init(allocator, 2, 2);
    defer a.deinit();
    // [[1, 2], [3, 4]] -> inv = [[-2, 1], [1.5, -0.5]]
    a.set(0, 0, 1.0);
    a.set(0, 1, 2.0);
    a.set(1, 0, 3.0);
    a.set(1, 1, 4.0);

    var inv = try a.invert();
    defer inv.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, -2.0), inv.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), inv.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), inv.get(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, -0.5), inv.get(1, 1), 1e-10);
}

test "Matrix.invert: A * A^-1 = I" {
    const allocator = std.testing.allocator;
    var a = try Matrix.init(allocator, 3, 3);
    defer a.deinit();
    a.set(0, 0, 2);
    a.set(0, 1, 1);
    a.set(0, 2, 1);
    a.set(1, 0, 4);
    a.set(1, 1, 3);
    a.set(1, 2, 3);
    a.set(2, 0, 8);
    a.set(2, 1, 7);
    a.set(2, 2, 9);

    var inv = try a.invert();
    defer inv.deinit();
    var prod = try Matrix.multiply(allocator, a, inv);
    defer prod.deinit();

    // Should be identity
    for (0..3) |i| {
        for (0..3) |j| {
            const expected: f64 = if (i == j) 1.0 else 0.0;
            try std.testing.expectApproxEqAbs(expected, prod.get(i, j), 1e-10);
        }
    }
}

test "Matrix.exp: identity * t = exp(0) = I" {
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 2, 2);
    defer m.deinit();
    // exp(0 * M) = I
    var result = try m.exp(0.0);
    defer result.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.get(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result.get(1, 1), 1e-10);
}

test "Matrix.exp: diagonal matrix" {
    const allocator = std.testing.allocator;
    // exp(t*diag(a,b)) = diag(exp(ta), exp(tb))
    var m = try Matrix.init(allocator, 2, 2);
    defer m.deinit();
    m.set(0, 0, 1.0);
    m.set(1, 1, 2.0);
    var result = try m.exp(1.0);
    defer result.deinit();
    try std.testing.expectApproxEqAbs(@exp(@as(f64, 1.0)), result.get(0, 0), 1e-8);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.get(0, 1), 1e-8);
    try std.testing.expectApproxEqAbs(@exp(@as(f64, 2.0)), result.get(1, 1), 1e-8);
}

test "Matrix.exp: known 2x2 case" {
    const allocator = std.testing.allocator;
    // exp([[0, 1], [0, 0]]) = [[1, 1], [0, 1]]
    var m = try Matrix.init(allocator, 2, 2);
    defer m.deinit();
    m.set(0, 1, 1.0);
    var result = try m.exp(1.0);
    defer result.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result.get(0, 0), 1e-8);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result.get(0, 1), 1e-8);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result.get(1, 0), 1e-8);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result.get(1, 1), 1e-8);
}

test "Matrix.isSymmetric" {
    const allocator = std.testing.allocator;
    var sym = try Matrix.initIdentity(allocator, 3);
    defer sym.deinit();
    sym.set(0, 1, 2.0);
    sym.set(1, 0, 2.0);
    try std.testing.expect(sym.isSymmetric(1e-10));

    var asym = try Matrix.initIdentity(allocator, 2);
    defer asym.deinit();
    asym.set(0, 1, 1.0);
    // asym[1,0] stays 0, so not symmetric
    try std.testing.expect(!asym.isSymmetric(1e-10));

    // Non-square is not symmetric
    var rect = try Matrix.init(allocator, 2, 3);
    defer rect.deinit();
    try std.testing.expect(!rect.isSymmetric(1e-10));
}

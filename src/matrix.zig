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

    /// Return a new matrix with all elements scaled by factor.
    pub fn scale(self: Matrix, factor: f64) !Matrix {
        const result = try self.allocator.dupe(f64, self.data);
        for (result) |*v| v.* *= factor;
        return Matrix{ .data = result, .rows = self.rows, .cols = self.cols, .allocator = self.allocator };
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
        var a_mat = try self.scale(sfactor);
        defer a_mat.deinit();

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
            var product = try Matrix.multiply(self.allocator, term, a_mat);
            term.deinit();
            term = try product.scale(1.0 / @as(f64, @floatFromInt(k)));
            product.deinit();

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

    /// Result type for eigenvalue decomposition.
    pub const EigenResult = struct {
        eigenvalues_re: []f64,
        eigenvalues_im: []f64,
        eigenvectors: Matrix,
        allocator: Allocator,

        pub fn deinit(self: *EigenResult) void {
            self.allocator.free(self.eigenvalues_re);
            self.allocator.free(self.eigenvalues_im);
            self.eigenvectors.deinit();
        }
    };

    /// Eigenvalue decomposition for a square symmetric matrix using the
    /// classical Jacobi eigenvalue algorithm.
    ///
    /// The Jacobi method iteratively applies Givens rotations to zero
    /// off-diagonal elements until the matrix is approximately diagonal.
    /// It is robust and well-suited for the small dense matrices (4x4,
    /// 20x20) typical in biological sequence analysis (rate matrices).
    ///
    /// The matrix must be symmetric; returns error.NotSymmetric otherwise.
    /// Returns eigenvalues (all real for symmetric matrices, so
    /// eigenvalues_im is all zeros) and the orthogonal eigenvector matrix V
    /// such that A = V * diag(eigenvalues) * V^T.
    ///
    /// Modelled after Easel's esl_dmx_Diagonalize (which wraps LAPACK dgeev),
    /// but uses a pure-Zig Jacobi implementation instead of requiring LAPACK.
    pub fn diagonalize(self: Matrix, allocator: Allocator) !EigenResult {
        if (self.rows != self.cols) return error.DimensionMismatch;
        if (!self.isSymmetric(1e-12)) return error.NotSymmetric;
        const n = self.rows;

        // Work on a copy that will converge to a diagonal matrix
        var a = try self.clone();
        defer a.deinit();

        // V accumulates the product of all Givens rotations (starts as identity)
        var v = try Matrix.initIdentity(allocator, n);
        errdefer v.deinit();

        const max_iterations: usize = 100 * n * n;
        var iteration: usize = 0;

        while (iteration < max_iterations) : (iteration += 1) {
            // Find the largest off-diagonal element
            var max_off: f64 = 0.0;
            var p: usize = 0;
            var q: usize = 1;
            for (0..n) |i| {
                for (i + 1..n) |j| {
                    const aij = @abs(a.get(i, j));
                    if (aij > max_off) {
                        max_off = aij;
                        p = i;
                        q = j;
                    }
                }
            }

            // Convergence check: off-diagonal elements are negligible
            if (max_off < 1e-15) break;

            // Compute the Jacobi rotation angle
            const app = a.get(p, p);
            const aqq = a.get(q, q);
            const apq = a.get(p, q);

            var cos: f64 = undefined;
            var sin: f64 = undefined;

            if (@abs(app - aqq) < 1e-15) {
                // theta = pi/4
                const inv_sqrt2 = 1.0 / @sqrt(2.0);
                cos = inv_sqrt2;
                sin = inv_sqrt2;
            } else {
                const tau = (aqq - app) / (2.0 * apq);
                // t = sign(tau) / (|tau| + sqrt(1 + tau^2))
                const t_val = if (tau >= 0)
                    1.0 / (tau + @sqrt(1.0 + tau * tau))
                else
                    -1.0 / (-tau + @sqrt(1.0 + tau * tau));
                cos = 1.0 / @sqrt(1.0 + t_val * t_val);
                sin = t_val * cos;
            }

            // Apply rotation to A: A' = G^T * A * G
            // Update rows/columns p and q
            for (0..n) |r| {
                if (r == p or r == q) continue;
                const arp = a.get(r, p);
                const arq = a.get(r, q);
                const new_rp = cos * arp - sin * arq;
                const new_rq = sin * arp + cos * arq;
                a.set(r, p, new_rp);
                a.set(p, r, new_rp);
                a.set(r, q, new_rq);
                a.set(q, r, new_rq);
            }

            // Update diagonal and off-diagonal for p,q
            const new_pp = cos * cos * app - 2.0 * sin * cos * apq + sin * sin * aqq;
            const new_qq = sin * sin * app + 2.0 * sin * cos * apq + cos * cos * aqq;
            a.set(p, p, new_pp);
            a.set(q, q, new_qq);
            a.set(p, q, 0.0);
            a.set(q, p, 0.0);

            // Accumulate eigenvectors: V' = V * G
            for (0..n) |i| {
                const vip = v.get(i, p);
                const viq = v.get(i, q);
                v.set(i, p, cos * vip - sin * viq);
                v.set(i, q, sin * vip + cos * viq);
            }
        }

        // Extract eigenvalues from the diagonal
        const eigenvalues_re = try allocator.alloc(f64, n);
        errdefer allocator.free(eigenvalues_re);
        const eigenvalues_im = try allocator.alloc(f64, n);
        errdefer allocator.free(eigenvalues_im);

        for (0..n) |i| {
            eigenvalues_re[i] = a.get(i, i);
            eigenvalues_im[i] = 0.0;
        }

        return EigenResult{
            .eigenvalues_re = eigenvalues_re,
            .eigenvalues_im = eigenvalues_im,
            .eigenvectors = v,
            .allocator = allocator,
        };
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

test "Matrix.scale returns new matrix" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 2);
    defer m.deinit();
    var s = try m.scale(3.0);
    defer s.deinit();
    // Scaled matrix has correct values.
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), s.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), s.get(0, 1), 1e-10);
    // Original is unchanged.
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), m.get(0, 0), 1e-10);
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

test "Matrix.diagonalize: 2x2 symmetric" {
    const allocator = std.testing.allocator;
    // A = [[2, 1], [1, 2]]
    // Eigenvalues: 3 and 1
    var a = try Matrix.init(allocator, 2, 2);
    defer a.deinit();
    a.set(0, 0, 2.0);
    a.set(0, 1, 1.0);
    a.set(1, 0, 1.0);
    a.set(1, 1, 2.0);

    var eigen = try a.diagonalize(allocator);
    defer eigen.deinit();

    // Sort eigenvalues for deterministic comparison
    const e0 = @min(eigen.eigenvalues_re[0], eigen.eigenvalues_re[1]);
    const e1 = @max(eigen.eigenvalues_re[0], eigen.eigenvalues_re[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), e0, 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), e1, 1e-10);

    // Imaginary parts should be zero
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), eigen.eigenvalues_im[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), eigen.eigenvalues_im[1], 1e-10);
}

test "Matrix.diagonalize: identity" {
    const allocator = std.testing.allocator;
    var a = try Matrix.initIdentity(allocator, 3);
    defer a.deinit();

    var eigen = try a.diagonalize(allocator);
    defer eigen.deinit();

    for (0..3) |i| {
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), eigen.eigenvalues_re[i], 1e-10);
        try std.testing.expectApproxEqAbs(@as(f64, 0.0), eigen.eigenvalues_im[i], 1e-10);
    }
}

test "Matrix.diagonalize: V * D * V^T reconstructs A" {
    const allocator = std.testing.allocator;
    // 3x3 symmetric matrix
    var a = try Matrix.init(allocator, 3, 3);
    defer a.deinit();
    a.set(0, 0, 4.0);
    a.set(0, 1, 1.0);
    a.set(0, 2, 0.5);
    a.set(1, 0, 1.0);
    a.set(1, 1, 3.0);
    a.set(1, 2, 0.25);
    a.set(2, 0, 0.5);
    a.set(2, 1, 0.25);
    a.set(2, 2, 2.0);

    var eigen = try a.diagonalize(allocator);
    defer eigen.deinit();

    // Reconstruct: A_reconstructed = V * diag(eigenvalues) * V^T
    var d = try Matrix.init(allocator, 3, 3);
    defer d.deinit();
    for (0..3) |i| d.set(i, i, eigen.eigenvalues_re[i]);

    var vd = try Matrix.multiply(allocator, eigen.eigenvectors, d);
    defer vd.deinit();

    var vt = try eigen.eigenvectors.transpose();
    defer vt.deinit();

    var reconstructed = try Matrix.multiply(allocator, vd, vt);
    defer reconstructed.deinit();

    // Check A == V * D * V^T
    for (0..3) |i| {
        for (0..3) |j| {
            try std.testing.expectApproxEqAbs(a.get(i, j), reconstructed.get(i, j), 1e-10);
        }
    }
}

test "Matrix.diagonalize: non-symmetric returns error" {
    const allocator = std.testing.allocator;
    var a = try Matrix.init(allocator, 2, 2);
    defer a.deinit();
    a.set(0, 0, 1.0);
    a.set(0, 1, 2.0);
    a.set(1, 0, 3.0); // asymmetric
    a.set(1, 1, 4.0);

    try std.testing.expectError(error.NotSymmetric, a.diagonalize(allocator));
}

test "Matrix.diagonalize: 4x4 rate-matrix-like" {
    const allocator = std.testing.allocator;
    // Symmetric 4x4 (typical size for DNA rate matrices after symmetrization)
    var a = try Matrix.init(allocator, 4, 4);
    defer a.deinit();
    const vals = [_]f64{
        5.0, 1.0, 0.5, 0.2,
        1.0, 4.0, 0.8, 0.3,
        0.5, 0.8, 3.0, 0.6,
        0.2, 0.3, 0.6, 2.0,
    };
    for (0..16) |idx| a.data[idx] = vals[idx];

    var eigen = try a.diagonalize(allocator);
    defer eigen.deinit();

    // Verify orthogonality of eigenvectors: V^T * V = I
    var vt = try eigen.eigenvectors.transpose();
    defer vt.deinit();
    var vtv = try Matrix.multiply(allocator, vt, eigen.eigenvectors);
    defer vtv.deinit();

    for (0..4) |i| {
        for (0..4) |j| {
            const expected: f64 = if (i == j) 1.0 else 0.0;
            try std.testing.expectApproxEqAbs(expected, vtv.get(i, j), 1e-10);
        }
    }
}

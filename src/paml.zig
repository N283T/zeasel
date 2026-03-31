// PAML amino acid exchangeability matrix reader.
//
// Reads the lower-triangular exchangeability matrix format used by
// PAML (Yang, 1997) and many evolutionary models (WAG, LG, JTT).
// Format: 190 exchangeability values (lower triangle of 20x20 matrix)
// followed by 20 stationary frequencies.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Matrix = @import("matrix.zig").Matrix;

/// Result of parsing a PAML rate matrix file.
pub const PamlModel = struct {
    /// Symmetric 20x20 exchangeability matrix S (S[i][j] = S[j][i]).
    exchangeability: Matrix,
    /// Stationary frequencies pi[0..19].
    frequencies: [20]f64,
    allocator: Allocator,

    pub fn deinit(self: *PamlModel) void {
        self.exchangeability.deinit();
    }
};

// PAML files use alphabetical-by-3-letter-code order: ARNDCQEGHILKMFPSTWYV
// zeasel (like Easel) uses alphabetical-by-1-letter-code: ACDEFGHIKLMNPQRSTVWY
// This table maps PAML index -> zeasel index.
const paml_to_zeasel: [20]usize = buildPermutation();

fn buildPermutation() [20]usize {
    const paml_order = "ARNDCQEGHILKMFPSTWYV";
    const zeasel_order = "ACDEFGHIKLMNPQRSTVWY";
    var perm: [20]usize = undefined;
    for (paml_order, 0..) |c, i| {
        for (zeasel_order, 0..) |z, j| {
            if (c == z) {
                perm[i] = j;
                break;
            }
        }
    }
    return perm;
}

/// Parse a PAML-format amino acid model from text.
/// Expects 190 exchangeability values (lower triangle, row by row)
/// followed by 20 frequency values.
/// Values are permuted from PAML order (ARNDCQEGHILKMFPSTWYV) to
/// zeasel alphabet order (ACDEFGHIKLMNPQRSTVWY).
pub fn parse(allocator: Allocator, data: []const u8) !PamlModel {
    var s = try Matrix.init(allocator, 20, 20);
    errdefer s.deinit();

    var it = std.mem.tokenizeAny(u8, data, " \t\n\r");

    // Read 190 lower-triangle values: S[i][j] for i=1..19, j=0..i-1
    // Permute from PAML order to zeasel order as we read.
    for (1..20) |i| {
        for (0..i) |j| {
            const token = it.next() orelse return error.InvalidFormat;
            const val = std.fmt.parseFloat(f64, token) catch return error.InvalidFormat;
            const pi = paml_to_zeasel[i];
            const pj = paml_to_zeasel[j];
            s.set(pi, pj, val);
            s.set(pj, pi, val); // symmetric
        }
    }

    // Read 20 frequencies, permuting from PAML order to zeasel order.
    var freq: [20]f64 = undefined;
    for (0..20) |i| {
        const token = it.next() orelse return error.InvalidFormat;
        freq[paml_to_zeasel[i]] = std.fmt.parseFloat(f64, token) catch return error.InvalidFormat;
    }

    return PamlModel{
        .exchangeability = s,
        .frequencies = freq,
        .allocator = allocator,
    };
}

// --- Tests ---

test "parse: minimal 4-residue-like format (using 20x20)" {
    const allocator = std.testing.allocator;

    // Build a valid 190 + 20 = 210 token input
    var buf: [4096]u8 = undefined;
    var stream = std.io.fixedBufferStream(&buf);
    const writer = stream.writer();

    // 190 exchangeability values
    for (0..190) |i| {
        try writer.print("{d:.1} ", .{@as(f64, @floatFromInt(i + 1))});
    }
    // 20 frequencies (uniform)
    for (0..20) |_| {
        try writer.print("0.05 ", .{});
    }

    const input = buf[0..stream.pos];

    var model = try parse(allocator, input);
    defer model.deinit();

    // S should be symmetric
    try std.testing.expect(model.exchangeability.isSymmetric(1e-10));
    // First value (1.0) is S_paml[1][0] = S(R,A).
    // PAML R -> zeasel index 14, PAML A -> zeasel index 0.
    // So S_zeasel[14][0] should be 1.0.
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), model.exchangeability.get(paml_to_zeasel[1], paml_to_zeasel[0]), 1e-10);
    // Should be symmetric
    try std.testing.expectApproxEqAbs(
        model.exchangeability.get(paml_to_zeasel[1], paml_to_zeasel[0]),
        model.exchangeability.get(paml_to_zeasel[0], paml_to_zeasel[1]),
        1e-10,
    );
    // Frequencies should sum to 1
    var sum: f64 = 0;
    for (model.frequencies) |f| sum += f;
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), sum, 1e-10);
}

test "parse: too few tokens returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(error.InvalidFormat, parse(allocator, "1.0 2.0 3.0"));
}

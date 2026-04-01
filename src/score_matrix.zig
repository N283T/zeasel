// Substitution score matrices for biological sequences.
//
// Canonical residue order (amino acid digital codes 0-19):
//   A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
//   M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19
//
// The full matrix is Kp x Kp (29 x 29 for amino acids), covering
// canonical residues, gap, degenerate codes (B, J, Z, O, U, X),
// nonresidue (*), and missing data (~).

const std = @import("std");
const math = std.math;
const alphabet_mod = @import("alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;
const composition = @import("composition.zig");
const Matrix = @import("matrix.zig").Matrix;
const Allocator = std.mem.Allocator;

/// Size of the full score matrix dimension (amino acid Kp).
pub const Kp: usize = 29;

/// Integer division with rounding to nearest (half away from zero).
fn roundedDiv(num: i32, den: i32) i32 {
    if (num >= 0) {
        return @divTrunc(num + @divTrunc(den, 2), den);
    } else {
        return @divTrunc(num - @divTrunc(den, 2), den);
    }
}

pub const ScoreMatrix = struct {
    /// Raw score table: data[i][j] for residue digital codes i, j (0..Kp-1).
    /// Positions beyond canonical residues hold degenerate/gap/special scores.
    data: *const [Kp][Kp]i16,
    /// Human-readable name (e.g. "BLOSUM62").
    name: []const u8,
    /// Associated alphabet.
    abc: *const Alphabet,
    /// If non-null, this matrix owns its data (allocated at runtime).
    owned: ?*[Kp][Kp]i16 = null,
    /// Allocator used for owned data (null for comptime matrices).
    allocator: ?Allocator = null,
    /// Allocated name storage (null for comptime matrices).
    owned_name: ?[]u8 = null,

    /// Return score for two digital codes.
    pub fn score(self: ScoreMatrix, i: u8, j: u8) i16 {
        std.debug.assert(i < Kp and j < Kp);
        return self.data[i][j];
    }

    /// Return score for two ASCII residue characters.
    /// Returns error.InvalidCharacter if either character is not recognized.
    pub fn scoreByChar(self: ScoreMatrix, a: u8, b: u8) !i16 {
        const i = try self.abc.encode(a);
        const j = try self.abc.encode(b);
        if (i >= Kp or j >= Kp) return error.InvalidCharacter;
        return self.data[i][j];
    }

    /// Check if the matrix is symmetric (S[i][j] == S[j][i] for all i,j in 0..K-1).
    pub fn isSymmetric(self: ScoreMatrix) bool {
        const k = self.abc.k;
        for (0..k) |i| {
            for (i + 1..k) |j| {
                if (self.data[i][j] != self.data[j][i]) return false;
            }
        }
        return true;
    }

    /// Return the maximum score value among canonical residue pairs.
    pub fn maxScore(self: ScoreMatrix) i16 {
        const k = self.abc.k;
        var m: i16 = self.data[0][0];
        for (0..k) |i| {
            for (0..k) |j| {
                if (self.data[i][j] > m) m = self.data[i][j];
            }
        }
        return m;
    }

    /// Return the minimum score value among canonical residue pairs.
    pub fn minScore(self: ScoreMatrix) i16 {
        const k = self.abc.k;
        var m: i16 = self.data[0][0];
        for (0..k) |i| {
            for (0..k) |j| {
                if (self.data[i][j] < m) m = self.data[i][j];
            }
        }
        return m;
    }

    /// Compute the expected score E = sum_ij fi[i] * S_ij * fj[j].
    /// fi and fj are the query and target background frequencies respectively.
    /// For symmetric scoring, pass the same array for both.
    pub fn expectedScore(self: ScoreMatrix, fi: *const [20]f64, fj: *const [20]f64) f64 {
        var e: f64 = 0;
        for (0..20) |i| {
            for (0..20) |j| {
                e += fi[i] * @as(f64, @floatFromInt(self.data[i][j])) * fj[j];
            }
        }
        return e;
    }

    /// Compute relative entropy H = sum_ij p_ij * log2(p_ij / (f_i * f_j))
    /// given background frequencies, where p_ij = f_i * f_j * exp(lambda * s_ij).
    /// Returns the relative entropy in bits.
    /// Returns error.InvalidProbabilities if the background frequencies do not
    /// sum to approximately 1.0.
    pub fn relativeEntropy(self: ScoreMatrix, bg: *const [20]f64, lambda: f64) !f64 {
        // Validate that background probabilities sum to ~1.0
        var bg_sum: f64 = 0;
        for (bg) |p| bg_sum += p;
        if (@abs(bg_sum - 1.0) > 1e-6) return error.InvalidProbabilities;

        // Compute the normalizing constant Z = sum_ij f_i * f_j * exp(lambda * s_ij).
        // If Z deviates from 1.0 by more than the tolerance, the lambda/background
        // combination is inconsistent and we return an error rather than silently
        // normalizing (matching Easel's validation behavior).
        var z: f64 = 0;
        for (0..20) |i| {
            for (0..20) |j| {
                z += bg[i] * bg[j] * @exp(lambda * @as(f64, @floatFromInt(self.data[i][j])));
            }
        }
        if (@abs(z - 1.0) > 0.05) return error.InvalidNormalization;

        var h: f64 = 0;
        for (0..20) |i| {
            for (0..20) |j| {
                const s: f64 = @floatFromInt(self.data[i][j]);
                const p_ij = bg[i] * bg[j] * @exp(lambda * s);
                if (p_ij > 0) {
                    h += p_ij * lambda * s * math.log2e;
                }
            }
        }
        return h;
    }

    /// Fill in scores for degenerate residue positions (indices >= K)
    /// by averaging over the constituent canonical residue pairs using
    /// the alphabet degeneracy bitmasks.
    ///
    /// For a degenerate code ip that maps to canonical set {a1, a2, ...}
    /// and degenerate code jp that maps to {b1, b2, ...}, the score is:
    ///   S[ip][jp] = round( sum(S[ai][bj]) / (n_ip * n_jp) )
    ///
    /// Gap (K), nonresidue (Kp-2), and missing (Kp-1) positions are set to 0.
    pub fn setDegenerateScores(self: *ScoreMatrix) void {
        const abc = self.abc;
        const k: usize = abc.k;
        const kp: usize = abc.kp;
        const ptr = self.owned orelse unreachable;

        // Fill canonical-to-degenerate: S[i][jp] for i < K, jp > K
        for (0..k) |i| {
            // Gap column
            ptr[i][k] = 0;
            // Degenerate columns: codes K+1 to Kp-4 (inclusive)
            var jp: usize = k + 1;
            while (jp + 2 < kp) : (jp += 1) {
                const mask_jp = abc.degen[jp];
                const n_jp: i32 = @intCast(abc.ndegen[jp]);
                if (n_jp == 0) {
                    ptr[i][jp] = 0;
                    continue;
                }
                var sum: i32 = 0;
                for (0..k) |j| {
                    if (mask_jp & (@as(u32, 1) << @intCast(j)) != 0) {
                        sum += ptr[i][j];
                    }
                }
                ptr[i][jp] = @intCast(roundedDiv(sum, n_jp));
            }
            // Nonresidue and missing columns
            ptr[i][kp - 2] = 0;
            ptr[i][kp - 1] = 0;
        }

        // Gap row: all zeros
        for (0..Kp) |j| {
            ptr[k][j] = 0;
        }

        // Degenerate rows: S[ip][j] for all j
        var ip: usize = k + 1;
        while (ip + 2 < kp) : (ip += 1) {
            const mask_ip = abc.degen[ip];
            const n_ip: i32 = @intCast(abc.ndegen[ip]);
            if (n_ip == 0) {
                for (0..Kp) |j| {
                    ptr[ip][j] = 0;
                }
                continue;
            }

            // Degenerate-to-canonical
            for (0..k) |j| {
                var sum: i32 = 0;
                for (0..k) |ci| {
                    if (mask_ip & (@as(u32, 1) << @intCast(ci)) != 0) {
                        sum += ptr[ci][j];
                    }
                }
                ptr[ip][j] = @intCast(roundedDiv(sum, n_ip));
            }

            // Gap
            ptr[ip][k] = 0;

            // Degenerate-to-degenerate
            var jp2: usize = k + 1;
            while (jp2 + 2 < kp) : (jp2 += 1) {
                const mask_jp = abc.degen[jp2];
                const n_jp: i32 = @intCast(abc.ndegen[jp2]);
                if (n_jp == 0) {
                    ptr[ip][jp2] = 0;
                    continue;
                }
                var sum: i32 = 0;
                for (0..k) |j| {
                    if (mask_jp & (@as(u32, 1) << @intCast(j)) != 0) {
                        sum += ptr[ip][j];
                    }
                }
                ptr[ip][jp2] = @intCast(roundedDiv(sum, n_jp));
            }

            // Nonresidue and missing
            ptr[ip][kp - 2] = 0;
            ptr[ip][kp - 1] = 0;
        }

        // Nonresidue and missing rows: all zeros
        for (0..Kp) |j| {
            ptr[kp - 2][j] = 0;
            ptr[kp - 1][j] = 0;
        }
    }

    /// Parse a BLAST/NCBI format score matrix from text data.
    ///
    /// Format: lines starting with '#' are comments. The first non-comment
    /// line is a header with single-letter residue labels. Subsequent rows
    /// optionally have a leading residue label, followed by integer scores.
    pub fn read(allocator: Allocator, abc: *const Alphabet, text: []const u8) !ScoreMatrix {
        const data_ptr = try allocator.create([Kp][Kp]i16);
        errdefer allocator.destroy(data_ptr);

        // Initialize to zero
        for (data_ptr) |*row| {
            @memset(row, 0);
        }

        // Split into lines
        var lines = std.mem.splitSequence(u8, text, "\n");

        // Skip comment lines and find header
        var header: ?[]const u8 = null;
        while (lines.next()) |line| {
            const trimmed = std.mem.trim(u8, line, " \t\r");
            if (trimmed.len == 0) continue;
            if (trimmed[0] == '#') continue;
            header = trimmed;
            break;
        }

        if (header == null) return error.InvalidFormat;

        // Parse header: single-character labels
        var col_map: [Kp]u8 = undefined; // maps column index -> digital code
        var nc: usize = 0;
        var hdr_iter = std.mem.tokenizeAny(u8, header.?, " \t");
        while (hdr_iter.next()) |tok| {
            if (tok.len != 1) return error.InvalidFormat;
            if (nc >= Kp) return error.InvalidFormat;
            const code = abc.encode(tok[0]) catch continue; // skip unknown chars like '*'
            if (code >= Kp) continue;
            col_map[nc] = code;
            nc += 1;
        }

        // Verify all canonical residues are present
        var have_canonical: [20]bool = .{false} ** 20;
        for (0..nc) |c| {
            if (col_map[c] < abc.k) {
                have_canonical[col_map[c]] = true;
            }
        }
        for (0..abc.k) |x| {
            if (!have_canonical[x]) return error.InvalidFormat;
        }

        // Parse data rows
        var row_idx: usize = 0;
        while (lines.next()) |line| {
            const trimmed = std.mem.trim(u8, line, " \t\r");
            if (trimmed.len == 0) continue;
            if (trimmed[0] == '#') continue;
            if (row_idx >= nc) return error.InvalidFormat;

            var col_idx: usize = 0;
            var tok_iter = std.mem.tokenizeAny(u8, trimmed, " \t");
            var skip_label = false;
            while (tok_iter.next()) |tok| {
                // Check if first token is a row label (single letter matching expected)
                if (col_idx == 0 and !skip_label and tok.len == 1) {
                    const maybe_code = abc.encode(tok[0]) catch null;
                    if (maybe_code) |code| {
                        if (code == col_map[row_idx]) {
                            skip_label = true;
                            continue;
                        }
                    }
                    // Not a label, try to parse as number
                }

                if (col_idx >= nc) return error.InvalidFormat;
                const val = std.fmt.parseInt(i16, tok, 10) catch return error.InvalidFormat;
                data_ptr[col_map[row_idx]][col_map[col_idx]] = val;
                col_idx += 1;
            }

            if (col_idx != nc) return error.InvalidFormat;
            row_idx += 1;
        }

        if (row_idx != nc) return error.InvalidFormat;

        return ScoreMatrix{
            .data = data_ptr,
            .name = "custom",
            .abc = abc,
            .owned = data_ptr,
            .allocator = allocator,
        };
    }

    /// Write the score matrix in BLAST/NCBI format.
    pub fn write(self: ScoreMatrix, dest: std.io.AnyWriter) !void {
        const k = self.abc.k;

        // Header line with column labels
        try dest.writeAll("  ");
        for (0..k) |j| {
            try dest.print("  {c} ", .{self.abc.decode(@intCast(j))});
        }
        try dest.writeByte('\n');

        // Data rows
        for (0..k) |i| {
            try dest.print("{c} ", .{self.abc.decode(@intCast(i))});
            for (0..k) |j| {
                try dest.print("{d:>3} ", .{self.data[i][j]});
            }
            try dest.writeByte('\n');
        }
    }

    /// Free runtime-allocated resources.
    pub fn deinit(self: *ScoreMatrix) void {
        if (self.allocator) |alloc| {
            if (self.owned) |ptr| {
                alloc.destroy(ptr);
                self.owned = null;
            }
            if (self.owned_name) |name_ptr| {
                alloc.free(name_ptr);
                self.owned_name = null;
            }
        }
    }
};

/// Look up a named score matrix. Returns null if name is not recognized.
pub fn getByName(name: []const u8) ?ScoreMatrix {
    if (std.ascii.eqlIgnoreCase(name, "BLOSUM62") or std.ascii.eqlIgnoreCase(name, "BL62")) return blosum62;
    if (std.ascii.eqlIgnoreCase(name, "BLOSUM45")) return blosum45;
    if (std.ascii.eqlIgnoreCase(name, "BLOSUM80")) return blosum80;
    if (std.ascii.eqlIgnoreCase(name, "PAM30")) return pam30;
    if (std.ascii.eqlIgnoreCase(name, "PAM70")) return pam70;
    if (std.ascii.eqlIgnoreCase(name, "PAM120")) return pam120;
    if (std.ascii.eqlIgnoreCase(name, "PAM250")) return pam250;
    if (std.ascii.eqlIgnoreCase(name, "IDENTITY")) return identity;
    return null;
}

// ---------------------------------------------------------------------------
// Helper: embed a 20x20 i8 matrix into a Kp x Kp i16 matrix at comptime
// ---------------------------------------------------------------------------

fn embedCanonical(comptime src: [20][20]i8) [Kp][Kp]i16 {
    var m: [Kp][Kp]i16 = .{.{0} ** Kp} ** Kp;
    var i: usize = 0;
    while (i < 20) : (i += 1) {
        var j: usize = 0;
        while (j < 20) : (j += 1) {
            m[i][j] = src[i][j];
        }
    }
    return m;
}

// ---------------------------------------------------------------------------
// BLOSUM62 data
// ---------------------------------------------------------------------------

const blosum62_canonical: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2 }, // A
    .{  0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2 }, // C
    .{ -2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3 }, // D
    .{ -1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2 }, // E
    .{ -2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3 }, // F
    .{  0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3 }, // G
    .{ -2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2 }, // H
    .{ -1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1 }, // I
    .{ -1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2 }, // K
    .{ -1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1 }, // L
    .{ -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1 }, // M
    .{ -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2 }, // N
    .{ -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3 }, // P
    .{ -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1 }, // Q
    .{ -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2 }, // R
    .{  1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2 }, // S
    .{  0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2 }, // T
    .{  0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1 }, // V
    .{ -3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2 }, // W
    .{ -2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7 }, // Y
};

const blosum62_data: [Kp][Kp]i16 = embedCanonical(blosum62_canonical);

pub const blosum62 = ScoreMatrix{
    .data = &blosum62_data,
    .name = "BLOSUM62",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// BLOSUM45 data
// ---------------------------------------------------------------------------

const blosum45_canonical: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  5, -1, -2, -1, -2,  0, -2, -1, -1, -1, -1, -1, -1, -1, -2,  1,  0,  0, -2, -2 }, // A
    .{ -1, 12, -3, -3, -2, -3, -3, -3, -3, -2, -2, -2, -4, -3, -3, -1, -1, -1, -5, -3 }, // C
    .{ -2, -3,  7,  2, -4, -1,  0, -4,  0, -3, -3,  2, -1,  0, -1,  0, -1, -3, -4, -2 }, // D
    .{ -1, -3,  2,  6, -3, -2,  0, -3,  1, -2, -2,  0,  0,  2,  0,  0, -1, -3, -3, -2 }, // E
    .{ -2, -2, -4, -3,  8, -3, -2,  0, -3,  1,  0, -2, -3, -4, -2, -2, -1,  0,  1,  3 }, // F
    .{  0, -3, -1, -2, -3,  7, -2, -4, -2, -3, -2,  0, -2, -2, -2,  0, -2, -3, -2, -3 }, // G
    .{ -2, -3,  0,  0, -2, -2, 10, -3, -1, -2,  0,  1, -2,  1,  0, -1, -2, -3, -3,  2 }, // H
    .{ -1, -3, -4, -3,  0, -4, -3,  5, -3,  2,  2, -2, -2, -2, -3, -2, -1,  3, -2,  0 }, // I
    .{ -1, -3,  0,  1, -3, -2, -1, -3,  5, -3, -1,  0, -1,  1,  3, -1, -1, -2, -2, -1 }, // K
    .{ -1, -2, -3, -2,  1, -3, -2,  2, -3,  5,  2, -3, -3, -2, -2, -3, -1,  1, -2,  0 }, // L
    .{ -1, -2, -3, -2,  0, -2,  0,  2, -1,  2,  6, -2, -2,  0, -1, -2, -1,  1, -2,  0 }, // M
    .{ -1, -2,  2,  0, -2,  0,  1, -2,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2 }, // N
    .{ -1, -4, -1,  0, -3, -2, -2, -2, -1, -3, -2, -2,  9, -1, -2, -1, -1, -3, -3, -3 }, // P
    .{ -1, -3,  0,  2, -4, -2,  1, -2,  1, -2,  0,  0, -1,  6,  1,  0, -1, -3, -2, -1 }, // Q
    .{ -2, -3, -1,  0, -2, -2,  0, -3,  3, -2, -1,  0, -2,  1,  7, -1, -1, -2, -2, -1 }, // R
    .{  1, -1,  0,  0, -2,  0, -1, -2, -1, -3, -2,  1, -1,  0, -1,  4,  2, -1, -4, -2 }, // S
    .{  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  2,  5,  0, -3, -1 }, // T
    .{  0, -1, -3, -3,  0, -3, -3,  3, -2,  1,  1, -3, -3, -3, -2, -1,  0,  5, -3, -1 }, // V
    .{ -2, -5, -4, -3,  1, -2, -3, -2, -2, -2, -2, -4, -3, -2, -2, -4, -3, -3, 15,  3 }, // W
    .{ -2, -3, -2, -2,  3, -3,  2,  0, -1,  0,  0, -2, -3, -1, -1, -2, -1, -1,  3,  8 }, // Y
};

const blosum45_data: [Kp][Kp]i16 = embedCanonical(blosum45_canonical);

pub const blosum45 = ScoreMatrix{
    .data = &blosum45_data,
    .name = "BLOSUM45",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// BLOSUM80 data
// ---------------------------------------------------------------------------

const blosum80_canonical: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  7, -1, -3, -2, -4,  0, -3, -3, -1, -3, -2, -3, -1, -2, -3,  2,  0, -1, -5, -4 }, // A
    .{ -1, 13, -7, -7, -4, -6, -7, -2, -6, -3, -3, -5, -6, -5, -6, -2, -2, -2, -5, -5 }, // C
    .{ -3, -7, 10,  2, -6, -3, -2, -5, -2, -7, -6,  2, -3, -1, -3, -1, -2, -6, -8, -6 }, // D
    .{ -2, -7,  2,  8, -6, -4,  0, -6,  1, -6, -4, -1, -2,  3, -1, -1, -2, -4, -6, -5 }, // E
    .{ -4, -4, -6, -6, 10, -6, -2, -1, -5,  2,  0, -6, -6, -5, -5, -4, -4, -2,  0,  4 }, // F
    .{  0, -6, -3, -4, -6,  9, -4, -7, -3, -7, -5, -1, -5, -4, -4, -1, -3, -6, -6, -6 }, // G
    .{ -3, -7, -2,  0, -2, -4, 12, -6, -1, -5, -4,  1, -4,  1,  0, -2, -3, -5, -4,  3 }, // H
    .{ -3, -2, -5, -6, -1, -7, -6,  7, -5,  2,  2, -6, -5, -5, -5, -4, -2,  4, -5, -3 }, // I
    .{ -1, -6, -2,  1, -5, -3, -1, -5,  8, -4, -3,  0, -2,  2,  3, -1, -1, -4, -6, -4 }, // K
    .{ -3, -3, -7, -6,  2, -7, -5,  2, -4,  6,  3, -6, -5, -4, -4, -4, -3,  1, -4, -2 }, // L
    .{ -2, -3, -6, -4,  0, -5, -4,  2, -3,  3,  9, -4, -4, -1, -3, -3, -1,  1, -3, -3 }, // M
    .{ -3, -5,  2, -1, -6, -1,  1, -6,  0, -6, -4,  9, -4,  0, -1,  1,  0, -5, -7, -4 }, // N
    .{ -1, -6, -3, -2, -6, -5, -4, -5, -2, -5, -4, -4, 12, -3, -3, -2, -3, -4, -7, -6 }, // P
    .{ -2, -5, -1,  3, -5, -4,  1, -5,  2, -4, -1,  0, -3,  9,  1, -1, -1, -4, -4, -3 }, // Q
    .{ -3, -6, -3, -1, -5, -4,  0, -5,  3, -4, -3, -1, -3,  1,  9, -2, -2, -4, -5, -4 }, // R
    .{  2, -2, -1, -1, -4, -1, -2, -4, -1, -4, -3,  1, -2, -1, -2,  7,  2, -3, -6, -3 }, // S
    .{  0, -2, -2, -2, -4, -3, -3, -2, -1, -3, -1,  0, -3, -1, -2,  2,  8, -1, -5, -3 }, // T
    .{ -1, -2, -6, -4, -2, -6, -5,  4, -4,  1,  1, -5, -4, -4, -4, -3, -1,  7, -5, -3 }, // V
    .{ -5, -5, -8, -6,  0, -6, -4, -5, -6, -4, -3, -7, -7, -4, -5, -6, -5, -5, 16,  3 }, // W
    .{ -4, -5, -6, -5,  4, -6,  3, -3, -4, -2, -3, -4, -6, -3, -4, -3, -3, -3,  3, 11 }, // Y
};

const blosum80_data: [Kp][Kp]i16 = embedCanonical(blosum80_canonical);

pub const blosum80 = ScoreMatrix{
    .data = &blosum80_data,
    .name = "BLOSUM80",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// PAM30 data
// ---------------------------------------------------------------------------

// Data from Easel esl_scorematrix.c (authoritative PAM30 values).
const pam30_canonical: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  6, -6, -3, -2, -8, -2, -7, -5, -7, -6, -5, -4, -2, -4, -7,  0, -1, -2,-13, -8 }, // A
    .{ -6, 10,-14,-14,-13, -9, -7, -6,-14,-15,-13,-11, -8,-14, -8, -3, -8, -6,-15, -4 }, // C
    .{ -3,-14,  8,  2,-15, -3, -4, -7, -4,-12,-11,  2, -8, -2,-10, -4, -5, -8,-15,-11 }, // D
    .{ -2,-14,  2,  8,-14, -4, -5, -5, -4, -9, -7, -2, -5,  1, -9, -4, -6, -6,-17, -8 }, // E
    .{ -8,-13,-15,-14,  9, -9, -6, -2,-14, -3, -4, -9,-10,-13, -9, -6, -9, -8, -4,  2 }, // F
    .{ -2, -9, -3, -4, -9,  6, -9,-11, -7,-10, -8, -3, -6, -7, -9, -2, -6, -5,-15,-14 }, // G
    .{ -7, -7, -4, -5, -6, -9,  9, -9, -6, -6,-10,  0, -4,  1, -2, -6, -7, -6, -7, -3 }, // H
    .{ -5, -6, -7, -5, -2,-11, -9,  8, -6, -1, -1, -5, -8, -8, -5, -7, -2,  2,-14, -6 }, // I
    .{ -7,-14, -4, -4,-14, -7, -6, -6,  7, -8, -2, -1, -6, -3,  0, -4, -3, -9,-12, -9 }, // K
    .{ -6,-15,-12, -9, -3,-10, -6, -1, -8,  7,  1, -7, -7, -5, -8, -8, -7, -2, -6, -7 }, // L
    .{ -5,-13,-11, -7, -4, -8,-10, -1, -2,  1, 11, -9, -8, -4, -4, -5, -4, -1,-13,-11 }, // M
    .{ -4,-11,  2, -2, -9, -3,  0, -5, -1, -7, -9,  8, -6, -3, -6,  0, -2, -8, -8, -4 }, // N
    .{ -2, -8, -8, -5,-10, -6, -4, -8, -6, -7, -8, -6,  8, -3, -4, -2, -4, -6,-14,-13 }, // P
    .{ -4,-14, -2,  1,-13, -7,  1, -8, -3, -5, -4, -3, -3,  8, -2, -5, -5, -7,-13,-12 }, // Q
    .{ -7, -8,-10, -9, -9, -9, -2, -5,  0, -8, -4, -6, -4, -2,  8, -3, -6, -8, -2,-10 }, // R
    .{  0, -3, -4, -4, -6, -2, -6, -7, -4, -8, -5,  0, -2, -5, -3,  6,  0, -6, -5, -7 }, // S
    .{ -1, -8, -5, -6, -9, -6, -7, -2, -3, -7, -4, -2, -4, -5, -6,  0,  7, -3,-13, -6 }, // T
    .{ -2, -6, -8, -6, -8, -5, -6,  2, -9, -2, -1, -8, -6, -7, -8, -6, -3,  7,-15, -7 }, // V
    .{-13,-15,-15,-17, -4,-15, -7,-14,-12, -6,-13, -8,-14,-13, -2, -5,-13,-15, 13, -5 }, // W
    .{ -8, -4,-11, -8,  2,-14, -3, -6, -9, -7,-11, -4,-13,-12,-10, -7, -6, -7, -5, 10 }, // Y
};

const pam30_data: [Kp][Kp]i16 = embedCanonical(pam30_canonical);

pub const pam30 = ScoreMatrix{
    .data = &pam30_data,
    .name = "PAM30",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// PAM70 data
// ---------------------------------------------------------------------------

// Data from Easel esl_scorematrix.c (authoritative PAM70 values).
const pam70_canonical: [20][20]i8 = .{
    //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    .{  5, -4, -1, -1, -6,  0, -4, -2, -4, -4, -3, -2,  0, -2, -4,  1,  1, -1, -9, -5 }, // A
    .{ -4,  9, -9, -9, -8, -6, -5, -4, -9,-10, -9, -7, -5, -9, -5, -1, -5, -4,-11, -2 }, // C
    .{ -1, -9,  6,  3,-10, -1, -1, -5, -2, -8, -7,  3, -4,  0, -6, -1, -2, -5,-10, -7 }, // D
    .{ -1, -9,  3,  6, -9, -2, -2, -4, -2, -6, -4,  0, -3,  2, -5, -2, -3, -4,-11, -6 }, // E
    .{ -6, -8,-10, -9,  8, -7, -4,  0, -9, -1, -2, -6, -7, -9, -7, -4, -6, -5, -2,  4 }, // F
    .{  0, -6, -1, -2, -7,  6, -6, -6, -5, -7, -6, -1, -3, -4, -6,  0, -3, -3,-10, -9 }, // G
    .{ -4, -5, -1, -2, -4, -6,  8, -6, -3, -4, -6,  1, -2,  2,  0, -3, -4, -4, -5, -1 }, // H
    .{ -2, -4, -5, -4,  0, -6, -6,  7, -4,  1,  1, -3, -5, -5, -3, -4, -1,  3, -9, -4 }, // I
    .{ -4, -9, -2, -2, -9, -5, -3, -4,  6, -5,  0,  0, -4, -1,  2, -2, -1, -6, -7, -7 }, // K
    .{ -4,-10, -8, -6, -1, -7, -4,  1, -5,  6,  2, -5, -5, -3, -6, -6, -4,  0, -4, -4 }, // L
    .{ -3, -9, -7, -4, -2, -6, -6,  1,  0,  2, 10, -5, -5, -2, -2, -3, -2,  0, -8, -7 }, // M
    .{ -2, -7,  3,  0, -6, -1,  1, -3,  0, -5, -5,  6, -3, -1, -3,  1,  0, -5, -6, -3 }, // N
    .{  0, -5, -4, -3, -7, -3, -2, -5, -4, -5, -5, -3,  7, -1, -2,  0, -2, -3, -9, -9 }, // P
    .{ -2, -9,  0,  2, -9, -4,  2, -5, -1, -3, -2, -1, -1,  7,  0, -3, -3, -4, -8, -8 }, // Q
    .{ -4, -5, -6, -5, -7, -6,  0, -3,  2, -6, -2, -3, -2,  0,  8, -1, -4, -5,  0, -7 }, // R
    .{  1, -1, -1, -2, -4,  0, -3, -4, -2, -6, -3,  1,  0, -3, -1,  5,  2, -3, -3, -5 }, // S
    .{  1, -5, -2, -3, -6, -3, -4, -1, -1, -4, -2,  0, -2, -3, -4,  2,  6, -1, -8, -4 }, // T
    .{ -1, -4, -5, -4, -5, -3, -4,  3, -6,  0,  0, -5, -3, -4, -5, -3, -1,  6,-10, -5 }, // V
    .{ -9,-11,-10,-11, -2,-10, -5, -9, -7, -4, -8, -6, -9, -8,  0, -3, -8,-10, 13, -3 }, // W
    .{ -5, -2, -7, -6,  4, -9, -1, -4, -7, -4, -7, -3, -9, -8, -7, -5, -4, -5, -3,  9 }, // Y
};

const pam70_data: [Kp][Kp]i16 = embedCanonical(pam70_canonical);

pub const pam70 = ScoreMatrix{
    .data = &pam70_data,
    .name = "PAM70",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// PAM120 data
// ---------------------------------------------------------------------------

const pam120_canonical: [20][20]i8 = .{
    .{  3, -3,  0,  0, -4,  1, -3, -1, -2, -3, -2, -1,  1, -1, -3,  1,  1,  0, -7, -4 }, // A
    .{ -3,  9, -7, -7, -6, -5, -4, -3, -7, -7, -6, -5, -4, -7, -4,  0, -3, -3, -8, -1 }, // C
    .{  0, -7,  5,  3, -7, -1,  0, -3, -1, -5, -4,  2, -3,  1, -3,  0, -1, -3, -8, -5 }, // D
    .{  0, -7,  3,  5, -7, -2, -1, -3, -1, -4, -3,  1, -2,  2, -3, -1, -2, -3, -8, -5 }, // E
    .{ -4, -6, -7, -7,  8, -5, -2,  0, -7,  0, -1, -4, -5, -6, -5, -3, -4, -3, -1,  4 }, // F
    .{  1, -5, -1, -2, -5,  5, -4, -4, -3, -5, -4,  0, -2, -3, -4,  1, -1, -2, -8, -6 }, // G
    .{ -3, -4,  0, -1, -2, -4,  7, -4, -2, -3, -4,  2, -1,  3,  1, -2, -3, -3, -3, -1 }, // H
    .{ -1, -3, -3, -3,  0, -4, -4,  6, -3,  1,  1, -2, -3, -3, -2, -2,  0,  3, -6, -2 }, // I
    .{ -2, -7, -1, -1, -7, -3, -2, -3,  5, -4,  0,  1, -2,  0,  2, -1, -1, -4, -5, -5 }, // K
    .{ -3, -7, -5, -4,  0, -5, -3,  1, -4,  5,  3, -4, -3, -2, -4, -4, -3,  1, -3, -2 }, // L
    .{ -2, -6, -4, -3, -1, -4, -4,  1,  0,  3,  8, -3, -3, -1, -1, -2, -1,  1, -6, -4 }, // M
    .{ -1, -5,  2,  1, -4,  0,  2, -2,  1, -4, -3,  4, -2,  0, -1,  1,  0, -3, -4, -2 }, // N
    .{  1, -4, -3, -2, -5, -2, -1, -3, -2, -3, -3, -2,  6, -1, -1,  1, -1, -2, -7, -6 }, // P
    .{ -1, -7,  1,  2, -6, -3,  3, -3,  0, -2, -1,  0, -1,  6,  1, -2, -2, -3, -6, -5 }, // Q
    .{ -3, -4, -3, -3, -5, -4,  1, -2,  2, -4, -1, -1, -1,  1,  6, -1, -2, -3,  1, -5 }, // R
    .{  1,  0,  0, -1, -3,  1, -2, -2, -1, -4, -2,  1,  1, -2, -1,  3,  2, -2, -2, -3 }, // S
    .{  1, -3, -1, -2, -4, -1, -3,  0, -1, -3, -1,  0, -1, -2, -2,  2,  4,  0, -6, -3 }, // T
    .{  0, -3, -3, -3, -3, -2, -3,  3, -4,  1,  1, -3, -2, -3, -3, -2,  0,  5, -8, -3 }, // V
    .{ -7, -8, -8, -8, -1, -8, -3, -6, -5, -3, -6, -4, -7, -6,  1, -2, -6, -8, 12, -2 }, // W
    .{ -4, -1, -5, -5,  4, -6, -1, -2, -5, -2, -4, -2, -6, -5, -5, -3, -3, -3, -2,  8 }, // Y
};

const pam120_data: [Kp][Kp]i16 = embedCanonical(pam120_canonical);

pub const pam120 = ScoreMatrix{
    .data = &pam120_data,
    .name = "PAM120",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// PAM250 data
// ---------------------------------------------------------------------------

const pam250_canonical: [20][20]i8 = .{
    .{  2, -2,  0,  0, -3,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3 }, // A
    .{ -2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0 }, // C
    .{  0, -5,  4,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4 }, // D
    .{  0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4 }, // E
    .{ -3, -4, -6, -5,  9, -5, -2,  1, -5,  2,  0, -3, -5, -5, -4, -3, -3, -1,  0,  7 }, // F
    .{  1, -3,  1,  0, -5,  5, -2, -3, -2, -4, -3,  0,  0, -1, -3,  1,  0, -1, -7, -5 }, // G
    .{ -1, -3,  1,  1, -2, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0 }, // H
    .{ -1, -2, -2, -2,  1, -3, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1 }, // I
    .{ -1, -5,  0,  0, -5, -2,  0, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4 }, // K
    .{ -2, -6, -4, -3,  2, -4, -2,  2, -3,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1 }, // L
    .{ -1, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6, -2, -2, -1,  0, -2, -1,  2, -4, -2 }, // M
    .{  0, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  2,  0,  1,  0,  1,  0, -2, -4, -2 }, // N
    .{  1, -3, -1, -1, -5,  0,  0, -2, -1, -3, -2,  0,  6,  0,  0,  1,  0, -1, -6, -5 }, // P
    .{  0, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4,  1, -1, -1, -2, -5, -4 }, // Q
    .{ -2, -4, -1, -1, -4, -3,  2, -2,  3, -3,  0,  0,  0,  1,  6,  0, -1, -2,  2, -4 }, // R
    .{  1,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1, -1,  0,  2,  1, -1, -2, -3 }, // S
    .{  1, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,  0, -5, -3 }, // T
    .{  0, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4, -6, -2 }, // V
    .{ -6, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17,  0 }, // W
    .{ -3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, 10 }, // Y
};

const pam250_data: [Kp][Kp]i16 = embedCanonical(pam250_canonical);

pub const pam250 = ScoreMatrix{
    .data = &pam250_data,
    .name = "PAM250",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// Identity matrix: +1 for match, -1 for mismatch.
// ---------------------------------------------------------------------------

const identity_data: [Kp][Kp]i16 = blk: {
    var m: [Kp][Kp]i16 = .{.{0} ** Kp} ** Kp;
    var i: usize = 0;
    while (i < 20) : (i += 1) {
        var j: usize = 0;
        while (j < 20) : (j += 1) {
            m[i][j] = if (i == j) 1 else -1;
        }
    }
    break :blk m;
};

pub const identity = ScoreMatrix{
    .data = &identity_data,
    .name = "IDENTITY",
    .abc = &alphabet_mod.amino,
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "blosum62: diagonal values" {
    try std.testing.expectEqual(@as(i16, 4), blosum62.score(0, 0));
    try std.testing.expectEqual(@as(i16, 9), blosum62.score(1, 1));
    try std.testing.expectEqual(@as(i16, 11), blosum62.score(18, 18));
}

test "blosum62: off-diagonal values" {
    try std.testing.expectEqual(@as(i16, 0), blosum62.score(0, 1));
}

test "blosum62: symmetry" {
    try std.testing.expect(blosum62.isSymmetric());
}

test "blosum62: scoreByChar valid" {
    try std.testing.expectEqual(@as(i16, 4), try blosum62.scoreByChar('A', 'A'));
    try std.testing.expectEqual(@as(i16, 11), try blosum62.scoreByChar('W', 'W'));
    try std.testing.expectEqual(@as(i16, 0), try blosum62.scoreByChar('A', 'C'));
}

test "blosum62: scoreByChar case-insensitive" {
    try std.testing.expectEqual(@as(i16, 4), try blosum62.scoreByChar('a', 'a'));
    try std.testing.expectEqual(@as(i16, 11), try blosum62.scoreByChar('w', 'W'));
}

test "blosum62: scoreByChar invalid character" {
    try std.testing.expectError(error.InvalidCharacter, blosum62.scoreByChar('!', 'A'));
    try std.testing.expectError(error.InvalidCharacter, blosum62.scoreByChar('A', '!'));
}

test "blosum62: scoreByChar degenerate residue" {
    // B is code 21 (D|N), should return 0 for comptime matrix (no degenerate scores set)
    try std.testing.expectEqual(@as(i16, 0), try blosum62.scoreByChar('B', 'A'));
}

test "identity: match and mismatch" {
    try std.testing.expectEqual(@as(i16, 1), identity.score(0, 0));
    try std.testing.expectEqual(@as(i16, -1), identity.score(0, 1));
    try std.testing.expectEqual(@as(i16, 1), identity.score(19, 19));
    try std.testing.expectEqual(@as(i16, -1), identity.score(0, 19));
}

test "identity: symmetry" {
    try std.testing.expect(identity.isSymmetric());
}

test "blosum62: known off-diagonal values" {
    try std.testing.expectEqual(@as(i16, 2), blosum62.score(2, 3));
    try std.testing.expectEqual(@as(i16, 3), blosum62.score(4, 19));
    try std.testing.expectEqual(@as(i16, -1), blosum62.score(0, 14));
    try std.testing.expectEqual(@as(i16, 8), blosum62.score(6, 6));
    try std.testing.expectEqual(@as(i16, 7), blosum62.score(12, 12));
}

test "maxScore and minScore" {
    try std.testing.expectEqual(@as(i16, 11), blosum62.maxScore());
    try std.testing.expectEqual(@as(i16, -4), blosum62.minScore());
}

test "expectedScore: BLOSUM62 with BL62 background is negative" {
    const e = blosum62.expectedScore(&composition.bl62, &composition.bl62);
    try std.testing.expect(e < 0);
}

test "getByName: known matrices" {
    try std.testing.expect(getByName("BLOSUM62") != null);
    try std.testing.expect(getByName("blosum62") != null);
    try std.testing.expect(getByName("BL62") != null);
    try std.testing.expect(getByName("BLOSUM45") != null);
    try std.testing.expect(getByName("BLOSUM80") != null);
    try std.testing.expect(getByName("PAM30") != null);
    try std.testing.expect(getByName("PAM70") != null);
    try std.testing.expect(getByName("PAM120") != null);
    try std.testing.expect(getByName("PAM250") != null);
    try std.testing.expect(getByName("IDENTITY") != null);
    try std.testing.expect(getByName("UNKNOWN") == null);
}

test "blosum45: symmetry" {
    try std.testing.expect(blosum45.isSymmetric());
}

test "setDegenerateScores: B averages D and N" {
    const allocator = std.testing.allocator;
    // Clone blosum62 into an owned mutable matrix
    const data_ptr = try allocator.create([Kp][Kp]i16);
    defer allocator.destroy(data_ptr);
    data_ptr.* = blosum62_data;

    var sm = ScoreMatrix{
        .data = data_ptr,
        .name = "test",
        .abc = &alphabet_mod.amino,
        .owned = data_ptr,
        .allocator = allocator,
    };

    sm.setDegenerateScores();

    // B = D|N (code 21). S[A][B] should be average of S[A][D] and S[A][N]
    // S[A][D] = -2, S[A][N] = -2 => avg = -2
    const b_code: u8 = 21;
    try std.testing.expectEqual(@as(i16, -2), sm.data[0][b_code]);

    // S[B][B] should be average of S[D|N][D|N] = (S[D][D]+S[D][N]+S[N][D]+S[N][N])/4
    // = (6 + 1 + 1 + 6) / 4 = 14/4 = 4 (rounded: (14 + 2) / 4 = 4)
    try std.testing.expectEqual(@as(i16, 4), sm.data[b_code][b_code]);
}

test "read: parse BLAST format matrix" {
    const allocator = std.testing.allocator;

    const matrix_text =
        \\# Test matrix
        \\   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
        \\A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
        \\C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
        \\D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
        \\E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
        \\F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
        \\G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
        \\H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
        \\I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
        \\K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
        \\L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
        \\M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
        \\N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
        \\P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
        \\Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
        \\R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
        \\S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
        \\T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
        \\V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
        \\W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
        \\Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
    ;

    var sm = try ScoreMatrix.read(allocator, &alphabet_mod.amino, matrix_text);
    defer sm.deinit();

    try std.testing.expectEqual(@as(i16, 4), sm.score(0, 0));
    try std.testing.expectEqual(@as(i16, 9), sm.score(1, 1));
    try std.testing.expectEqual(@as(i16, 11), sm.score(18, 18));
    try std.testing.expectEqual(@as(i16, 0), sm.score(0, 1));
    try std.testing.expect(sm.isSymmetric());
}

test "write: round-trip" {
    const allocator = std.testing.allocator;
    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);

    try blosum62.write(buf.writer(allocator).any());

    // Parse it back
    var sm = try ScoreMatrix.read(allocator, &alphabet_mod.amino, buf.items);
    defer sm.deinit();

    // Verify canonical scores match
    for (0..20) |i| {
        for (0..20) |j| {
            try std.testing.expectEqual(blosum62.data[i][j], sm.data[i][j]);
        }
    }
}

test "relativeEntropy: BLOSUM62 is positive" {
    const h = try blosum62.relativeEntropy(&composition.bl62, 0.3466);
    try std.testing.expect(h > 0);
}

test "relativeEntropy: invalid probabilities" {
    var bad_bg: [20]f64 = undefined;
    for (&bad_bg) |*v| v.* = 0.1; // sums to 2.0
    try std.testing.expectError(error.InvalidProbabilities, blosum62.relativeEntropy(&bad_bg, 0.3466));
}

test "relativeEntropy: invalid normalization with wrong lambda" {
    // A very wrong lambda should cause Z to deviate far from 1.0
    try std.testing.expectError(error.InvalidNormalization, blosum62.relativeEntropy(&composition.bl62, 5.0));
}

// ---------------------------------------------------------------------------
// Probabilistic score matrix derivation
// ---------------------------------------------------------------------------

/// Derive an integer score matrix from target probabilities and background
/// frequencies.
///
/// For each canonical residue pair (i, j):
///   s[i][j] = round( log(pij[i][j] / (fi[i] * fj[j])) / lambda )
///
/// Rounding uses "half away from zero" (standard rounding).
/// Only the K=20 canonical amino acid portion is filled; degenerate codes
/// are set to zero.
///
/// Reference: Easel `esl_scorematrix_SetFromProbs()`.
pub fn setFromProbs(
    allocator: Allocator,
    lambda: f64,
    pij: Matrix,
    fi: []const f64,
    fj: []const f64,
) !ScoreMatrix {
    const k: usize = 20;
    std.debug.assert(pij.rows >= k and pij.cols >= k);
    std.debug.assert(fi.len >= k and fj.len >= k);
    std.debug.assert(lambda > 0);

    const data_ptr = try allocator.create([Kp][Kp]i16);
    errdefer allocator.destroy(data_ptr);

    // Initialize to zero (covers degenerate, gap, missing positions).
    for (data_ptr) |*row| {
        @memset(row, 0);
    }

    for (0..k) |i| {
        for (0..k) |j| {
            const sc: f64 = @log(pij.get(i, j) / (fi[i] * fj[j])) / lambda;
            // Round to nearest integer, half away from zero.
            data_ptr[i][j] = @intFromFloat(@round(sc));
        }
    }

    return ScoreMatrix{
        .data = data_ptr,
        .name = "from_probs",
        .abc = &alphabet_mod.amino,
        .owned = data_ptr,
        .allocator = allocator,
    };
}

/// Result of probifyGivenBG: the solved scale parameter and target
/// probability matrix.
pub const ProbifyResult = struct {
    lambda: f64,
    pij: Matrix,

    pub fn deinit(self: *ProbifyResult) void {
        self.pij.deinit();
    }
};

/// Given a score matrix and background frequencies, solve for the scale
/// parameter lambda and recover the implicit target probabilities.
///
/// Uses Newton-Raphson to find lambda such that:
///   sum_ij fi[i] * fj[j] * exp(lambda * s[i][j]) = 1.0
///
/// Then computes:
///   pij[i][j] = fi[i] * fj[j] * exp(lambda * s[i][j])
///
/// Reference: Easel `esl_scorematrix_ProbifyGivenBG()`.
pub fn probifyGivenBG(
    allocator: Allocator,
    sm: ScoreMatrix,
    fi: []const f64,
    fj: []const f64,
) !ProbifyResult {
    const k: usize = 20;
    std.debug.assert(fi.len >= k and fj.len >= k);

    // Find the maximum score to set the initial lambda guess.
    const s_max: f64 = @floatFromInt(sm.maxScore());
    if (s_max <= 0) return error.InvalidScoreMatrix;

    // Bracket the root: start from 1/s_max and double until f(lambda) > 0.
    // f(lambda) = sum_ij fi[i]*fj[j]*exp(lambda*s[i][j]) - 1.0
    var lambda_guess: f64 = 1.0 / s_max;
    var fx: f64 = -1.0;
    while (lambda_guess < 50.0) {
        fx = evalLambdaFunc(sm, fi, fj, k, lambda_guess);
        if (fx > 0) break;
        lambda_guess *= 2.0;
    }
    if (fx <= 0) return error.FailedToBracket;

    // Newton-Raphson iteration.
    const max_iter: usize = 100;
    const abs_tol: f64 = 1e-15;
    const rel_tol: f64 = 1e-15;

    var lambda: f64 = lambda_guess;
    fx = evalLambdaFunc(sm, fi, fj, k, lambda);
    var dfx: f64 = evalLambdaDeriv(sm, fi, fj, k, lambda);

    for (0..max_iter) |_| {
        const x0 = lambda;
        lambda = lambda - fx / dfx;

        fx = evalLambdaFunc(sm, fi, fj, k, lambda);
        dfx = evalLambdaDeriv(sm, fi, fj, k, lambda);

        if (fx == 0.0) break;
        if (@abs(lambda - x0) < abs_tol + rel_tol * lambda or @abs(fx) < abs_tol) break;
    } else {
        return error.DidNotConverge;
    }

    // Compute target probabilities: pij[i][j] = fi[i] * fj[j] * exp(lambda * s[i][j])
    var pij = try Matrix.init(allocator, k, k);
    errdefer pij.deinit();

    for (0..k) |i| {
        for (0..k) |j| {
            const s: f64 = @floatFromInt(sm.data[i][j]);
            pij.set(i, j, fi[i] * fj[j] * @exp(lambda * s));
        }
    }

    return ProbifyResult{ .lambda = lambda, .pij = pij };
}

/// Evaluate f(lambda) = sum_ij fi[i]*fj[j]*exp(lambda*s[i][j]) - 1.0
fn evalLambdaFunc(
    sm: ScoreMatrix,
    fi: []const f64,
    fj: []const f64,
    k: usize,
    lambda: f64,
) f64 {
    var total: f64 = 0;
    for (0..k) |i| {
        for (0..k) |j| {
            const s: f64 = @floatFromInt(sm.data[i][j]);
            total += fi[i] * fj[j] * @exp(lambda * s);
        }
    }
    return total - 1.0;
}

/// Evaluate f'(lambda) = sum_ij fi[i]*fj[j]*s[i][j]*exp(lambda*s[i][j])
fn evalLambdaDeriv(
    sm: ScoreMatrix,
    fi: []const f64,
    fj: []const f64,
    k: usize,
    lambda: f64,
) f64 {
    var total: f64 = 0;
    for (0..k) |i| {
        for (0..k) |j| {
            const s: f64 = @floatFromInt(sm.data[i][j]);
            total += fi[i] * fj[j] * s * @exp(lambda * s);
        }
    }
    return total;
}

test "setFromProbs: round-trip with BLOSUM62" {
    // Use probifyGivenBG to get lambda and pij from BLOSUM62,
    // then reconstruct the score matrix and verify it matches.
    const allocator = std.testing.allocator;

    var result = try probifyGivenBG(allocator, blosum62, &composition.bl62, &composition.bl62);
    defer result.deinit();

    var sm = try setFromProbs(allocator, result.lambda, result.pij, &composition.bl62, &composition.bl62);
    defer sm.deinit();

    // The reconstructed scores should match the original for all canonical pairs.
    for (0..20) |i| {
        for (0..20) |j| {
            try std.testing.expectEqual(blosum62.data[i][j], sm.data[i][j]);
        }
    }
}

test "probifyGivenBG: lambda is positive and probabilities sum to ~1.0" {
    const allocator = std.testing.allocator;

    var result = try probifyGivenBG(allocator, blosum62, &composition.bl62, &composition.bl62);
    defer result.deinit();

    // Lambda must be positive for a valid score matrix.
    try std.testing.expect(result.lambda > 0);

    // Target probabilities must sum to approximately 1.0.
    var total: f64 = 0;
    for (0..20) |i| {
        for (0..20) |j| {
            const p = result.pij.get(i, j);
            try std.testing.expect(p > 0);
            total += p;
        }
    }
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), total, 1e-10);
}

test "probifyGivenBG: known BLOSUM62 lambda" {
    const allocator = std.testing.allocator;

    var result = try probifyGivenBG(allocator, blosum62, &composition.bl62, &composition.bl62);
    defer result.deinit();

    // Lambda for BLOSUM62 should be in a reasonable range (0.3-0.4).
    // The exact value depends on the background frequency vector used.
    try std.testing.expect(result.lambda > 0.3 and result.lambda < 0.4);
}

test "setFromProbs: scores are symmetric when inputs are symmetric" {
    const allocator = std.testing.allocator;

    // Build a symmetric pij from BLOSUM62 background.
    var pij = try Matrix.init(allocator, 20, 20);
    defer pij.deinit();

    const lambda: f64 = 0.3;
    // Create symmetric target probs from symmetric scores.
    var z: f64 = 0;
    for (0..20) |i| {
        for (0..20) |j| {
            const s: f64 = @floatFromInt(blosum62.data[i][j]);
            const v = composition.bl62[i] * composition.bl62[j] * @exp(lambda * s);
            pij.set(i, j, v);
            z += v;
        }
    }
    // Normalize.
    for (0..20) |i| {
        for (0..20) |j| {
            pij.set(i, j, pij.get(i, j) / z);
        }
    }

    var sm = try setFromProbs(allocator, lambda, pij, &composition.bl62, &composition.bl62);
    defer sm.deinit();

    // BLOSUM62 is symmetric, so the result should be too.
    for (0..20) |i| {
        for (i + 1..20) |j| {
            try std.testing.expectEqual(sm.data[i][j], sm.data[j][i]);
        }
    }
}

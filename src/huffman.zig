// Huffman coding for digitized alphabets.
//
// Builds a canonical Huffman code from symbol frequencies and provides
// encode/decode operations. Used for compressed sequence storage.

const std = @import("std");
const Allocator = std.mem.Allocator;

/// Maximum code length supported (u5 max = 31).
const MAX_CODE_LEN = 31;

pub const HuffmanCode = struct {
    /// Code word for each symbol (bit pattern stored in u32).
    codes: []u32,
    /// Code length (bits) for each symbol.
    lengths: []u5,
    /// Number of symbols.
    n_symbols: usize,
    /// Decoding table indexed by code length (1..MAX_CODE_LEN).
    /// For each length, stores (first_code, symbol_array) for O(1) lookup.
    decode_first_code: [MAX_CODE_LEN + 1]u32,
    decode_symbols: [MAX_CODE_LEN + 1]?[]u8,
    allocator: Allocator,

    /// Build a Huffman code from symbol frequencies.
    /// freq[i] = frequency of symbol i, for i in 0..n_symbols-1.
    pub fn build(allocator: Allocator, freq: []const f64) !HuffmanCode {
        const n = freq.len;
        if (n == 0) return error.InvalidInput;

        // Count symbols with non-zero frequency
        var n_active: usize = 0;
        for (freq) |f| {
            if (f > 0) n_active += 1;
        }
        if (n_active == 0) return error.InvalidInput;

        // Allocate output arrays; zero-frequency symbols get code 0, length 0
        const lengths = try allocator.alloc(u5, n);
        errdefer allocator.free(lengths);
        @memset(lengths, 0);
        const codes = try allocator.alloc(u32, n);
        errdefer allocator.free(codes);
        @memset(codes, 0);

        if (n_active == 1) {
            var sole_symbol: usize = 0;
            for (0..n) |i| {
                if (freq[i] > 0) {
                    lengths[i] = 1;
                    codes[i] = 0;
                    sole_symbol = i;
                    break;
                }
            }
            // Build decode table for the single-symbol case.
            var dt_first: [MAX_CODE_LEN + 1]u32 = @splat(0);
            var dt_syms: [MAX_CODE_LEN + 1]?[]u8 = @splat(null);
            const sym_arr = try allocator.alloc(u8, 1);
            sym_arr[0] = @intCast(sole_symbol);
            dt_first[1] = 0;
            dt_syms[1] = sym_arr;
            return .{
                .codes = codes,
                .lengths = lengths,
                .n_symbols = n,
                .decode_first_code = dt_first,
                .decode_symbols = dt_syms,
                .allocator = allocator,
            };
        }

        // Build indices of active symbols only, sorted by frequency (ascending)
        const indices = try allocator.alloc(usize, n_active);
        defer allocator.free(indices);
        {
            var idx: usize = 0;
            for (0..n) |i| {
                if (freq[i] > 0) {
                    indices[idx] = i;
                    idx += 1;
                }
            }
        }

        std.mem.sort(usize, indices, freq, struct {
            fn lessThan(f: []const f64, a: usize, b: usize) bool {
                return f[a] < f[b];
            }
        }.lessThan);

        // Build tree bottom-up to determine lengths
        const tree_size = 2 * n_active - 1;
        const tree_freq = try allocator.alloc(f64, tree_size);
        defer allocator.free(tree_freq);
        const tree_left = try allocator.alloc(?usize, tree_size);
        defer allocator.free(tree_left);
        const tree_right = try allocator.alloc(?usize, tree_size);
        defer allocator.free(tree_right);

        // Initialize leaves
        for (0..n_active) |i| {
            tree_freq[i] = freq[indices[i]];
            tree_left[i] = null;
            tree_right[i] = null;
        }

        // Priority queue using sorted merge
        var active = try allocator.alloc(bool, tree_size);
        defer allocator.free(active);
        @memset(active, false);
        for (0..n_active) |i| active[i] = true;

        var next_internal: usize = n_active;
        for (0..n_active - 1) |_| {
            // Find two smallest active nodes
            var m1: usize = tree_size;
            var m2: usize = tree_size;
            var m1_freq: f64 = std.math.inf(f64);
            var m2_freq: f64 = std.math.inf(f64);

            for (0..next_internal) |i| {
                if (!active[i]) continue;
                if (tree_freq[i] < m1_freq) {
                    m2 = m1;
                    m2_freq = m1_freq;
                    m1 = i;
                    m1_freq = tree_freq[i];
                } else if (tree_freq[i] < m2_freq) {
                    m2 = i;
                    m2_freq = tree_freq[i];
                }
            }

            // Create internal node
            tree_freq[next_internal] = m1_freq + m2_freq;
            tree_left[next_internal] = m1;
            tree_right[next_internal] = m2;
            active[next_internal] = true;
            active[m1] = false;
            active[m2] = false;
            next_internal += 1;
        }

        // Traverse tree to assign lengths
        const depth = try allocator.alloc(u5, tree_size);
        defer allocator.free(depth);
        @memset(depth, 0);

        // BFS from root (last created node)
        var ii: usize = next_internal;
        while (ii > 0) {
            ii -= 1;
            if (tree_left[ii]) |l| depth[l] = depth[ii] +| 1;
            if (tree_right[ii]) |r| depth[r] = depth[ii] +| 1;
        }

        // Map leaf depths back to original symbols (active symbols only)
        for (0..n_active) |i| {
            lengths[indices[i]] = if (depth[i] > 0) depth[i] else 1;
        }

        // Generate canonical codes from lengths (active symbols only)
        const code_order = try allocator.alloc(usize, n_active);
        defer allocator.free(code_order);
        for (0..n_active) |i| code_order[i] = indices[i];
        std.mem.sort(usize, code_order, lengths, struct {
            fn lessThan(l: []const u5, a: usize, b: usize) bool {
                if (l[a] != l[b]) return l[a] < l[b];
                return a < b;
            }
        }.lessThan);

        var code: u32 = 0;
        var prev_len: u5 = lengths[code_order[0]];
        codes[code_order[0]] = 0;
        for (1..n_active) |i| {
            code += 1;
            const curr_len = lengths[code_order[i]];
            if (curr_len > prev_len) {
                code <<= @intCast(curr_len - prev_len);
            }
            codes[code_order[i]] = code;
            prev_len = curr_len;
        }

        // Build decoding table: for each code length, store first_code and
        // a symbol array so that decode can do O(1) lookup per length.
        var dt_first: [MAX_CODE_LEN + 1]u32 = @splat(0);
        var dt_syms: [MAX_CODE_LEN + 1]?[]u8 = @splat(null);

        // code_order is sorted by (length, symbol). Walk it to find runs of same length.
        var run_start: usize = 0;
        while (run_start < n_active) {
            const run_len_val = lengths[code_order[run_start]];
            var run_end: usize = run_start + 1;
            while (run_end < n_active and lengths[code_order[run_end]] == run_len_val) {
                run_end += 1;
            }
            const run_count = run_end - run_start;
            const sym_arr = try allocator.alloc(u8, run_count);
            for (0..run_count) |j| {
                sym_arr[j] = @intCast(code_order[run_start + j]);
            }
            dt_first[run_len_val] = codes[code_order[run_start]];
            dt_syms[run_len_val] = sym_arr;
            run_start = run_end;
        }

        return .{
            .codes = codes,
            .lengths = lengths,
            .n_symbols = n,
            .decode_first_code = dt_first,
            .decode_symbols = dt_syms,
            .allocator = allocator,
        };
    }

    /// Average code length (bits per symbol) given frequencies.
    pub fn avgLength(self: HuffmanCode, freq: []const f64) f64 {
        var total_freq: f64 = 0;
        var weighted: f64 = 0;
        for (0..self.n_symbols) |i| {
            total_freq += freq[i];
            weighted += freq[i] * @as(f64, @floatFromInt(self.lengths[i]));
        }
        if (total_freq == 0) return 0;
        return weighted / total_freq;
    }

    /// Encode a sequence of symbols into a packed bitstream.
    /// Returns the packed bytes and the total number of valid bits.
    pub fn encode(self: HuffmanCode, allocator: Allocator, symbols: []const u8) !struct { data: []u8, nbits: usize } {
        // Compute total bits needed
        var total_bits: usize = 0;
        for (symbols) |s| {
            if (s >= self.n_symbols) return error.InvalidSymbol;
            if (self.lengths[s] == 0) return error.InvalidSymbol;
            total_bits += @as(usize, self.lengths[s]);
        }
        const nbytes = (total_bits + 7) / 8;
        const data = try allocator.alloc(u8, nbytes);
        @memset(data, 0);

        var bit_pos: usize = 0;
        for (symbols) |s| {
            const code = self.codes[s];
            const len = @as(usize, self.lengths[s]);
            // Write bits MSB-first
            for (0..len) |b| {
                const bit: u1 = @intCast((code >> @intCast(len - 1 - b)) & 1);
                if (bit == 1) {
                    data[bit_pos / 8] |= @as(u8, 1) << @intCast(7 - (bit_pos % 8));
                }
                bit_pos += 1;
            }
        }

        return .{ .data = data, .nbits = total_bits };
    }

    /// Decode a packed bitstream back into symbols.
    /// `nbits` is the number of valid bits in the stream.
    /// Returns the decoded symbol sequence.
    /// Uses the precomputed decoding table for O(Lmax) per symbol instead of O(n*Lmax).
    pub fn decode(self: HuffmanCode, allocator: Allocator, data: []const u8, nbits: usize) ![]u8 {
        var result: std.ArrayList(u8) = .empty;
        errdefer result.deinit(allocator);

        var bit_pos: usize = 0;
        while (bit_pos < nbits) {
            var code: u32 = 0;
            var found = false;
            for (1..MAX_CODE_LEN + 1) |len| {
                if (bit_pos >= nbits) break;
                const bit: u32 = @intCast((data[bit_pos / 8] >> @intCast(7 - (bit_pos % 8))) & 1);
                code = (code << 1) | bit;
                bit_pos += 1;
                // Look up in the decoding table for this code length.
                if (self.decode_symbols[len]) |symbols| {
                    const first = self.decode_first_code[len];
                    if (code >= first and code - first < symbols.len) {
                        try result.append(allocator, symbols[code - first]);
                        found = true;
                        break;
                    }
                }
            }
            if (!found) return error.InvalidBitstream;
        }

        return result.toOwnedSlice(allocator);
    }

    pub fn deinit(self: *HuffmanCode) void {
        for (&self.decode_symbols) |*slot| {
            if (slot.*) |s| {
                self.allocator.free(s);
                slot.* = null;
            }
        }
        self.allocator.free(self.codes);
        self.allocator.free(self.lengths);
    }
};

// --- Tests ---

test "build: uniform frequencies give equal-ish lengths" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 1.0, 1.0, 1.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    try std.testing.expectEqual(@as(usize, 4), hc.n_symbols);
    // All lengths should be 2 for 4 uniform symbols
    for (hc.lengths) |l| {
        try std.testing.expectEqual(@as(u5, 2), l);
    }
}

test "build: skewed frequencies give different lengths" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 100.0, 1.0, 1.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    // Most frequent symbol should have shortest code
    try std.testing.expect(hc.lengths[0] <= hc.lengths[1]);
    try std.testing.expect(hc.lengths[0] <= hc.lengths[2]);
    try std.testing.expect(hc.lengths[0] <= hc.lengths[3]);
}

test "build: 2 symbols" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 3.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    try std.testing.expectEqual(@as(u5, 1), hc.lengths[0]);
    try std.testing.expectEqual(@as(u5, 1), hc.lengths[1]);
    // Codes should be 0 and 1
    try std.testing.expect(hc.codes[0] != hc.codes[1]);
}

test "encode/decode: round trip" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 10.0, 5.0, 2.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    const symbols = [_]u8{ 0, 1, 2, 3, 0, 0, 1, 3, 2 };
    const encoded = try hc.encode(allocator, &symbols);
    defer allocator.free(encoded.data);

    const decoded = try hc.decode(allocator, encoded.data, encoded.nbits);
    defer allocator.free(decoded);

    try std.testing.expectEqualSlices(u8, &symbols, decoded);
}

test "encode/decode: single symbol repeated" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 100.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    const symbols = [_]u8{ 0, 0, 0, 0 };
    const encoded = try hc.encode(allocator, &symbols);
    defer allocator.free(encoded.data);

    const decoded = try hc.decode(allocator, encoded.data, encoded.nbits);
    defer allocator.free(decoded);

    try std.testing.expectEqualSlices(u8, &symbols, decoded);
}

test "build: zero-frequency symbols excluded from tree" {
    const allocator = std.testing.allocator;
    // Symbols 1 and 3 have zero frequency
    const freq = [_]f64{ 10.0, 0.0, 5.0, 0.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    try std.testing.expectEqual(@as(usize, 5), hc.n_symbols);
    // Zero-frequency symbols get length 0 and code 0
    try std.testing.expectEqual(@as(u5, 0), hc.lengths[1]);
    try std.testing.expectEqual(@as(u32, 0), hc.codes[1]);
    try std.testing.expectEqual(@as(u5, 0), hc.lengths[3]);
    try std.testing.expectEqual(@as(u32, 0), hc.codes[3]);
    // Non-zero symbols get positive lengths
    try std.testing.expect(hc.lengths[0] > 0);
    try std.testing.expect(hc.lengths[2] > 0);
    try std.testing.expect(hc.lengths[4] > 0);

    // Round-trip with active symbols only
    const symbols = [_]u8{ 0, 2, 4, 0, 2 };
    const encoded = try hc.encode(allocator, &symbols);
    defer allocator.free(encoded.data);
    const decoded = try hc.decode(allocator, encoded.data, encoded.nbits);
    defer allocator.free(decoded);
    try std.testing.expectEqualSlices(u8, &symbols, decoded);

    // Encoding a zero-frequency symbol should fail
    const bad = [_]u8{1};
    try std.testing.expectError(error.InvalidSymbol, hc.encode(allocator, &bad));
}

test "decode table: direct lookup matches encode for 8-symbol alphabet" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 50.0, 30.0, 10.0, 5.0, 3.0, 1.0, 0.5, 0.5 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    // Verify decode table entries exist for each active code length.
    for (0..hc.n_symbols) |s| {
        const len = hc.lengths[s];
        if (len == 0) continue;
        const symbols = hc.decode_symbols[len] orelse {
            return error.InvalidBitstream;
        };
        const first = hc.decode_first_code[len];
        const offset = hc.codes[s] - first;
        try std.testing.expectEqual(@as(u8, @intCast(s)), symbols[offset]);
    }

    // Round-trip all symbols.
    const symbols = [_]u8{ 0, 1, 2, 3, 4, 5, 6, 7, 7, 6, 5, 4, 3, 2, 1, 0 };
    const encoded = try hc.encode(allocator, &symbols);
    defer allocator.free(encoded.data);
    const decoded = try hc.decode(allocator, encoded.data, encoded.nbits);
    defer allocator.free(decoded);
    try std.testing.expectEqualSlices(u8, &symbols, decoded);
}

test "avgLength: uniform 4 symbols = 2 bits" {
    const allocator = std.testing.allocator;
    const freq = [_]f64{ 1.0, 1.0, 1.0, 1.0 };
    var hc = try HuffmanCode.build(allocator, &freq);
    defer hc.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, 2.0), hc.avgLength(&freq), 1e-9);
}

// Genetic code and DNA-to-protein translation.
//
// Codons are indexed as base1*16 + base2*4 + base3, using the zeacel DNA
// digital encoding: A=0, C=1, G=2, T=3.
//
// Amino acid digital codes follow the zeacel amino alphabet:
// A=0 C=1 D=2 E=3 F=4 G=5 H=6 I=7 K=8 L=9 M=10 N=11 P=12 Q=13 R=14 S=15
// T=16 V=17 W=18 Y=19
// Stop codon is represented as 255.

const std = @import("std");
const Allocator = std.mem.Allocator;

/// Sentinel value meaning "stop codon" in a codon table.
pub const STOP_CODON: u8 = 255;

/// A codon table mapping 64 codons to amino acid digital codes.
/// Codons are encoded as (base1 * 16 + base2 * 4 + base3) using DNA digital
/// codes (A=0, C=1, G=2, T=3).
pub const GeneticCode = struct {
    /// Amino acid digital code for each of the 64 codons.
    /// STOP_CODON (255) marks stop codons.
    codon_table: [64]u8,
    /// Which codons are start codons.
    is_start: [64]bool,
    name: []const u8,

    /// Translate a single codon given three DNA digital codes.
    /// Returns the amino acid digital code, or null for stop codons or
    /// degenerate/gap bases (values >= 4).
    pub fn translate(self: GeneticCode, c1: u8, c2: u8, c3: u8) ?u8 {
        if (c1 >= 4 or c2 >= 4 or c3 >= 4) return null;
        const idx = @as(usize, c1) * 16 + @as(usize, c2) * 4 + @as(usize, c3);
        const aa = self.codon_table[idx];
        if (aa == STOP_CODON) return null;
        return aa;
    }

    /// Return true if the codon formed by the three DNA digital codes is a
    /// start codon.  Degenerate/gap bases always return false.
    pub fn isStart(self: GeneticCode, c1: u8, c2: u8, c3: u8) bool {
        if (c1 >= 4 or c2 >= 4 or c3 >= 4) return false;
        const idx = @as(usize, c1) * 16 + @as(usize, c2) * 4 + @as(usize, c3);
        return self.is_start[idx];
    }

    /// Translate a DNA digital sequence to a protein digital sequence.
    /// Reads in frame starting at offset 0.  Translation stops at the first
    /// stop codon or when fewer than 3 bases remain.
    ///
    /// The caller owns the returned slice and must free it with the same
    /// allocator.
    pub fn translateSeq(
        self: GeneticCode,
        allocator: Allocator,
        dna_dsq: []const u8,
    ) ![]u8 {
        const n_codons = dna_dsq.len / 3;
        var result: std.ArrayList(u8) = .empty;
        errdefer result.deinit(allocator);
        for (0..n_codons) |i| {
            const aa = self.translate(
                dna_dsq[i * 3],
                dna_dsq[i * 3 + 1],
                dna_dsq[i * 3 + 2],
            );
            if (aa) |a| {
                try result.append(allocator, a);
            } else {
                break; // stop codon or degenerate base
            }
        }
        return result.toOwnedSlice(allocator);
    }

    /// Translate in all 3 forward frames.
    /// Returns an array of 3 protein sequences; the caller owns each slice and
    /// must free them individually with the same allocator.
    pub fn translateThreeFrames(
        self: GeneticCode,
        allocator: Allocator,
        dna_dsq: []const u8,
    ) ![3][]u8 {
        var frames: [3][]u8 = undefined;
        var built: usize = 0;
        errdefer for (0..built) |k| allocator.free(frames[k]);
        for (0..3) |frame| {
            if (frame < dna_dsq.len) {
                frames[frame] = try self.translateSeq(allocator, dna_dsq[frame..]);
            } else {
                frames[frame] = try allocator.alloc(u8, 0);
            }
            built += 1;
        }
        return frames;
    }
};

// ---------------------------------------------------------------------------
// Comptime construction helpers
// ---------------------------------------------------------------------------

/// Return the index into a 64-entry codon table for three DNA digital codes.
fn codonIndex(b1: u8, b2: u8, b3: u8) usize {
    return @as(usize, b1) * 16 + @as(usize, b2) * 4 + @as(usize, b3);
}

/// Build a comptime ASCII-to-amino-digital-code map.
/// Non-amino characters map to 255.
fn buildAaMap() [128]u8 {
    // zeacel amino alphabet canonical symbols in order:
    const amino_symbols = "ACDEFGHIKLMNPQRSTVWY";
    var map: [128]u8 = .{255} ** 128;
    for (amino_symbols, 0..) |c, i| {
        map[c] = @intCast(i);
    }
    return map;
}

/// Map a DNA base ASCII character (A/C/G/T) to its zeacel digital code.
fn baseCharToDigital(c: u8) u8 {
    return switch (c) {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        else => unreachable,
    };
}

// ---------------------------------------------------------------------------
// Standard genetic code (NCBI translation table 1)
// ---------------------------------------------------------------------------
//
// The NCBI codon string is:
//   "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
//
// Its position i encodes the codon:
//   base1 = "TCAG"[i / 16]
//   base2 = "TCAG"[(i / 4) % 4]
//   base3 = "TCAG"[i % 4]
//
// We convert each base to our digital code (T=3, C=1, A=0, G=2) and place
// the amino acid digital code at the corresponding index.

/// Standard genetic code (NCBI translation table 1).
pub const standard: GeneticCode = blk: {
    const ncbi_order = "TCAG";
    const ncbi_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    const aa_map = buildAaMap();

    var table: [64]u8 = undefined;
    var starts: [64]bool = .{false} ** 64;

    // Populate codon table from NCBI string.
    var i: usize = 0;
    while (i < 64) : (i += 1) {
        const b1 = baseCharToDigital(ncbi_order[i / 16]);
        const b2 = baseCharToDigital(ncbi_order[(i / 4) % 4]);
        const b3 = baseCharToDigital(ncbi_order[i % 4]);
        const idx = codonIndex(b1, b2, b3);
        const aa_char = ncbi_aa[i];
        if (aa_char == '*') {
            table[idx] = STOP_CODON;
        } else {
            table[idx] = aa_map[aa_char];
        }
    }

    // ATG is the canonical start codon.
    // A=0, T=3, G=2 -> index = 0*16 + 3*4 + 2 = 14
    starts[codonIndex(0, 3, 2)] = true; // ATG

    break :blk GeneticCode{
        .codon_table = table,
        .is_start = starts,
        .name = "Standard",
    };
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "translate ATG -> M (10)" {
    // ATG: A=0, T=3, G=2
    const aa = standard.translate(0, 3, 2);
    try std.testing.expectEqual(@as(?u8, 10), aa); // M=10
}

test "translate TAA -> null (stop)" {
    // TAA: T=3, A=0, A=0
    const aa = standard.translate(3, 0, 0);
    try std.testing.expectEqual(@as(?u8, null), aa);
}

test "translate TAG -> null (stop)" {
    // TAG: T=3, A=0, G=2
    const aa = standard.translate(3, 0, 2);
    try std.testing.expectEqual(@as(?u8, null), aa);
}

test "translate TGA -> null (stop)" {
    // TGA: T=3, G=2, A=0
    const aa = standard.translate(3, 2, 0);
    try std.testing.expectEqual(@as(?u8, null), aa);
}

test "translate TTT -> F (4)" {
    // TTT: T=3, T=3, T=3
    const aa = standard.translate(3, 3, 3);
    try std.testing.expectEqual(@as(?u8, 4), aa); // F=4
}

test "translate GCT -> A (0)" {
    // GCT: G=2, C=1, T=3
    const aa = standard.translate(2, 1, 3);
    try std.testing.expectEqual(@as(?u8, 0), aa); // A=0
}

test "translate degenerate base -> null" {
    // N=15 is a degenerate base (>= 4)
    const aa = standard.translate(15, 0, 0);
    try std.testing.expectEqual(@as(?u8, null), aa);
}

test "isStart ATG -> true" {
    try std.testing.expect(standard.isStart(0, 3, 2)); // ATG
}

test "isStart TTT -> false" {
    try std.testing.expect(!standard.isStart(3, 3, 3)); // TTT
}

test "isStart TAA -> false (stop, not start)" {
    try std.testing.expect(!standard.isStart(3, 0, 0)); // TAA
}

test "translateSeq ATGGCTTAT -> [M, A, Y]" {
    // ATG=M(10), GCT=A(0), TAT=Y(19)
    const dna_text = "ATGGCTTAT";
    const allocator = std.testing.allocator;
    var dna_dsq: [9]u8 = undefined;
    const dna_alpha_encode = [128]u8{
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 0,   255, 1,   255, 255, 255, 2,   255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 3,   255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    };
    for (dna_text, 0..) |c, i| {
        dna_dsq[i] = dna_alpha_encode[c];
    }

    const protein = try standard.translateSeq(allocator, &dna_dsq);
    defer allocator.free(protein);

    try std.testing.expectEqualSlices(u8, &[_]u8{ 10, 0, 19 }, protein);
}

test "translateSeq stops at stop codon ATGTAAATG -> [M]" {
    // ATG=M(10), TAA=stop, ATG=M(10)  -- should stop after first M
    // Digital: A=0,T=3,G=2, T=3,A=0,A=0, A=0,T=3,G=2
    const dna_dsq = [_]u8{ 0, 3, 2, 3, 0, 0, 0, 3, 2 };
    const allocator = std.testing.allocator;
    const protein = try standard.translateSeq(allocator, &dna_dsq);
    defer allocator.free(protein);

    try std.testing.expectEqualSlices(u8, &[_]u8{10}, protein);
}

test "translateSeq empty input -> empty output" {
    const allocator = std.testing.allocator;
    const protein = try standard.translateSeq(allocator, &[_]u8{});
    defer allocator.free(protein);
    try std.testing.expectEqual(@as(usize, 0), protein.len);
}

test "translateSeq shorter than one codon -> empty output" {
    const allocator = std.testing.allocator;
    const protein = try standard.translateSeq(allocator, &[_]u8{ 0, 3 });
    defer allocator.free(protein);
    try std.testing.expectEqual(@as(usize, 0), protein.len);
}

test "translateThreeFrames produces three distinct translations" {
    // ATGGCTTATGCTGCGTTA
    // Frame 0: ATG GCT TAT GCT GCG TTA -> M A Y A A L
    // Frame 1: TGG CTT ATG CTG CGT TA  -> W L M L R (partial)
    // Frame 2: GGC TTA TGC TGC GTT A   -> G L C C V (partial)
    //
    // Digital encode: A=0, C=1, G=2, T=3
    // A T G G C T T A T G C T G C G T T A
    // 0 3 2 2 1 3 3 0 3 2 1 3 2 1 2 3 3 0
    const dna_dsq = [_]u8{ 0, 3, 2, 2, 1, 3, 3, 0, 3, 2, 1, 3, 2, 1, 2, 3, 3, 0 };
    const allocator = std.testing.allocator;

    const frames = try standard.translateThreeFrames(allocator, &dna_dsq);
    defer for (frames) |f| allocator.free(f);

    // Frame 0 should start with M (10)
    try std.testing.expect(frames[0].len > 0);
    try std.testing.expectEqual(@as(u8, 10), frames[0][0]); // M

    // All three frames should be slices (may differ in length and content)
    // Frame 1 starts at offset 1: TGG... -> W=18
    try std.testing.expect(frames[1].len > 0);
    try std.testing.expectEqual(@as(u8, 18), frames[1][0]); // W

    // Frame 2 starts at offset 2: GGC... -> G=5
    try std.testing.expect(frames[2].len > 0);
    try std.testing.expectEqual(@as(u8, 5), frames[2][0]); // G
}

test "all 64 codon table entries are set (no undefined)" {
    for (standard.codon_table) |aa| {
        // Each entry is either a valid amino digital code (0..19) or STOP_CODON.
        const valid = (aa < 20) or (aa == STOP_CODON);
        try std.testing.expect(valid);
    }
}

test "exactly 3 stop codons in standard code" {
    var stop_count: usize = 0;
    for (standard.codon_table) |aa| {
        if (aa == STOP_CODON) stop_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 3), stop_count);
}

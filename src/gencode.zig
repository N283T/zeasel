// Genetic code and DNA-to-protein translation.
//
// Codons are indexed as base1*16 + base2*4 + base3, using the zeasel DNA
// digital encoding: A=0, C=1, G=2, T=3.
//
// Amino acid digital codes follow the zeasel amino alphabet:
// A=0 C=1 D=2 E=3 F=4 G=5 H=6 I=7 K=8 L=9 M=10 N=11 P=12 Q=13 R=14 S=15
// T=16 V=17 W=18 Y=19
// Stop codon is represented as 255.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("alphabet.zig").Alphabet;

/// Sentinel value meaning "stop codon" in a codon table.
pub const STOP_CODON: u8 = 255;

/// A codon table mapping 64 codons to amino acid digital codes.
/// Codons are encoded as (base1 * 16 + base2 * 4 + base3) using DNA digital
/// codes (A=0, C=1, G=2, T=3).
pub const GeneticCode = struct {
    /// NCBI translation table ID.
    id: u8,
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

    /// Translate a single codon, resolving degenerate nucleotides using the
    /// DNA alphabet's degeneracy bitmasks.  If every possible unambiguous
    /// expansion of the codon maps to the same amino acid (or stop), that
    /// result is returned.  Otherwise returns null.
    pub fn translateAmbiguous(self: GeneticCode, abc: *const Alphabet, c1: u8, c2: u8, c3: u8) ?u8 {
        // Fast path: canonical bases
        if (c1 < 4 and c2 < 4 and c3 < 4) {
            const idx = @as(usize, c1) * 16 + @as(usize, c2) * 4 + @as(usize, c3);
            return self.codon_table[idx];
        }

        // Reject codes outside the alphabet's symbol range
        if (c1 >= abc.kp or c2 >= abc.kp or c3 >= abc.kp) return null;

        const mask1 = abc.degen[c1];
        const mask2 = abc.degen[c2];
        const mask3 = abc.degen[c3];

        // Zero mask means gap/nonresidue/missing — cannot translate
        if (mask1 == 0 or mask2 == 0 or mask3 == 0) return null;

        var consensus: ?u8 = null;

        var b1: u2 = 0;
        while (true) : (b1 += 1) {
            if (mask1 & (@as(u32, 1) << b1) != 0) {
                var b2: u2 = 0;
                while (true) : (b2 += 1) {
                    if (mask2 & (@as(u32, 1) << b2) != 0) {
                        var b3: u2 = 0;
                        while (true) : (b3 += 1) {
                            if (mask3 & (@as(u32, 1) << b3) != 0) {
                                const idx = @as(usize, b1) * 16 + @as(usize, b2) * 4 + @as(usize, b3);
                                const aa = self.codon_table[idx];
                                if (consensus) |prev| {
                                    if (prev != aa) return null;
                                } else {
                                    consensus = aa;
                                }
                            }
                            if (b3 == 3) break;
                        }
                    }
                    if (b2 == 3) break;
                }
            }
            if (b1 == 3) break;
        }

        return consensus;
    }

    /// Return true if the codon formed by the three DNA digital codes is a
    /// start codon.  For degenerate bases, returns true only if ALL
    /// possible expansions are start codons (matching Easel's esl_gencode_IsInitiator).
    pub fn isStart(self: GeneticCode, c1: u8, c2: u8, c3: u8) bool {
        if (c1 >= 4 or c2 >= 4 or c3 >= 4) return self.isStartAmbiguous(c1, c2, c3);
        const idx = @as(usize, c1) * 16 + @as(usize, c2) * 4 + @as(usize, c3);
        return self.is_start[idx];
    }

    /// Check if a degenerate codon is a start codon by expanding all possible
    /// unambiguous codons and verifying they are all initiators.
    fn isStartAmbiguous(self: GeneticCode, c1: u8, c2: u8, c3: u8) bool {
        const dna_abc = &@import("alphabet.zig").dna;
        if (c1 >= dna_abc.kp or c2 >= dna_abc.kp or c3 >= dna_abc.kp) return false;

        const mask1 = dna_abc.degen[c1];
        const mask2 = dna_abc.degen[c2];
        const mask3 = dna_abc.degen[c3];
        if (mask1 == 0 or mask2 == 0 or mask3 == 0) return false;

        var b1: u2 = 0;
        while (true) : (b1 += 1) {
            if (mask1 & (@as(u32, 1) << b1) != 0) {
                var b2: u2 = 0;
                while (true) : (b2 += 1) {
                    if (mask2 & (@as(u32, 1) << b2) != 0) {
                        var b3: u2 = 0;
                        while (true) : (b3 += 1) {
                            if (mask3 & (@as(u32, 1) << b3) != 0) {
                                const idx = @as(usize, b1) * 16 + @as(usize, b2) * 4 + @as(usize, b3);
                                if (!self.is_start[idx]) return false;
                            }
                            if (b3 == 3) break;
                        }
                    }
                    if (b2 == 3) break;
                }
            }
            if (b1 == 3) break;
        }
        return true;
    }

    /// Translate a DNA digital sequence to a protein digital sequence.
    /// Reads in frame starting at offset 0.  Degenerate bases are resolved
    /// via translateAmbiguous; unresolvable codons become the unknown amino
    /// acid (X). Stop codons are included as STOP_CODON markers (the caller
    /// decides how to handle them). Translation continues until fewer than
    /// 3 bases remain.
    ///
    /// The caller owns the returned slice and must free it with the same
    /// allocator.
    pub fn translateSeq(
        self: GeneticCode,
        allocator: Allocator,
        dna_dsq: []const u8,
    ) ![]u8 {
        const abc = &@import("alphabet.zig").dna;
        const n_codons = dna_dsq.len / 3;
        var result: std.ArrayList(u8) = .empty;
        errdefer result.deinit(allocator);
        for (0..n_codons) |i| {
            const c1 = dna_dsq[i * 3];
            const c2 = dna_dsq[i * 3 + 1];
            const c3 = dna_dsq[i * 3 + 2];
            // Try canonical translation first, then ambiguous resolution
            const aa = self.translate(c1, c2, c3) orelse
                self.translateAmbiguous(abc, c1, c2, c3) orelse
                @import("alphabet.zig").amino.unknownCode(); // X for unresolvable
            try result.append(allocator, aa);
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

    /// Write the translation table in a human-readable format showing
    /// codon -> amino acid mappings, grouped by first base.
    /// Example output line: "  AAA -> K (Lys)"
    pub fn writeTable(self: GeneticCode, dest: std.io.AnyWriter) !void {
        const bases = "ACGT";
        const amino_symbols = "ACDEFGHIKLMNPQRSTVWY";
        // Three-letter amino acid names indexed by digital code (0..19).
        const amino_names = [20][]const u8{
            "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
            "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr",
        };

        try dest.print("# Genetic Code: {s} (ID {d})\n", .{ self.name, self.id });

        for (0..4) |b1| {
            for (0..4) |b2| {
                for (0..4) |b3| {
                    const idx = b1 * 16 + b2 * 4 + b3;
                    const aa = self.codon_table[idx];
                    const start_marker: u8 = if (self.is_start[idx]) 'i' else ' ';
                    if (aa == STOP_CODON) {
                        try dest.print("  {c}{c}{c} -> *   (Stop){c}\n", .{
                            bases[b1], bases[b2], bases[b3], start_marker,
                        });
                    } else {
                        try dest.print("  {c}{c}{c} -> {c}   ({s}){c}\n", .{
                            bases[b1], bases[b2], bases[b3],
                            amino_symbols[aa],
                            amino_names[aa],
                            start_marker,
                        });
                    }
                }
            }
        }
    }

    /// Return a mutable copy of this genetic code.
    /// Since GeneticCode is a value type with only arrays and a slice,
    /// a simple copy suffices.
    pub fn clone(self: GeneticCode) GeneticCode {
        return self;
    }

    /// Mark all sense (non-stop) codons as initiators.
    /// Equivalent to Easel's esl_gencode_SetInitiatorAny.
    pub fn setInitiatorAny(self: *GeneticCode) void {
        for (0..64) |idx| {
            self.is_start[idx] = self.codon_table[idx] != STOP_CODON;
        }
    }

    /// Clear all initiators and set only ATG as the start codon.
    /// Equivalent to Easel's esl_gencode_SetInitiatorOnlyAUG.
    pub fn setInitiatorOnlyAUG(self: *GeneticCode) void {
        self.is_start = .{false} ** 64;
        // ATG: A=0, T=3, G=2 -> index = 0*16 + 3*4 + 2 = 14
        self.is_start[codonIndex(0, 3, 2)] = true;
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
    // zeasel amino alphabet canonical symbols in order:
    const amino_symbols = "ACDEFGHIKLMNPQRSTVWY";
    var map: [128]u8 = .{255} ** 128;
    for (amino_symbols, 0..) |c, i| {
        map[c] = @intCast(i);
    }
    return map;
}

/// Map a DNA base ASCII character (A/C/G/T) to its zeasel digital code.
fn baseCharToDigital(c: u8) u8 {
    return switch (c) {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        else => unreachable,
    };
}

/// Build a GeneticCode from Easel-format arrays where amino acid codes use
/// the Easel indexing (0..19 = ACDEFGHIKLMNPQRSTVWY, 27 = stop).
/// The arrays are in zeasel codon order: position i = base1*16 + base2*4 + base3
/// with A=0, C=1, G=2, T=3.
fn buildFromEaselArrays(
    id: u8,
    name: []const u8,
    easel_codons: [64]u8,
    easel_starts: [64]u8,
) GeneticCode {
    var table: [64]u8 = undefined;
    var starts: [64]bool = undefined;

    for (0..64) |i| {
        table[i] = if (easel_codons[i] == 27) STOP_CODON else easel_codons[i];
        starts[i] = easel_starts[i] != 0;
    }

    return GeneticCode{
        .id = id,
        .codon_table = table,
        .is_start = starts,
        .name = name,
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

    // Standard code initiator codons (NCBI translation table 1):
    // ATG, CTG, TTG
    starts[codonIndex(0, 3, 2)] = true; // ATG: A=0, T=3, G=2
    starts[codonIndex(1, 3, 2)] = true; // CTG: C=1, T=3, G=2
    starts[codonIndex(3, 3, 2)] = true; // TTG: T=3, T=3, G=2

    break :blk GeneticCode{
        .id = 1,
        .codon_table = table,
        .is_start = starts,
        .name = "Standard",
    };
};

// ---------------------------------------------------------------------------
// Non-standard NCBI translation tables (2..25)
// ---------------------------------------------------------------------------
// Easel-format arrays: amino acid codes 0..19 = ACDEFGHIKLMNPQRSTVWY, 27 = stop.
// Array order: codon index = base1*16 + base2*4 + base3 (A=0, C=1, G=2, T=3).

pub const vertebrate_mito: GeneticCode = buildFromEaselArrays(2, "Vertebrate Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 27, 15, 27, 15, 10, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const yeast_mito: GeneticCode = buildFromEaselArrays(3, "Yeast Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 10, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 16, 16, 16, 16, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const mold_protozoan_mito: GeneticCode = buildFromEaselArrays(4, "Mold/Protozoan/Coelenterate Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0 });

pub const invertebrate_mito: GeneticCode = buildFromEaselArrays(5, "Invertebrate Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 15, 15, 15, 15, 10, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 });

pub const ciliate_nuclear: GeneticCode = buildFromEaselArrays(6, "Ciliate/Dasycladacean/Hexamita Nuclear", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 13, 19, 13, 19, 15, 15, 15, 15, 27, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const echinoderm_mito: GeneticCode = buildFromEaselArrays(9, "Echinoderm/Flatworm Mitochondrial", .{ 11, 11, 8, 11, 16, 16, 16, 16, 15, 15, 15, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const euplotid_nuclear: GeneticCode = buildFromEaselArrays(10, "Euplotid Nuclear", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 1, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const bacterial: GeneticCode = buildFromEaselArrays(11, "Bacterial/Archaeal/Plant Plastid", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 27, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 });

pub const alt_yeast: GeneticCode = buildFromEaselArrays(12, "Alternative Yeast Nuclear", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 15, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 27, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const ascidian_mito: GeneticCode = buildFromEaselArrays(13, "Ascidian Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 5, 15, 5, 15, 10, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 });

pub const alt_flatworm_mito: GeneticCode = buildFromEaselArrays(14, "Alternative Flatworm Mitochondrial", .{ 11, 11, 8, 11, 16, 16, 16, 16, 15, 15, 15, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 19, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const chlorophycean_mito: GeneticCode = buildFromEaselArrays(16, "Chlorophycean Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 9, 19, 15, 15, 15, 15, 27, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const trematode_mito: GeneticCode = buildFromEaselArrays(21, "Trematode Mitochondrial", .{ 11, 11, 8, 11, 16, 16, 16, 16, 15, 15, 15, 15, 10, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const scenedesmus_mito: GeneticCode = buildFromEaselArrays(22, "Scenedesmus obliquus Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 9, 19, 27, 15, 15, 15, 27, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const thraustochytrium_mito: GeneticCode = buildFromEaselArrays(23, "Thraustochytrium Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 27, 1, 18, 1, 27, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });

pub const pterobranchia_mito: GeneticCode = buildFromEaselArrays(24, "Pterobranchia Mitochondrial", .{ 8, 11, 8, 11, 16, 16, 16, 16, 15, 15, 8, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 18, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 });

pub const gracilibacteria: GeneticCode = buildFromEaselArrays(25, "Candidate Division SR1 and Gracilibacteria", .{ 8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, 13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, 3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, 27, 19, 27, 19, 15, 15, 15, 15, 5, 1, 18, 1, 9, 4, 9, 4 }, .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 });

// ---------------------------------------------------------------------------
// ORF (Open Reading Frame) finding
// ---------------------------------------------------------------------------

/// Represents a single Open Reading Frame found in a DNA sequence.
pub const Orf = struct {
    /// Start position in the original DNA sequence (0-based, inclusive).
    start: usize,
    /// End position in the original DNA sequence (0-based, exclusive, points one past last base of stop codon or sequence end).
    end: usize,
    /// Reading frame (0, 1, or 2).
    frame: u2,
    /// Which strand this ORF was found on.
    strand: Strand,
    /// Translated protein sequence (amino acid digital codes). Caller owns this memory.
    protein: []u8,

    pub const Strand = enum { forward, reverse };
};

/// Translate a single reading frame of a DNA digital sequence.
///
/// Starting at the given frame offset (0, 1, or 2), translates consecutive
/// codons through the entire sequence. Stop codons are included in the output
/// as STOP_CODON markers (255). Translation continues past stop codons.
///
/// The caller owns the returned slice.
pub fn translateFrame(allocator: Allocator, gc: *const GeneticCode, dsq: []const u8, frame: u2) ![]u8 {
    if (dsq.len < @as(usize, frame) + 3) {
        return allocator.alloc(u8, 0);
    }
    const usable = dsq.len - @as(usize, frame);
    const n_codons = usable / 3;
    var result = try std.ArrayList(u8).initCapacity(allocator, n_codons);
    errdefer result.deinit(allocator);

    const dna_abc = &@import("alphabet.zig").dna;
    var pos: usize = frame;
    while (pos + 3 <= dsq.len) {
        const c1 = dsq[pos];
        const c2 = dsq[pos + 1];
        const c3 = dsq[pos + 2];
        const aa = gc.translate(c1, c2, c3) orelse
            gc.translateAmbiguous(dna_abc, c1, c2, c3) orelse
            @import("alphabet.zig").amino.unknownCode(); // X for unresolvable degenerate codons
        try result.append(allocator, aa);
        pos += 3;
    }

    return result.toOwnedSlice(allocator);
}

/// Find all Open Reading Frames (ORFs) in a digital DNA sequence.
///
/// Scans all 3 reading frames on both the forward and reverse strands.
/// An ORF begins at a start codon (or the beginning of the sequence if no
/// start codon precedes the first stop) and ends at a stop codon or the end
/// of the sequence. Only ORFs whose translated protein length is >= min_length
/// amino acids are returned.
///
/// The caller owns the returned slice and each Orf's protein slice.
pub fn findOrfs(
    allocator: Allocator,
    gc: *const GeneticCode,
    abc: *const Alphabet,
    dsq: []const u8,
    min_length: usize,
) ![]Orf {
    var orfs = std.ArrayList(Orf){};
    errdefer {
        for (orfs.items) |orf| allocator.free(orf.protein);
        orfs.deinit(allocator);
    }

    // Forward strand
    try scanStrandOrfs(allocator, gc, dsq, .forward, min_length, &orfs);

    // Reverse strand: make a reverse-complemented copy
    if (abc.complement) |comp| {
        const rc = try allocator.dupe(u8, dsq);
        defer allocator.free(rc);

        // In-place reverse complement
        var i: usize = 0;
        var j: usize = rc.len;
        while (i < j) {
            j -= 1;
            const tmp = comp[rc[i]];
            rc[i] = comp[rc[j]];
            rc[j] = tmp;
            i += 1;
        }

        try scanStrandOrfs(allocator, gc, rc, .reverse, min_length, &orfs);
    }

    return orfs.toOwnedSlice(allocator);
}

/// Scan one strand for ORFs in all 3 reading frames.
fn scanStrandOrfs(
    allocator: Allocator,
    gc: *const GeneticCode,
    dsq: []const u8,
    strand: Orf.Strand,
    min_length: usize,
    orfs: *std.ArrayList(Orf),
) !void {
    for (0..3) |frame_idx| {
        const frame: u2 = @intCast(frame_idx);
        var pos: usize = frame_idx;

        // Accumulate protein residues for the current ORF candidate
        var protein_buf = std.ArrayList(u8){};
        errdefer protein_buf.deinit(allocator);
        var orf_start: usize = pos;
        var seen_start = false;

        while (pos + 3 <= dsq.len) {
            const c1 = dsq[pos];
            const c2 = dsq[pos + 1];
            const c3 = dsq[pos + 2];

            // Handle degenerate bases: try to resolve via ambiguity expansion.
            if (c1 >= 4 or c2 >= 4 or c3 >= 4) {
                const dna_abc = &@import("alphabet.zig").dna;
                const resolved = gc.translateAmbiguous(dna_abc, c1, c2, c3);
                if (resolved) |aa_r| {
                    if (aa_r == STOP_CODON) {
                        // Resolved to stop — end current ORF
                        if (seen_start and protein_buf.items.len >= min_length) {
                            const protein = try protein_buf.toOwnedSlice(allocator);
                            errdefer allocator.free(protein);
                            try orfs.append(allocator, .{
                                .start = orf_start,
                                .end = pos + 3,
                                .frame = frame,
                                .strand = strand,
                                .protein = protein,
                            });
                        } else {
                            protein_buf.clearRetainingCapacity();
                        }
                        seen_start = false;
                        pos += 3;
                        continue;
                    }
                    // Resolved to an amino acid
                    if (seen_start) {
                        try protein_buf.append(allocator, aa_r);
                    }
                } else {
                    // Unresolvable — use X (unknown amino acid) if in ORF
                    if (seen_start) {
                        try protein_buf.append(allocator, @import("alphabet.zig").amino.unknownCode());
                    }
                }
                pos += 3;
                continue;
            }

            const idx = @as(usize, c1) * 16 + @as(usize, c2) * 4 + @as(usize, c3);
            const aa = gc.codon_table[idx];

            if (aa == STOP_CODON) {
                // End of ORF
                if (seen_start and protein_buf.items.len >= min_length) {
                    const protein = try protein_buf.toOwnedSlice(allocator);
                    errdefer allocator.free(protein);
                    try orfs.append(allocator, .{
                        .start = orf_start,
                        .end = pos + 3,
                        .frame = frame,
                        .strand = strand,
                        .protein = protein,
                    });
                } else {
                    protein_buf.clearRetainingCapacity();
                }
                // Reset for next ORF
                seen_start = false;
                pos += 3;
                continue;
            }

            // Check for start codon
            if (!seen_start and gc.is_start[idx]) {
                seen_start = true;
                orf_start = pos;
                protein_buf.clearRetainingCapacity();
                // Start codons translate to Met (M=10)
                try protein_buf.append(allocator, 10);
                pos += 3;
                continue;
            }

            // Regular codon
            if (seen_start) {
                try protein_buf.append(allocator, aa);
            }
            pos += 3;
        }

        // End of sequence: emit trailing ORF if long enough
        if (seen_start and protein_buf.items.len >= min_length) {
            const protein = try protein_buf.toOwnedSlice(allocator);
            errdefer allocator.free(protein);
            try orfs.append(allocator, .{
                .start = orf_start,
                .end = dsq.len,
                .frame = frame,
                .strand = strand,
                .protein = protein,
            });
        } else {
            protein_buf.deinit(allocator);
        }
    }
}

// ---------------------------------------------------------------------------
// Lookup by NCBI translation table ID
// ---------------------------------------------------------------------------

/// All supported genetic codes, indexed for lookup.
const all_codes = [_]*const GeneticCode{
    &standard,
    &vertebrate_mito,
    &yeast_mito,
    &mold_protozoan_mito,
    &invertebrate_mito,
    &ciliate_nuclear,
    &echinoderm_mito,
    &euplotid_nuclear,
    &bacterial,
    &alt_yeast,
    &ascidian_mito,
    &alt_flatworm_mito,
    &chlorophycean_mito,
    &trematode_mito,
    &scenedesmus_mito,
    &thraustochytrium_mito,
    &pterobranchia_mito,
    &gracilibacteria,
};

/// Look up a genetic code by its NCBI translation table ID.
/// Returns null if the ID is not supported.
pub fn byId(id: u8) ?*const GeneticCode {
    for (&all_codes) |gc| {
        if (gc.id == id) return gc;
    }
    return null;
}

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

test "translateSeq includes stop codon ATGTAAATG -> [M, STOP, M]" {
    // ATG=M(10), TAA=stop(255), ATG=M(10) — stop codons are now included
    // Digital: A=0,T=3,G=2, T=3,A=0,A=0, A=0,T=3,G=2
    const dna_dsq = [_]u8{ 0, 3, 2, 3, 0, 0, 0, 3, 2 };
    const allocator = std.testing.allocator;
    const protein = try standard.translateSeq(allocator, &dna_dsq);
    defer allocator.free(protein);

    try std.testing.expectEqualSlices(u8, &[_]u8{ 10, STOP_CODON, 10 }, protein);
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

test "byId returns correct tables" {
    // Standard code
    const std_code = byId(1);
    try std.testing.expect(std_code != null);
    try std.testing.expectEqual(@as(u8, 1), std_code.?.id);
    try std.testing.expectEqualStrings("Standard", std_code.?.name);

    // Vertebrate mito
    const vmito = byId(2);
    try std.testing.expect(vmito != null);
    try std.testing.expectEqual(@as(u8, 2), vmito.?.id);

    // Bacterial
    const bact = byId(11);
    try std.testing.expect(bact != null);
    try std.testing.expectEqual(@as(u8, 11), bact.?.id);

    // Non-existent table
    try std.testing.expectEqual(@as(?*const GeneticCode, null), byId(0));
    try std.testing.expectEqual(@as(?*const GeneticCode, null), byId(7));
    try std.testing.expectEqual(@as(?*const GeneticCode, null), byId(99));
}

test "vertebrate mito: TGA -> W (18) instead of stop" {
    // In the standard code, TGA (T=3, G=2, A=0) is a stop codon.
    // In vertebrate mitochondrial code, TGA encodes Trp (W=18).
    const std_aa = standard.translate(3, 2, 0);
    try std.testing.expectEqual(@as(?u8, null), std_aa); // stop in standard

    const vmito_aa = vertebrate_mito.translate(3, 2, 0);
    try std.testing.expectEqual(@as(?u8, 18), vmito_aa); // W in vert mito

    // AGA -> stop in vertebrate mito (not Arg as in standard)
    // AGA: A=0, G=2, A=0
    const vmito_aga = vertebrate_mito.translate(0, 2, 0);
    try std.testing.expectEqual(@as(?u8, null), vmito_aga); // stop

    const std_aga = standard.translate(0, 2, 0);
    try std.testing.expectEqual(@as(?u8, 14), std_aga); // R=14 in standard
}

test "translateAmbiguous with canonical bases" {
    const abc = &@import("alphabet.zig").dna;

    // Canonical bases should work the same as translate
    const aa = standard.translateAmbiguous(abc, 0, 3, 2); // ATG -> M
    try std.testing.expectEqual(@as(?u8, 10), aa);
}

test "translateAmbiguous with degenerate bases" {
    const abc = &@import("alphabet.zig").dna;

    // GCN (N=any) -> all four GCx codons encode Ala (A=0)
    // N is code 15 in DNA alphabet
    const aa = standard.translateAmbiguous(abc, 2, 1, 15); // GC + N
    try std.testing.expectEqual(@as(?u8, 0), aa); // A=0

    // ATR (R=A|G, code 5): ATA=I(7), ATG=M(10) -> ambiguous
    const ambig = standard.translateAmbiguous(abc, 0, 3, 5); // AT + R
    try std.testing.expectEqual(@as(?u8, null), ambig);

    // GTN (N=any) -> all four GTx codons encode Val (V=17)
    const val = standard.translateAmbiguous(abc, 2, 3, 15); // GT + N
    try std.testing.expectEqual(@as(?u8, 17), val); // V=17
}

test "translateAmbiguous with stop codons" {
    const abc = &@import("alphabet.zig").dna;

    // TAA is stop, TAR (R=A|G) includes TAA=stop and TAG=stop -> both stop -> STOP_CODON
    const stop = standard.translateAmbiguous(abc, 3, 0, 5); // TA + R
    try std.testing.expectEqual(@as(?u8, STOP_CODON), stop);

    // TAY (Y=C|T): TAC=Y(19), TAT=Y(19) -> Y
    const tyr = standard.translateAmbiguous(abc, 3, 0, 6); // TA + Y
    try std.testing.expectEqual(@as(?u8, 19), tyr);
}

test "writeTable: produces 64 codon lines" {
    const allocator = std.testing.allocator;
    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);
    try standard.writeTable(buf.writer(allocator).any());

    // Count lines starting with "  " (codon lines)
    var codon_lines: usize = 0;
    var lines = std.mem.splitScalar(u8, buf.items, '\n');
    while (lines.next()) |line| {
        if (std.mem.startsWith(u8, line, "  ")) codon_lines += 1;
    }
    try std.testing.expectEqual(@as(usize, 64), codon_lines);

    // Should contain the header with the code name
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "Standard") != null);

    // Should contain ATG -> M mapping
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "ATG -> M") != null);
}

test "translateFrame: frame 0 translates all codons including stops" {
    const allocator = std.testing.allocator;
    // ATG TAA GCT = M(10) STOP(255) A(0)
    // Digital: A=0 T=3 G=2 T=3 A=0 A=0 G=2 C=1 T=3
    const dsq = [_]u8{ 0, 3, 2, 3, 0, 0, 2, 1, 3 };
    const result = try translateFrame(allocator, &standard, &dsq, 0);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 3), result.len);
    try std.testing.expectEqual(@as(u8, 10), result[0]); // M
    try std.testing.expectEqual(@as(u8, STOP_CODON), result[1]); // stop
    try std.testing.expectEqual(@as(u8, 0), result[2]); // A
}

test "translateFrame: frame 1 shifts start by 1" {
    const allocator = std.testing.allocator;
    // Sequence: x ATG GCT ...
    // Digital: A=0  A=0 T=3 G=2  G=2 C=1 T=3
    const dsq = [_]u8{ 0, 0, 3, 2, 2, 1, 3 };
    const result = try translateFrame(allocator, &standard, &dsq, 1);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 2), result.len);
    try std.testing.expectEqual(@as(u8, 10), result[0]); // M (ATG)
    try std.testing.expectEqual(@as(u8, 0), result[1]); // A (GCT)
}

test "translateFrame: empty for short sequence" {
    const allocator = std.testing.allocator;
    const dsq = [_]u8{ 0, 3 };
    const result = try translateFrame(allocator, &standard, &dsq, 0);
    defer allocator.free(result);
    try std.testing.expectEqual(@as(usize, 0), result.len);
}

test "findOrfs: single forward ORF with ATG...TAA" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // ATG GCT TAT TAA = M A Y stop
    // Digital: A=0 T=3 G=2 G=2 C=1 T=3 T=3 A=0 T=3 T=3 A=0 A=0
    const dsq = [_]u8{ 0, 3, 2, 2, 1, 3, 3, 0, 3, 3, 0, 0 };
    const orfs = try findOrfs(allocator, &standard, abc, &dsq, 1);
    defer {
        for (orfs) |orf| allocator.free(orf.protein);
        allocator.free(orfs);
    }

    // Should find at least one forward ORF
    var found_forward = false;
    for (orfs) |orf| {
        if (orf.strand == .forward and orf.frame == 0) {
            found_forward = true;
            try std.testing.expectEqual(@as(usize, 0), orf.start);
            try std.testing.expectEqual(@as(usize, 12), orf.end); // past stop codon
            try std.testing.expectEqual(@as(usize, 3), orf.protein.len);
            try std.testing.expectEqual(@as(u8, 10), orf.protein[0]); // M
            try std.testing.expectEqual(@as(u8, 0), orf.protein[1]); // A
            try std.testing.expectEqual(@as(u8, 19), orf.protein[2]); // Y
        }
    }
    try std.testing.expect(found_forward);
}

test "findOrfs: min_length filters short ORFs" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // ATG TAA = M stop (protein length 1)
    // Digital: A=0 T=3 G=2 T=3 A=0 A=0
    const dsq = [_]u8{ 0, 3, 2, 3, 0, 0 };

    // min_length=2 should filter it out (forward frame 0 only has 1 aa)
    const orfs = try findOrfs(allocator, &standard, abc, &dsq, 2);
    defer {
        for (orfs) |orf| allocator.free(orf.protein);
        allocator.free(orfs);
    }

    var has_frame0_fwd = false;
    for (orfs) |orf| {
        if (orf.strand == .forward and orf.frame == 0) has_frame0_fwd = true;
    }
    try std.testing.expect(!has_frame0_fwd);
}

test "findOrfs: reverse strand ORFs are found" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // Forward: TTA TAG GCC CAT (no start codon in frame 0)
    // Reverse complement: ATG GGC CTA TAA -> M G L stop
    // Digital forward: T=3 T=3 A=0 T=3 A=0 G=2 G=2 C=1 C=1 C=1 A=0 T=3
    const dsq = [_]u8{ 3, 3, 0, 3, 0, 2, 2, 1, 1, 1, 0, 3 };
    const orfs = try findOrfs(allocator, &standard, abc, &dsq, 1);
    defer {
        for (orfs) |orf| allocator.free(orf.protein);
        allocator.free(orfs);
    }

    var found_reverse = false;
    for (orfs) |orf| {
        if (orf.strand == .reverse and orf.frame == 0) {
            found_reverse = true;
            try std.testing.expectEqual(@as(usize, 3), orf.protein.len);
            try std.testing.expectEqual(@as(u8, 10), orf.protein[0]); // M
        }
    }
    try std.testing.expect(found_reverse);
}

test "findOrfs: ORF at end of sequence without stop codon" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // ATG GCT GCT (no stop codon) -> M A A (length 3)
    // Digital: A=0 T=3 G=2 G=2 C=1 T=3 G=2 C=1 T=3
    const dsq = [_]u8{ 0, 3, 2, 2, 1, 3, 2, 1, 3 };
    const orfs = try findOrfs(allocator, &standard, abc, &dsq, 1);
    defer {
        for (orfs) |orf| allocator.free(orf.protein);
        allocator.free(orfs);
    }

    var found = false;
    for (orfs) |orf| {
        if (orf.strand == .forward and orf.frame == 0) {
            found = true;
            try std.testing.expectEqual(@as(usize, 0), orf.start);
            try std.testing.expectEqual(@as(usize, 9), orf.end); // end of sequence
            try std.testing.expectEqual(@as(usize, 3), orf.protein.len);
            try std.testing.expectEqual(@as(u8, 10), orf.protein[0]); // M
        }
    }
    try std.testing.expect(found);
}

test "findOrfs: multiple ORFs in different frames" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // Build a sequence that has ORFs in frame 0 and frame 1
    // Frame 0: ATG GCT TAA ... (M A stop)
    // Frame 1:  TGG CTT AA... (no start in frame 1)
    // We need something starting with ATG in frame 1 too.
    // So: x ATG GCT TAA where x is a single base
    // Frame 0 from pos 0: xAT GGC TTA A -> no start
    // Frame 1 from pos 1: ATG GCT TAA   -> M A stop
    // Let x = G:
    // G=2, A=0, T=3, G=2, G=2, C=1, T=3, T=3, A=0, A=0
    // Frame 0: GAT(D=2) GGC(G=5) TTA(L=9) A -- no ATG start
    // Frame 1: ATG(M=10) GCT(A=0) TAA(stop)
    const dsq = [_]u8{ 2, 0, 3, 2, 2, 1, 3, 3, 0, 0 };
    const orfs = try findOrfs(allocator, &standard, abc, &dsq, 1);
    defer {
        for (orfs) |orf| allocator.free(orf.protein);
        allocator.free(orfs);
    }

    var found_frame1 = false;
    for (orfs) |orf| {
        if (orf.strand == .forward and orf.frame == 1) {
            found_frame1 = true;
            try std.testing.expectEqual(@as(u8, 10), orf.protein[0]); // M
        }
    }
    try std.testing.expect(found_frame1);
}

test "findOrfs: start codon translates to Met" {
    const allocator = std.testing.allocator;
    const abc = &@import("alphabet.zig").dna;

    // CTG is a start codon in standard code but normally translates to Leu.
    // In ORF finding, start codons should translate to Met (M=10).
    // CTG GCT TAA
    // Digital: C=1 T=3 G=2 G=2 C=1 T=3 T=3 A=0 A=0
    const dsq = [_]u8{ 1, 3, 2, 2, 1, 3, 3, 0, 0 };
    const orfs = try findOrfs(allocator, &standard, abc, &dsq, 1);
    defer {
        for (orfs) |orf| allocator.free(orf.protein);
        allocator.free(orfs);
    }

    var found = false;
    for (orfs) |orf| {
        if (orf.strand == .forward and orf.frame == 0) {
            found = true;
            try std.testing.expectEqual(@as(u8, 10), orf.protein[0]); // M, not L
        }
    }
    try std.testing.expect(found);
}

test "clone returns independent mutable copy" {
    var code = standard.clone();
    // Mutating the clone should not affect the original (it's a value copy).
    code.is_start[0] = !standard.is_start[0];
    try std.testing.expect(code.is_start[0] != standard.is_start[0]);
}

test "setInitiatorAny: all sense codons are starts" {
    var code = standard.clone();
    code.setInitiatorAny();
    // CTG (which is normally a start) should still be start
    try std.testing.expect(code.isStart(1, 3, 2));
    // AAA (Lys, not normally a start) should now be start
    try std.testing.expect(code.isStart(0, 0, 0));
    // TAA (stop) should NOT be start
    try std.testing.expect(!code.isStart(3, 0, 0));
    // TAG (stop) should NOT be start
    try std.testing.expect(!code.isStart(3, 0, 2));
    // TGA (stop) should NOT be start
    try std.testing.expect(!code.isStart(3, 2, 0));
}

test "setInitiatorOnlyAUG: only ATG is start" {
    var code = standard.clone();
    code.setInitiatorOnlyAUG();
    // ATG should be start
    try std.testing.expect(code.isStart(0, 3, 2));
    // CTG (normally a start) should NOT be start
    try std.testing.expect(!code.isStart(1, 3, 2));
    // TTG (normally a start) should NOT be start
    try std.testing.expect(!code.isStart(3, 3, 2));
}

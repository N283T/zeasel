// Biological sequence with digital encoding and metadata.

const std = @import("std");
const Allocator = std.mem.Allocator;
const alphabet_mod = @import("alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;

pub const Sequence = struct {
    name: []const u8,
    accession: ?[]const u8,
    description: ?[]const u8,
    taxonomy_id: ?i32,
    dsq: []u8,
    secondary_structure: ?[]const u8,
    source: ?Source,
    abc: *const Alphabet,
    allocator: Allocator,

    pub const Source = struct {
        name: []const u8,
        start: i64, // 1-indexed; start > end means reverse strand
        end: i64,
        full_length: i64,
    };

    /// Create a new sequence from digital data (copies dsq and name).
    pub fn init(allocator: Allocator, abc: *const Alphabet, name: []const u8, dsq: []const u8) !Sequence {
        const name_copy = try allocator.dupe(u8, name);
        errdefer allocator.free(name_copy);
        const dsq_copy = try allocator.dupe(u8, dsq);
        errdefer allocator.free(dsq_copy);

        return Sequence{
            .name = name_copy,
            .accession = null,
            .description = null,
            .taxonomy_id = null,
            .dsq = dsq_copy,
            .secondary_structure = null,
            .source = null,
            .abc = abc,
            .allocator = allocator,
        };
    }

    /// Create from text sequence (digitizes internally, copies name).
    pub fn fromText(allocator: Allocator, abc: *const Alphabet, name: []const u8, text: []const u8) !Sequence {
        const dsq = try abc.digitize(allocator, text);
        errdefer allocator.free(dsq);
        const name_copy = try allocator.dupe(u8, name);
        errdefer allocator.free(name_copy);

        return Sequence{
            .name = name_copy,
            .accession = null,
            .description = null,
            .taxonomy_id = null,
            .dsq = dsq,
            .secondary_structure = null,
            .source = null,
            .abc = abc,
            .allocator = allocator,
        };
    }

    /// Length of the digital sequence.
    pub fn len(self: Sequence) usize {
        return self.dsq.len;
    }

    /// Deep copy.
    pub fn clone(self: Sequence) !Sequence {
        const name_copy = try self.allocator.dupe(u8, self.name);
        errdefer self.allocator.free(name_copy);
        const dsq_copy = try self.allocator.dupe(u8, self.dsq);
        errdefer self.allocator.free(dsq_copy);

        var accession_copy: ?[]const u8 = null;
        if (self.accession) |acc| {
            accession_copy = try self.allocator.dupe(u8, acc);
        }
        errdefer if (accession_copy) |acc| self.allocator.free(acc);

        var description_copy: ?[]const u8 = null;
        if (self.description) |desc| {
            description_copy = try self.allocator.dupe(u8, desc);
        }
        errdefer if (description_copy) |desc| self.allocator.free(desc);

        var ss_copy: ?[]const u8 = null;
        if (self.secondary_structure) |ss| {
            ss_copy = try self.allocator.dupe(u8, ss);
        }
        errdefer if (ss_copy) |ss| self.allocator.free(ss);

        var source_copy: ?Source = null;
        if (self.source) |src| {
            const src_name = try self.allocator.dupe(u8, src.name);
            source_copy = Source{
                .name = src_name,
                .start = src.start,
                .end = src.end,
                .full_length = src.full_length,
            };
        }
        errdefer if (source_copy) |src| self.allocator.free(src.name);

        return Sequence{
            .name = name_copy,
            .accession = accession_copy,
            .description = description_copy,
            .taxonomy_id = self.taxonomy_id,
            .dsq = dsq_copy,
            .secondary_structure = ss_copy,
            .source = source_copy,
            .abc = self.abc,
            .allocator = self.allocator,
        };
    }

    /// Get text representation (caller owns returned memory).
    pub fn toText(self: Sequence) ![]u8 {
        return self.abc.textize(self.allocator, self.dsq);
    }

    /// Reverse complement in place (DNA/RNA only).
    /// Swaps source start/end coordinates and invalidates secondary structure.
    pub fn reverseComplement(self: *Sequence) !void {
        try self.abc.reverseComplement(self.dsq);
        // Swap source coordinates to reflect reverse strand
        if (self.source) |src| {
            self.source = Source{
                .name = src.name,
                .start = src.end,
                .end = src.start,
                .full_length = src.full_length,
            };
        }
        // Secondary structure is invalidated by reverse complement
        if (self.secondary_structure) |ss| {
            self.allocator.free(ss);
            self.secondary_structure = null;
        }
    }

    /// Extract subsequence [start..end) (0-indexed, half-open). Returns new Sequence with Source tracking.
    pub fn subseq(self: Sequence, start: usize, end: usize) !Sequence {
        if (start > end or end > self.dsq.len) return error.OutOfBounds;

        const dsq_copy = try self.allocator.dupe(u8, self.dsq[start..end]);
        errdefer self.allocator.free(dsq_copy);
        const src_name = try self.allocator.dupe(u8, self.name);
        errdefer self.allocator.free(src_name);

        const src = Source{
            .name = src_name,
            .start = @intCast(start + 1),
            .end = @intCast(end),
            .full_length = @intCast(self.dsq.len),
        };

        const seq_name = try self.allocator.dupe(u8, self.name);
        errdefer self.allocator.free(seq_name);

        return Sequence{
            .name = seq_name,
            .accession = null,
            .description = null,
            .taxonomy_id = null,
            .dsq = dsq_copy,
            .secondary_structure = null,
            .source = src,
            .abc = self.abc,
            .allocator = self.allocator,
        };
    }

    // --- Metadata setters ---

    pub fn setName(self: *Sequence, s: []const u8) !void {
        const new_name = try self.allocator.dupe(u8, s);
        self.allocator.free(self.name);
        self.name = new_name;
    }

    pub fn setAccession(self: *Sequence, s: []const u8) !void {
        if (self.accession) |old| self.allocator.free(old);
        self.accession = try self.allocator.dupe(u8, s);
    }

    pub fn setDescription(self: *Sequence, s: []const u8) !void {
        if (self.description) |old| self.allocator.free(old);
        self.description = try self.allocator.dupe(u8, s);
    }

    pub fn setSource(self: *Sequence, name: []const u8, start: i64, end: i64, full_length: i64) !void {
        if (self.source) |old| self.allocator.free(old.name);
        const src_name = try self.allocator.dupe(u8, name);
        self.source = Source{
            .name = src_name,
            .start = start,
            .end = end,
            .full_length = full_length,
        };
    }

    /// Append additional residues to an existing sequence by digitizing the
    /// text and extending dsq. Enables incremental building from multiple lines.
    pub fn appendDigital(self: *Sequence, text: []const u8) !void {
        const new_codes = try self.abc.digitize(self.allocator, text);
        errdefer self.allocator.free(new_codes);

        const old_len = self.dsq.len;
        const combined = try self.allocator.alloc(u8, old_len + new_codes.len);
        @memcpy(combined[0..old_len], self.dsq);
        @memcpy(combined[old_len..], new_codes);
        self.allocator.free(self.dsq);
        self.allocator.free(new_codes);
        self.dsq = combined;
    }

    // --- Utilities ---

    /// Compute a checksum of the digital sequence data using Jenkins one-at-a-time hash.
    /// Compatible with Easel's esl_sq_Checksum().
    pub fn checksum(self: Sequence) u32 {
        var val: u32 = 0;
        for (self.dsq) |code| {
            val +%= @as(u32, code);
            val +%= val << 10;
            val ^= val >> 6;
        }
        val +%= val << 3;
        val ^= val >> 11;
        val +%= val << 15;
        return val;
    }

    /// Count residue frequencies. Returns array of length abc.K (canonical residues only).
    /// Degenerate residues are distributed equally across their constituent canonical
    /// residues using the degeneracy bitmask. Gaps, missing data, and nonresidue
    /// symbols are ignored. Caller owns the returned slice.
    pub fn countResidues(self: Sequence) ![]f64 {
        const k: usize = @intCast(self.abc.k);
        const counts = try self.allocator.alloc(f64, k);
        @memset(counts, 0.0);
        for (self.dsq) |code| {
            if (self.abc.isCanonical(code)) {
                // Canonical residue: count directly
                counts[code] += 1.0;
            } else if (self.abc.isGap(code) or self.abc.isMissing(code) or self.abc.isNonresidue(code)) {
                // Skip gaps, missing data, nonresidue
                continue;
            } else if (code < self.abc.kp) {
                // Degenerate or unknown: distribute weight equally across constituents
                const ndegen = self.abc.ndegen[code];
                if (ndegen == 0) continue;
                const mask = self.abc.degen[code];
                const wt: f64 = 1.0 / @as(f64, @floatFromInt(ndegen));
                for (0..k) |y| {
                    if (mask & (@as(u32, 1) << @intCast(y)) != 0) {
                        counts[y] += wt;
                    }
                }
            }
        }
        return counts;
    }

    /// Convert all degenerate residue codes to the unknown symbol (X/N).
    /// Modifies the sequence in place.
    pub fn convertDegen2X(self: *Sequence) void {
        const unknown = self.abc.unknownCode();
        for (self.dsq) |*code| {
            if (self.abc.isDegenerate(code.*)) {
                code.* = unknown;
            }
        }
    }

    /// Guess the alphabet type from the sequence content.
    pub fn guessAlphabet(self: Sequence) ?alphabet_mod.AlphabetType {
        const text = self.abc.textize(self.allocator, self.dsq) catch return null;
        defer self.allocator.free(text);
        return alphabet_mod.guessType(text);
    }

    /// Compare two sequences for equality (name and digital sequence data).
    pub fn eql(self: Sequence, other: Sequence) bool {
        if (!std.mem.eql(u8, self.name, other.name)) return false;
        if (!std.mem.eql(u8, self.dsq, other.dsq)) return false;
        return true;
    }

    /// Validate internal consistency.
    pub fn isValid(self: Sequence) bool {
        if (self.name.len == 0) return false;
        for (self.dsq) |code| {
            if (code >= self.abc.kp) return false;
        }
        return true;
    }

    /// Free all owned memory.
    pub fn deinit(self: *Sequence) void {
        self.allocator.free(self.name);
        self.allocator.free(self.dsq);
        if (self.accession) |acc| self.allocator.free(acc);
        if (self.description) |desc| self.allocator.free(desc);
        if (self.secondary_structure) |ss| self.allocator.free(ss);
        if (self.source) |src| self.allocator.free(src.name);
    }
};

/// A block of sequences for batch processing in threaded pipelines.
/// Sequences are read in blocks and handed to worker threads.
pub const SequenceBlock = struct {
    seqs: []Sequence,
    count: usize, // Number of valid sequences (may be less than seqs.len for last block)
    allocator: Allocator,

    /// Create a block with pre-allocated capacity.
    pub fn initCapacity(allocator: Allocator, capacity: usize) !SequenceBlock {
        const seqs = try allocator.alloc(Sequence, capacity);
        return SequenceBlock{
            .seqs = seqs,
            .count = 0,
            .allocator = allocator,
        };
    }

    /// Add a sequence to the block. Takes ownership of the sequence.
    pub fn add(self: *SequenceBlock, seq: Sequence) !void {
        if (self.count >= self.seqs.len) {
            // Grow the block
            const new_cap = if (self.seqs.len == 0) 8 else self.seqs.len * 2;
            const new_seqs = try self.allocator.alloc(Sequence, new_cap);
            @memcpy(new_seqs[0..self.count], self.seqs[0..self.count]);
            self.allocator.free(self.seqs);
            self.seqs = new_seqs;
        }
        self.seqs[self.count] = seq;
        self.count += 1;
    }

    /// Reset the block for reuse, freeing all contained sequences.
    pub fn reset(self: *SequenceBlock) void {
        for (0..self.count) |i| {
            self.seqs[i].deinit();
        }
        self.count = 0;
    }

    /// Get the valid sequences in this block.
    pub fn items(self: SequenceBlock) []Sequence {
        return self.seqs[0..self.count];
    }

    /// Free the block and all contained sequences.
    pub fn deinit(self: *SequenceBlock) void {
        for (0..self.count) |i| {
            self.seqs[i].deinit();
        }
        self.allocator.free(self.seqs);
    }
};

// --- Tests ---

test "init: basic properties" {
    const allocator = std.testing.allocator;
    const dsq = &[_]u8{ 0, 1, 2, 3 };
    var seq = try Sequence.init(allocator, &alphabet_mod.dna, "seq1", dsq);
    defer seq.deinit();

    try std.testing.expectEqualStrings("seq1", seq.name);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, seq.dsq);
    try std.testing.expectEqual(@as(usize, 4), seq.len());
}

test "init: optional fields are null" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.init(allocator, &alphabet_mod.dna, "seq1", &[_]u8{0});
    defer seq.deinit();

    try std.testing.expectEqual(@as(?[]const u8, null), seq.accession);
    try std.testing.expectEqual(@as(?[]const u8, null), seq.description);
    try std.testing.expectEqual(@as(?[]const u8, null), seq.secondary_structure);
    try std.testing.expectEqual(@as(?i32, null), seq.taxonomy_id);
    try std.testing.expectEqual(@as(?Sequence.Source, null), seq.source);
}

test "fromText: ACGT digitizes to [0,1,2,3]" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, seq.dsq);
    try std.testing.expectEqual(@as(usize, 4), seq.len());
}

test "fromText: invalid text returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(
        error.InvalidCharacter,
        Sequence.fromText(allocator, &alphabet_mod.dna, "bad", "ACGZ"),
    );
}

test "deinit: no memory leaks" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.init(allocator, &alphabet_mod.dna, "seq1", &[_]u8{ 0, 1, 2 });
    seq.deinit();
    // std.testing.allocator detects leaks automatically
}

test "clone: deep copy independence" {
    const allocator = std.testing.allocator;
    var original = try Sequence.fromText(allocator, &alphabet_mod.dna, "orig", "ACGT");
    defer original.deinit();

    var cloned = try original.clone();
    defer cloned.deinit();

    // Modify original dsq
    original.dsq[0] = 3;

    // Clone should be unaffected
    try std.testing.expectEqual(@as(u8, 0), cloned.dsq[0]);
    try std.testing.expectEqualStrings("orig", cloned.name);
}

test "toText: round trip ACGT" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    const text = try seq.toText();
    defer allocator.free(text);

    try std.testing.expectEqualStrings("ACGT", text);
}

test "reverseComplement: AACG -> CGTT" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "AACG");
    defer seq.deinit();

    try seq.reverseComplement();

    const text = try seq.toText();
    defer allocator.free(text);

    try std.testing.expectEqualStrings("CGTT", text);
}

test "reverseComplement: amino acid returns error" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.amino, "prot", "ACDEF");
    defer seq.deinit();

    try std.testing.expectError(error.NoComplement, seq.reverseComplement());
}

test "reverseComplement: swaps source start and end" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "AACG");
    defer seq.deinit();

    try seq.setSource("chr1", 100, 200, 50000);
    try seq.reverseComplement();

    // start and end should be swapped
    try std.testing.expectEqual(@as(i64, 200), seq.source.?.start);
    try std.testing.expectEqual(@as(i64, 100), seq.source.?.end);
    try std.testing.expectEqual(@as(i64, 50000), seq.source.?.full_length);
}

test "reverseComplement: invalidates secondary structure" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "AACG");
    defer seq.deinit();

    // Set secondary structure
    seq.secondary_structure = try allocator.dupe(u8, "<<>>");

    try seq.reverseComplement();

    // Secondary structure should be null after revcomp
    try std.testing.expectEqual(@as(?[]const u8, null), seq.secondary_structure);
}

test "subseq: basic extraction with source tracking" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "full", "ACGTACGT");
    defer seq.deinit();

    var sub = try seq.subseq(2, 6);
    defer sub.deinit();

    // dsq[2..6] = G,T,A,C
    const text = try sub.toText();
    defer allocator.free(text);
    try std.testing.expectEqualStrings("GTAC", text);

    // Source tracking: start is 1-indexed so start=3, end=6, full_length=8
    try std.testing.expect(sub.source != null);
    try std.testing.expectEqual(@as(i64, 3), sub.source.?.start);
    try std.testing.expectEqual(@as(i64, 6), sub.source.?.end);
    try std.testing.expectEqual(@as(i64, 8), sub.source.?.full_length);
    try std.testing.expectEqualStrings("full", sub.source.?.name);
}

test "subseq: out of bounds returns error" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    try std.testing.expectError(error.OutOfBounds, seq.subseq(3, 6));
    try std.testing.expectError(error.OutOfBounds, seq.subseq(3, 2));
}

// --- New tests for added functionality ---

test "setName: replaces name" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "old", "ACGT");
    defer seq.deinit();

    try seq.setName("new_name");
    try std.testing.expectEqualStrings("new_name", seq.name);
}

test "setAccession and setDescription" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    try seq.setAccession("P12345");
    try std.testing.expectEqualStrings("P12345", seq.accession.?);

    try seq.setDescription("A test sequence");
    try std.testing.expectEqualStrings("A test sequence", seq.description.?);

    // Replace existing
    try seq.setAccession("Q99999");
    try std.testing.expectEqualStrings("Q99999", seq.accession.?);
}

test "setSource" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    try seq.setSource("chr1", 100, 200, 50000);
    try std.testing.expectEqualStrings("chr1", seq.source.?.name);
    try std.testing.expectEqual(@as(i64, 100), seq.source.?.start);
    try std.testing.expectEqual(@as(i64, 200), seq.source.?.end);
}

test "checksum: deterministic" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    const c1 = seq.checksum();
    const c2 = seq.checksum();
    try std.testing.expectEqual(c1, c2);
}

test "checksum: different sequences differ" {
    const allocator = std.testing.allocator;
    var s1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "a", "ACGT");
    defer s1.deinit();
    var s2 = try Sequence.fromText(allocator, &alphabet_mod.dna, "b", "TGCA");
    defer s2.deinit();

    try std.testing.expect(s1.checksum() != s2.checksum());
}

test "countResidues: AACG -> 2A, 1C, 1G (K-length array)" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "AACG");
    defer seq.deinit();

    const counts = try seq.countResidues();
    defer allocator.free(counts);

    // Array length should be K (4 for DNA), not Kp
    try std.testing.expectEqual(@as(usize, 4), counts.len);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), counts[0], 1e-9); // A
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), counts[1], 1e-9); // C
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), counts[2], 1e-9); // G
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), counts[3], 1e-9); // T
}

test "countResidues: degenerate R distributes to A and G" {
    const allocator = std.testing.allocator;
    // R = A or G, so should distribute 0.5 to each
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "RR");
    defer seq.deinit();

    const counts = try seq.countResidues();
    defer allocator.free(counts);

    try std.testing.expectApproxEqAbs(@as(f64, 1.0), counts[0], 1e-9); // A: 2 * 0.5
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), counts[1], 1e-9); // C
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), counts[2], 1e-9); // G: 2 * 0.5
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), counts[3], 1e-9); // T
}

test "convertDegen2X: R becomes N" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACRG");
    defer seq.deinit();

    seq.convertDegen2X();

    const text = try seq.toText();
    defer allocator.free(text);
    try std.testing.expectEqualStrings("ACNG", text);
}

test "eql: identical sequences" {
    const allocator = std.testing.allocator;
    var s1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer s1.deinit();
    var s2 = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer s2.deinit();

    try std.testing.expect(s1.eql(s2));
}

test "eql: different sequences" {
    const allocator = std.testing.allocator;
    var s1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer s1.deinit();
    var s2 = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq2", "ACGT");
    defer s2.deinit();

    try std.testing.expect(!s1.eql(s2));
}

test "isValid: valid sequence" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    try std.testing.expect(seq.isValid());
}

test "SequenceBlock: basic add and iterate" {
    const allocator = std.testing.allocator;
    var block = try SequenceBlock.initCapacity(allocator, 4);
    defer block.deinit();

    const s1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "s1", "ACGT");
    try block.add(s1);
    const s2 = try Sequence.fromText(allocator, &alphabet_mod.dna, "s2", "TGCA");
    try block.add(s2);

    try std.testing.expectEqual(@as(usize, 2), block.count);
    try std.testing.expectEqualStrings("s1", block.items()[0].name);
    try std.testing.expectEqualStrings("s2", block.items()[1].name);
}

test "SequenceBlock: reset and reuse" {
    const allocator = std.testing.allocator;
    var block = try SequenceBlock.initCapacity(allocator, 4);
    defer block.deinit();

    const s1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "s1", "ACGT");
    try block.add(s1);
    try std.testing.expectEqual(@as(usize, 1), block.count);

    block.reset();
    try std.testing.expectEqual(@as(usize, 0), block.count);

    // Reuse
    const s2 = try Sequence.fromText(allocator, &alphabet_mod.dna, "s2", "TGCA");
    try block.add(s2);
    try std.testing.expectEqual(@as(usize, 1), block.count);
    try std.testing.expectEqualStrings("s2", block.items()[0].name);
}

test "SequenceBlock: auto-grow beyond initial capacity" {
    const allocator = std.testing.allocator;
    var block = try SequenceBlock.initCapacity(allocator, 2);
    defer block.deinit();

    for (0..5) |i| {
        var buf: [8]u8 = undefined;
        const name = std.fmt.bufPrint(&buf, "s{d}", .{i}) catch unreachable;
        const seq = try Sequence.fromText(allocator, &alphabet_mod.dna, name, "ACGT");
        try block.add(seq);
    }

    try std.testing.expectEqual(@as(usize, 5), block.count);
}

test "appendDigital: incremental sequence building" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    try seq.appendDigital("TGCA");

    try std.testing.expectEqual(@as(usize, 8), seq.len());

    const text = try seq.toText();
    defer allocator.free(text);
    try std.testing.expectEqualStrings("ACGTTGCA", text);
}

test "appendDigital: empty append is no-op" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    try seq.appendDigital("");

    try std.testing.expectEqual(@as(usize, 4), seq.len());
}

test "appendDigital: multiple appends" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "AC");
    defer seq.deinit();

    try seq.appendDigital("GT");
    try seq.appendDigital("AA");

    try std.testing.expectEqual(@as(usize, 6), seq.len());

    const text = try seq.toText();
    defer allocator.free(text);
    try std.testing.expectEqualStrings("ACGTAA", text);
}

test "appendDigital: invalid text returns error" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    try std.testing.expectError(error.InvalidCharacter, seq.appendDigital("ACGZ"));
    // Original sequence should be unchanged.
    try std.testing.expectEqual(@as(usize, 4), seq.len());
}

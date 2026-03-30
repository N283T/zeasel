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
    pub fn reverseComplement(self: *Sequence) !void {
        try self.abc.reverseComplement(self.dsq);
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

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

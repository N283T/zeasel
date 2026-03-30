// Multiple Sequence Alignment (MSA) in digital encoding.

const std = @import("std");
const Allocator = std.mem.Allocator;
const alphabet_mod = @import("alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;
const Sequence = @import("sequence.zig").Sequence;

pub const Msa = struct {
    names: [][]const u8, // Sequence names [nseq], each is allocator-owned
    seqs: [][]u8, // Digital sequences [nseq][alen], each row is allocator-owned
    alen: usize, // Alignment length (columns)
    abc: *const Alphabet,
    allocator: Allocator,

    // Optional metadata (all allocator-owned if non-null)
    name: ?[]const u8 = null,
    accession: ?[]const u8 = null,
    description: ?[]const u8 = null,
    weights: ?[]f64 = null,

    pub fn nseq(self: Msa) usize {
        return self.names.len;
    }

    /// Create MSA from text sequences (all must be same length).
    /// This is the primary constructor — the Stockholm parser will use this.
    pub fn fromText(
        allocator: Allocator,
        abc: *const Alphabet,
        names: []const []const u8,
        text_seqs: []const []const u8,
    ) !Msa {
        if (names.len != text_seqs.len) return error.InvalidInput;
        if (names.len == 0) return error.InvalidInput;

        const alen = text_seqs[0].len;
        for (text_seqs) |s| {
            if (s.len != alen) return error.InvalidInput;
        }

        const n = names.len;

        var owned_names = try allocator.alloc([]const u8, n);
        var names_done: usize = 0;
        errdefer {
            for (0..names_done) |i| allocator.free(owned_names[i]);
            allocator.free(owned_names);
        }

        var owned_seqs = try allocator.alloc([]u8, n);
        var seqs_done: usize = 0;
        errdefer {
            for (0..seqs_done) |i| allocator.free(owned_seqs[i]);
            allocator.free(owned_seqs);
        }

        for (0..n) |i| {
            owned_names[i] = try allocator.dupe(u8, names[i]);
            names_done += 1;
            owned_seqs[i] = try abc.digitize(allocator, text_seqs[i]);
            seqs_done += 1;
        }

        return Msa{
            .names = owned_names,
            .seqs = owned_seqs,
            .alen = alen,
            .abc = abc,
            .allocator = allocator,
        };
    }

    /// Extract a single sequence from the alignment, removing gap columns.
    pub fn extractSeq(self: Msa, idx: usize) !Sequence {
        if (idx >= self.names.len) return error.OutOfBounds;

        // Count non-gap residues
        var count: usize = 0;
        for (self.seqs[idx]) |code| {
            if (!self.abc.isGap(code)) count += 1;
        }

        // Build ungapped sequence
        const dsq = try self.allocator.alloc(u8, count);
        errdefer self.allocator.free(dsq);
        var j: usize = 0;
        for (self.seqs[idx]) |code| {
            if (!self.abc.isGap(code)) {
                dsq[j] = code;
                j += 1;
            }
        }

        const name = try self.allocator.dupe(u8, self.names[idx]);
        errdefer self.allocator.free(name);

        return Sequence{
            .name = name,
            .accession = null,
            .description = null,
            .taxonomy_id = null,
            .dsq = dsq,
            .secondary_structure = null,
            .source = null,
            .abc = self.abc,
            .allocator = self.allocator,
        };
    }

    /// Select a subset of columns by boolean mask.
    /// mask.len must equal self.alen. Returns a new Msa with only true columns.
    pub fn selectColumns(self: Msa, mask: []const bool) !Msa {
        if (mask.len != self.alen) return error.InvalidInput;

        // Count selected columns
        var new_alen: usize = 0;
        for (mask) |m| {
            if (m) new_alen += 1;
        }

        const n = self.names.len;

        var new_names = try self.allocator.alloc([]const u8, n);
        var names_done: usize = 0;
        errdefer {
            for (0..names_done) |i| self.allocator.free(new_names[i]);
            self.allocator.free(new_names);
        }

        var new_seqs = try self.allocator.alloc([]u8, n);
        var seqs_done: usize = 0;
        errdefer {
            for (0..seqs_done) |i| self.allocator.free(new_seqs[i]);
            self.allocator.free(new_seqs);
        }

        for (0..n) |i| {
            new_names[i] = try self.allocator.dupe(u8, self.names[i]);
            names_done += 1;

            new_seqs[i] = try self.allocator.alloc(u8, new_alen);
            seqs_done += 1;

            var col: usize = 0;
            for (0..self.alen) |j| {
                if (mask[j]) {
                    new_seqs[i][col] = self.seqs[i][j];
                    col += 1;
                }
            }
        }

        return Msa{
            .names = new_names,
            .seqs = new_seqs,
            .alen = new_alen,
            .abc = self.abc,
            .allocator = self.allocator,
        };
    }

    /// Free all owned memory.
    pub fn deinit(self: *Msa) void {
        for (self.names) |n| self.allocator.free(n);
        self.allocator.free(self.names);
        for (self.seqs) |s| self.allocator.free(s);
        self.allocator.free(self.seqs);
        if (self.name) |n| self.allocator.free(n);
        if (self.accession) |a| self.allocator.free(a);
        if (self.description) |d| self.allocator.free(d);
        if (self.weights) |w| self.allocator.free(w);
    }
};

// --- Tests ---

test "fromText: basic construction" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2", "seq3" };
    const seqs = [_][]const u8{ "ACGT-", "A-GT-", "AC-T-" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 3), msa.nseq());
    try std.testing.expectEqual(@as(usize, 5), msa.alen);
}

test "fromText: mismatched sequence lengths return error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGT", "ACG" }; // different lengths

    try std.testing.expectError(
        error.InvalidInput,
        Msa.fromText(allocator, abc, &names, &seqs),
    );
}

test "fromText: empty input returns error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{};
    const seqs = [_][]const u8{};

    try std.testing.expectError(
        error.InvalidInput,
        Msa.fromText(allocator, abc, &names, &seqs),
    );
}

test "extractSeq: gaps are removed" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"seq1"};
    const seqs = [_][]const u8{"AC-GT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var seq = try msa.extractSeq(0);
    defer seq.deinit();

    const text = try seq.toText();
    defer allocator.free(text);

    try std.testing.expectEqualStrings("ACGT", text);
    try std.testing.expectEqualStrings("seq1", seq.name);
}

test "extractSeq: out of bounds returns error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"seq1"};
    const seqs = [_][]const u8{"ACGT-"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expectError(error.OutOfBounds, msa.extractSeq(1));
}

test "selectColumns: selects correct columns" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGTA", "TGCAT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    // Select columns 0, 2, 4 (true, false, true, false, true)
    const mask = [_]bool{ true, false, true, false, true };
    var sub = try msa.selectColumns(&mask);
    defer sub.deinit();

    try std.testing.expectEqual(@as(usize, 2), sub.nseq());
    try std.testing.expectEqual(@as(usize, 3), sub.alen);

    const text0 = try sub.abc.textize(allocator, sub.seqs[0]);
    defer allocator.free(text0);
    try std.testing.expectEqualStrings("AGA", text0);

    const text1 = try sub.abc.textize(allocator, sub.seqs[1]);
    defer allocator.free(text1);
    try std.testing.expectEqualStrings("TCT", text1);
}

test "selectColumns: mask length mismatch returns error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"seq1"};
    const seqs = [_][]const u8{"ACGTA"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    const mask = [_]bool{ true, false, true }; // wrong length
    try std.testing.expectError(error.InvalidInput, msa.selectColumns(&mask));
}

test "deinit: no memory leaks" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2", "seq3" };
    const seqs = [_][]const u8{ "ACGT-", "A-GT-", "AC-T-" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    msa.deinit();
    // std.testing.allocator detects leaks automatically
}

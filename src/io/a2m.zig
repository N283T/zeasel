// A2M (UCSC SAM) MSA format parser/writer.
//
// A2M is similar to aligned FASTA but uses case to distinguish
// consensus (uppercase) from insert (lowercase) columns.
// Periods and whitespace are ignored. Only '-' is a gap.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("../msa.zig").Msa;
const Alphabet = @import("../alphabet.zig").Alphabet;

/// Parse an A2M format alignment.
/// Insert columns (lowercase) are converted to gap in the consensus view.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    var name_list: std.ArrayList([]const u8) = .empty;
    defer name_list.deinit(allocator);
    var seq_list: std.ArrayList(std.ArrayList(u8)) = .empty;
    defer {
        for (seq_list.items) |*s| s.deinit(allocator);
        seq_list.deinit(allocator);
    }

    var lines = std.mem.splitScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        if (trimmed[0] == '>') {
            // Header line
            const name_start: usize = 1;
            var name_end = name_start;
            while (name_end < trimmed.len and trimmed[name_end] != ' ' and trimmed[name_end] != '\t') : (name_end += 1) {}
            const name = try allocator.dupe(u8, trimmed[name_start..name_end]);
            errdefer allocator.free(name);
            try name_list.append(allocator, name);

            var new_seq: std.ArrayList(u8) = .empty;
            try seq_list.append(allocator, new_seq);
            _ = &new_seq;
        } else if (seq_list.items.len > 0) {
            // Sequence line — append to current sequence
            const idx = seq_list.items.len - 1;
            for (trimmed) |ch| {
                // Skip whitespace and periods (A2M ignores them)
                if (ch == ' ' or ch == '\t' or ch == '.' or ch == 'O' or ch == 'o') continue;
                // Lowercase = insert column, convert to uppercase for digitization
                const out = if (ch >= 'a' and ch <= 'z') ch - 32 else ch;
                try seq_list.items[idx].append(allocator, out);
            }
        }
    }

    if (name_list.items.len == 0) return error.InvalidFormat;

    // Pad sequences to same length
    var max_len: usize = 0;
    for (seq_list.items) |s| {
        if (s.items.len > max_len) max_len = s.items.len;
    }
    for (seq_list.items) |*s| {
        while (s.items.len < max_len) {
            try s.append(allocator, '-');
        }
    }

    // Convert to text slices
    const n = name_list.items.len;
    var text_seqs = try allocator.alloc([]const u8, n);
    defer allocator.free(text_seqs);
    for (0..n) |i| {
        text_seqs[i] = seq_list.items[i].items;
    }

    const result = try Msa.fromText(allocator, abc, name_list.items, text_seqs);

    // Free names (fromText duped them)
    for (name_list.items) |name| allocator.free(name);

    return result;
}

/// Write an MSA in A2M format (uppercase for all columns, no insert distinction).
pub fn write(msa: Msa, dest: std.io.AnyWriter) !void {
    for (0..msa.nseq()) |i| {
        try dest.writeByte('>');
        try dest.writeAll(msa.names[i]);
        try dest.writeByte('\n');
        for (msa.seqs[i]) |code| {
            try dest.writeByte(msa.abc.decode(code));
        }
        try dest.writeByte('\n');
    }
}

// --- Tests ---

test "parse: simple A2M" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\>seq1
        \\ACGT
        \\>seq2
        \\TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 4), msa.alen);
    try std.testing.expectEqualStrings("seq1", msa.names[0]);
}

test "parse: ignores periods and whitespace" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\>seq1
        \\A.C G.T
        \\>seq2
        \\T.G C.A
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 4), msa.alen);
}

test "write: round trip" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };
    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try write(msa, buf.writer(allocator).any());

    var msa2 = try parse(allocator, abc, buf.items);
    defer msa2.deinit();

    try std.testing.expect(msa.compare(msa2));
}

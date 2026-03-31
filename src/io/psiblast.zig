// PSI-BLAST MSA format parser/writer.
//
// PSI-BLAST format is a simple interleaved alignment where:
// - First sequence is the query (no name on first block's first line in some variants)
// - Each line is: name  sequence_chunk
// - Blocks are separated by blank lines
// - Only '-' allowed for gaps; 'O' residues disallowed

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("../msa.zig").Msa;
const Alphabet = @import("../alphabet.zig").Alphabet;

/// Parse a PSI-BLAST format alignment.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    var name_list: std.ArrayList([]const u8) = .empty;
    defer {
        for (name_list.items) |n| allocator.free(n);
        name_list.deinit(allocator);
    }
    var seq_bufs: std.ArrayList(std.ArrayList(u8)) = .empty;
    defer {
        for (seq_bufs.items) |*s| s.deinit(allocator);
        seq_bufs.deinit(allocator);
    }

    var name_to_idx = std.StringHashMap(usize).init(allocator);
    defer name_to_idx.deinit();

    var lines = std.mem.splitScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        // Split into name and sequence
        var it = std.mem.tokenizeAny(u8, trimmed, " \t");
        const name = it.next() orelse continue;
        const seq_part = it.rest();
        if (seq_part.len == 0) continue;

        if (name_to_idx.get(name)) |idx| {
            // Existing sequence — append
            for (seq_part) |ch| {
                if (ch != ' ' and ch != '\t') {
                    try seq_bufs.items[idx].append(allocator, ch);
                }
            }
        } else {
            // New sequence
            const name_copy = try allocator.dupe(u8, name);
            errdefer allocator.free(name_copy);
            const idx = name_list.items.len;
            try name_to_idx.put(name_copy, idx);
            try name_list.append(allocator, name_copy);

            var buf: std.ArrayList(u8) = .empty;
            for (seq_part) |ch| {
                if (ch != ' ' and ch != '\t') {
                    try buf.append(allocator, ch);
                }
            }
            try seq_bufs.append(allocator, buf);
        }
    }

    if (name_list.items.len == 0) return error.InvalidFormat;

    const n = name_list.items.len;
    var text_seqs = try allocator.alloc([]const u8, n);
    defer allocator.free(text_seqs);
    for (0..n) |i| {
        text_seqs[i] = seq_bufs.items[i].items;
    }

    const result = try Msa.fromText(allocator, abc, name_list.items, text_seqs);

    // fromText duped names, free our copies
    for (name_list.items) |nm| allocator.free(nm);
    name_list.items.len = 0; // prevent double-free in defer

    return result;
}

/// Write an MSA in PSI-BLAST format.
pub fn write(msa: Msa, dest: std.io.AnyWriter) !void {
    const block_size: usize = 60;
    var col: usize = 0;

    while (col < msa.alen) {
        const end = @min(col + block_size, msa.alen);
        for (0..msa.nseq()) |i| {
            try dest.writeAll(msa.names[i]);
            try dest.writeAll("  ");
            for (col..end) |j| {
                try dest.writeByte(msa.abc.decode(msa.seqs[i][j]));
            }
            try dest.writeByte('\n');
        }
        try dest.writeByte('\n');
        col = end;
    }
}

// --- Tests ---

test "parse: simple PSI-BLAST" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").amino;

    const data =
        \\query  ACDEF
        \\hit1   ACDEG
        \\
        \\query  GHIKL
        \\hit1   GHIKM
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 10), msa.alen);
    try std.testing.expectEqualStrings("query", msa.names[0]);
}

test "write: round trip" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").amino;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACDEF", "GHIKL" };
    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try write(msa, buf.writer(allocator).any());

    var msa2 = try parse(allocator, abc, buf.items);
    defer msa2.deinit();

    try std.testing.expect(msa.compare(msa2));
}

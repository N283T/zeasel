// SELEX MSA format parser/writer.
//
// SELEX is a block-based format where each line is: name  sequence
// Spaces in sequence data are treated as gaps (mapped to '-').
// Annotation lines: #=RF, #=CS, #=SS, #=SA can appear in blocks.
// Blocks are separated by blank lines.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("../msa.zig").Msa;
const Alphabet = @import("../alphabet.zig").Alphabet;

/// Parse a SELEX format alignment.
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
        // Skip comment/annotation lines
        if (trimmed[0] == '#') continue;

        // Find where name ends — first whitespace
        var name_end: usize = 0;
        while (name_end < trimmed.len and trimmed[name_end] != ' ' and trimmed[name_end] != '\t') : (name_end += 1) {}
        if (name_end == 0) continue;
        const name = trimmed[0..name_end];

        // Rest is sequence (spaces become gaps)
        const seq_part = if (name_end < trimmed.len) trimmed[name_end..] else "";

        // Extract sequence: skip leading whitespace, then map spaces to gaps
        const seq_trimmed = std.mem.trimLeft(u8, seq_part, " \t");

        if (name_to_idx.get(name)) |idx| {
            for (seq_trimmed) |ch| {
                if (ch == '\t') continue;
                const out: u8 = if (ch == ' ') '-' else ch;
                try seq_bufs.items[idx].append(allocator, out);
            }
        } else {
            const name_copy = try allocator.dupe(u8, name);
            errdefer allocator.free(name_copy);
            const idx = name_list.items.len;
            try name_to_idx.put(name_copy, idx);
            try name_list.append(allocator, name_copy);

            var buf: std.ArrayList(u8) = .empty;
            for (seq_trimmed) |ch| {
                if (ch == '\t') continue;
                const out: u8 = if (ch == ' ') '-' else ch;
                try buf.append(allocator, out);
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

    for (name_list.items) |nm| allocator.free(nm);
    name_list.items.len = 0;

    return result;
}

/// Write an MSA in SELEX format (single block).
pub fn write(msa: Msa, dest: std.io.AnyWriter) !void {
    // Find max name length for alignment
    var max_name: usize = 0;
    for (msa.names) |name| {
        if (name.len > max_name) max_name = name.len;
    }

    for (0..msa.nseq()) |i| {
        try dest.writeAll(msa.names[i]);
        for (0..max_name - msa.names[i].len + 2) |_| try dest.writeByte(' ');
        for (msa.seqs[i]) |code| {
            try dest.writeByte(msa.abc.decode(code));
        }
        try dest.writeByte('\n');
    }
}

// --- Tests ---

test "parse: simple SELEX" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\seq1  ACGT
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 4), msa.alen);
}

test "parse: spaces as gaps" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\seq1  AC-GT
        \\seq2  TGCAA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 5), msa.alen);
}

test "parse: skips comment lines" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\# comment
        \\seq1  ACGT
        \\#=RF  xxxx
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
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

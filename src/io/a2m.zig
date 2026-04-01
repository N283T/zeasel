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
///
/// A2M uses case to distinguish consensus columns (uppercase/'-') from
/// insert columns (lowercase/'.'). The parsed MSA contains only consensus
/// columns. Insert residues (lowercase) are counted per-sequence to
/// determine the maximum insert length before each consensus column, then
/// gaps are inserted for sequences with fewer inserts at that position.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    var name_list: std.ArrayList([]const u8) = .empty;
    defer {
        for (name_list.items) |n| allocator.free(n);
        name_list.deinit(allocator);
    }

    // First pass: collect raw character sequences (preserving case).
    var raw_seqs: std.ArrayList(std.ArrayList(u8)) = .empty;
    defer {
        for (raw_seqs.items) |*s| s.deinit(allocator);
        raw_seqs.deinit(allocator);
    }

    var lines = std.mem.splitScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        if (trimmed[0] == '>') {
            const name_start: usize = 1;
            var name_end = name_start;
            while (name_end < trimmed.len and trimmed[name_end] != ' ' and trimmed[name_end] != '\t') : (name_end += 1) {}
            const name = try allocator.dupe(u8, trimmed[name_start..name_end]);
            errdefer allocator.free(name);
            try name_list.append(allocator, name);

            var new_seq: std.ArrayList(u8) = .empty;
            try raw_seqs.append(allocator, new_seq);
            _ = &new_seq;
        } else if (raw_seqs.items.len > 0) {
            const idx = raw_seqs.items.len - 1;
            for (trimmed) |ch| {
                if (ch == ' ' or ch == '\t') continue;
                try raw_seqs.items[idx].append(allocator, ch);
            }
        }
    }

    if (name_list.items.len == 0) return error.InvalidFormat;

    // Determine the number of consensus columns and max inserts before each.
    // Consensus characters: uppercase letters, '-', '*'
    // Insert characters: lowercase letters, '.'
    const nseq = name_list.items.len;

    // Count consensus columns from the first sequence to determine ncols.
    var ncols: usize = 0;
    for (raw_seqs.items[0].items) |ch| {
        if (isConsensusChar(ch)) ncols += 1;
    }

    // Compute max insert length before each consensus column (ncols + 1 positions:
    // before col 0, between col 0 and 1, ..., after last col).
    var max_inserts = try allocator.alloc(usize, ncols + 1);
    defer allocator.free(max_inserts);
    @memset(max_inserts, 0);

    for (raw_seqs.items) |raw| {
        var col_idx: usize = 0;
        var ins_count: usize = 0;
        for (raw.items) |ch| {
            if (isConsensusChar(ch)) {
                if (ins_count > max_inserts[col_idx]) max_inserts[col_idx] = ins_count;
                ins_count = 0;
                col_idx += 1;
            } else if (isInsertChar(ch)) {
                ins_count += 1;
            }
        }
        // Trailing inserts after last consensus column
        if (ins_count > max_inserts[col_idx]) max_inserts[col_idx] = ins_count;
    }

    // Compute total alignment length
    var alen: usize = ncols;
    for (max_inserts) |mi| alen += mi;

    // Second pass: build padded alignment.
    var seq_list = try allocator.alloc(std.ArrayList(u8), nseq);
    defer {
        for (seq_list) |*s| s.deinit(allocator);
        allocator.free(seq_list);
    }
    for (seq_list) |*s| s.* = .empty;

    for (0..nseq) |si| {
        const raw = raw_seqs.items[si].items;
        var col_idx: usize = 0;
        var ins_count: usize = 0;
        for (raw) |ch| {
            if (isConsensusChar(ch)) {
                // Pad remaining insert slots with gaps
                while (ins_count < max_inserts[col_idx]) : (ins_count += 1) {
                    try seq_list[si].append(allocator, '-');
                }
                // Write the consensus character (uppercase)
                const out = if (ch >= 'a' and ch <= 'z') ch - 32 else ch;
                try seq_list[si].append(allocator, out);
                col_idx += 1;
                ins_count = 0;
            } else if (isInsertChar(ch)) {
                // Write insert character as uppercase
                const out = if (ch >= 'a' and ch <= 'z') ch - 32 else ch;
                // '.' becomes '-'
                const final = if (out == '.') @as(u8, '-') else out;
                try seq_list[si].append(allocator, final);
                ins_count += 1;
            }
        }
        // Pad trailing inserts
        while (ins_count < max_inserts[col_idx]) : (ins_count += 1) {
            try seq_list[si].append(allocator, '-');
        }
    }

    // Convert to text slices
    var text_seqs = try allocator.alloc([]const u8, nseq);
    defer allocator.free(text_seqs);
    for (0..nseq) |i| {
        text_seqs[i] = seq_list[i].items;
    }

    const result = try Msa.fromText(allocator, abc, name_list.items, text_seqs);

    // Free names (fromText duped them). Clear list to prevent double-free in defer.
    for (name_list.items) |name| allocator.free(name);
    name_list.items.len = 0;

    return result;
}

fn isConsensusChar(ch: u8) bool {
    // 'O' is a free-insertion module marker in A2M, not pyrrolysine — skip it.
    if (ch == 'O') return false;
    return (ch >= 'A' and ch <= 'Z') or ch == '-' or ch == '*';
}

fn isInsertChar(ch: u8) bool {
    // 'o' is also skipped (lowercase of FIM marker 'O').
    if (ch == 'o') return false;
    return (ch >= 'a' and ch <= 'z') or ch == '.';
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

test "parse: insert columns handled correctly" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    // A2M format: uppercase/'-' = consensus, lowercase/'.' = insert
    // seq1: consensus A, insert "ac", consensus G, consensus T
    // seq2: consensus C, insert "..",  consensus G, consensus A
    // Both have 3 consensus columns. Max insert before col 1 = 2.
    // Final alignment: 3 consensus + 2 insert = 5 columns
    const data =
        \\>seq1
        \\AacGT
        \\>seq2
        \\C..GA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    // 3 consensus cols + 2 max inserts = 5
    try std.testing.expectEqual(@as(usize, 5), msa.alen);
}

test "parse: periods become gaps in insert columns" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    // '.' is an insert gap character in A2M
    const data =
        \\>seq1
        \\AaC
        \\>seq2
        \\A.C
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 3), msa.alen); // 2 consensus + 1 insert
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

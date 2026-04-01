// Clustal format parser and writer for multiple sequence alignments.
//
// Format:
//   - First line starts with "CLUSTAL" (header)
//   - Blank lines separate blocks
//   - Sequence lines: <name> <sequence_fragment>
//   - Conservation lines (starting with spaces) are skipped
//   - Interleaved: each sequence name appears once per block

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const Msa = @import("../msa.zig").Msa;

/// Parse a Clustal format alignment from a byte buffer.
/// Returns an Msa owning all allocated memory.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    var lines = std.mem.splitScalar(u8, data, '\n');

    // First non-blank line must start with "CLUSTAL" and contain "alignment" (case-insensitive).
    var found_header = false;
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, "\r");
        if (trimmed.len == 0) continue;
        if (!std.mem.startsWith(u8, trimmed, "CLUSTAL")) return error.InvalidFormat;
        // Validate that the rest of the header line contains "alignment" (case-insensitive).
        const rest = trimmed["CLUSTAL".len..];
        if (!containsCaseInsensitive(rest, "alignment")) return error.InvalidFormat;
        found_header = true;
        break;
    }
    if (!found_header) return error.InvalidFormat;

    // Ordered list of sequence names (first-seen order).
    var name_list = std.ArrayList([]const u8){};
    defer {
        for (name_list.items) |n| allocator.free(n);
        name_list.deinit(allocator);
    }

    // name_index maps name to its position in name_list / seq_parts.
    var name_index = std.StringHashMap(usize).init(allocator);
    defer name_index.deinit();

    // seq_parts[i] collects all text fragments for sequence i across blocks.
    var seq_parts = std.ArrayList(std.ArrayList(u8)){};
    defer {
        for (seq_parts.items) |*parts| parts.deinit(allocator);
        seq_parts.deinit(allocator);
    }

    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, "\r");

        // Skip blank lines (block separators).
        if (trimmed.len == 0) continue;

        // Skip conservation lines (start with space).
        if (trimmed[0] == ' ' or trimmed[0] == '\t') continue;

        // Parse: <name> <whitespace> <sequence_fragment>
        // Find the end of the name token.
        var name_end: usize = 0;
        while (name_end < trimmed.len and trimmed[name_end] != ' ' and trimmed[name_end] != '\t') {
            name_end += 1;
        }
        if (name_end == 0) continue;

        const seq_name = trimmed[0..name_end];

        // Skip whitespace between name and sequence.
        var seq_start = name_end;
        while (seq_start < trimmed.len and (trimmed[seq_start] == ' ' or trimmed[seq_start] == '\t')) {
            seq_start += 1;
        }

        // The sequence fragment may be followed by an optional column count;
        // we stop at the first whitespace after the sequence characters.
        var seq_end = seq_start;
        while (seq_end < trimmed.len and trimmed[seq_end] != ' ' and trimmed[seq_end] != '\t') {
            seq_end += 1;
        }

        const seq_frag = trimmed[seq_start..seq_end];

        // Look up or register the sequence name.
        if (name_index.get(seq_name)) |idx| {
            try seq_parts.items[idx].appendSlice(allocator, seq_frag);
        } else {
            const name_copy = try allocator.dupe(u8, seq_name);
            errdefer allocator.free(name_copy);

            const idx = name_list.items.len;
            try name_list.append(allocator, name_copy);
            errdefer _ = name_list.pop();

            // Use the name slice that is owned by name_list.
            try name_index.put(name_list.items[idx], idx);

            var parts = std.ArrayList(u8){};
            try parts.appendSlice(allocator, seq_frag);
            try seq_parts.append(allocator, parts);
        }
    }

    if (name_list.items.len == 0) return error.InvalidFormat;

    // Validate that all sequences have the same length.
    const alen = seq_parts.items[0].items.len;
    for (seq_parts.items) |parts| {
        if (parts.items.len != alen) return error.InvalidInput;
    }

    // Build text_seqs slice for Msa.fromText.
    const n = name_list.items.len;
    const text_seqs = try allocator.alloc([]const u8, n);
    defer allocator.free(text_seqs);
    for (0..n) |i| {
        text_seqs[i] = seq_parts.items[i].items;
    }

    return Msa.fromText(allocator, abc, name_list.items, text_seqs);
}

/// Write an MSA in Clustal format.
/// Uses a fixed block width of 60 residues per line.
pub fn write(dest: std.io.AnyWriter, m: Msa) !void {
    const block_width: usize = 60;

    try dest.writeAll("CLUSTAL W (1.83) multiple sequence alignment\n");

    // Find longest name for padding.
    var max_name_len: usize = 0;
    for (m.names) |name| {
        if (name.len > max_name_len) max_name_len = name.len;
    }

    const col_sep: usize = 6; // spaces between name and sequence
    var col: usize = 0;
    while (col < m.alen) {
        const end = @min(col + block_width, m.alen);

        try dest.writeByte('\n');

        for (0..m.nseq()) |i| {
            // Write name padded to max_name_len + col_sep.
            try dest.writeAll(m.names[i]);
            const pad = max_name_len + col_sep - m.names[i].len;
            for (0..pad) |_| try dest.writeByte(' ');

            // Write sequence fragment for this block.
            const frag = m.seqs[i][col..end];
            const text = try m.abc.textize(m.allocator, frag);
            defer m.allocator.free(text);
            try dest.writeAll(text);
            try dest.writeByte('\n');
        }

        // Write a blank conservation line (all spaces) to match real Clustal output.
        for (0..max_name_len + col_sep) |_| try dest.writeByte(' ');
        try dest.writeByte('\n');

        col = end;
    }
}

/// Check if haystack contains needle, case-insensitive.
fn containsCaseInsensitive(haystack: []const u8, needle: []const u8) bool {
    if (needle.len == 0) return true;
    if (haystack.len < needle.len) return false;
    const limit = haystack.len - needle.len + 1;
    for (0..limit) |i| {
        var match = true;
        for (0..needle.len) |j| {
            const h = std.ascii.toLower(haystack[i + j]);
            const n = std.ascii.toLower(needle[j]);
            if (h != n) {
                match = false;
                break;
            }
        }
        if (match) return true;
    }
    return false;
}

// --- Tests ---

const alphabet_mod = @import("../alphabet.zig");

test "parse: simple two-sequence alignment" {
    const allocator = std.testing.allocator;
    const data =
        \\CLUSTAL W (1.83) multiple sequence alignment
        \\
        \\seq1      ACGT
        \\seq2      ACGT
        \\
        \\
    ;

    var msa = try parse(allocator, &alphabet_mod.amino, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 4), msa.alen);
    try std.testing.expectEqualStrings("seq1", msa.names[0]);
    try std.testing.expectEqualStrings("seq2", msa.names[1]);
}

test "parse: interleaved two blocks" {
    const allocator = std.testing.allocator;
    const data =
        \\CLUSTAL W (1.83) multiple sequence alignment
        \\
        \\seq1      ACDE
        \\seq2      ACDE
        \\          ****
        \\
        \\seq1      FGHI
        \\seq2      FGHI
        \\          ****
        \\
    ;

    var msa = try parse(allocator, &alphabet_mod.amino, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 8), msa.alen);

    const text0 = try msa.abc.textize(allocator, msa.seqs[0]);
    defer allocator.free(text0);
    try std.testing.expectEqualStrings("ACDEFGHI", text0);
}

test "parse: conservation line is skipped" {
    const allocator = std.testing.allocator;
    const data =
        \\CLUSTAL W (1.83) multiple sequence alignment
        \\
        \\seq1      ACDE
        \\seq2      ACDE
        \\          ****
        \\
    ;

    var msa = try parse(allocator, &alphabet_mod.amino, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
}

test "parse: missing CLUSTAL header returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(
        error.InvalidFormat,
        parse(allocator, &alphabet_mod.amino, "seq1  ACDE\n"),
    );
}

test "parse: CLUSTAL header without 'alignment' keyword returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(
        error.InvalidFormat,
        parse(allocator, &alphabet_mod.amino, "CLUSTAL W (1.83)\n\nseq1  ACDE\n"),
    );
}

test "parse: CLUSTAL header with uppercase ALIGNMENT accepted" {
    const allocator = std.testing.allocator;
    const data =
        \\CLUSTAL W (1.83) multiple sequence ALIGNMENT
        \\
        \\seq1      ACDE
        \\seq2      ACDE
        \\
    ;

    var msa = try parse(allocator, &alphabet_mod.amino, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
}

test "parse: CLUSTAL O header with 'alignment' accepted" {
    const allocator = std.testing.allocator;
    const data =
        \\CLUSTAL O(1.2.4) multiple sequence alignment
        \\
        \\seq1      ACDE
        \\seq2      ACDE
        \\
    ;

    var msa = try parse(allocator, &alphabet_mod.amino, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
}

test "write: basic output" {
    const allocator = std.testing.allocator;
    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACDEFGHIKLMNPQRS", "ACDEFGHIKLMNPQRS" };

    var msa = try Msa.fromText(allocator, &alphabet_mod.amino, &names, &seqs);
    defer msa.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa);

    // Must start with CLUSTAL header.
    try std.testing.expect(std.mem.startsWith(u8, buf.items, "CLUSTAL"));
    // Must contain the sequence name.
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "seq1") != null);
}

test "round trip" {
    const allocator = std.testing.allocator;
    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACDE-FGHI", "ACD--FGHI" };

    var msa = try Msa.fromText(allocator, &alphabet_mod.amino, &names, &seqs);
    defer msa.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa);

    var msa2 = try parse(allocator, &alphabet_mod.amino, buf.items);
    defer msa2.deinit();

    try std.testing.expectEqual(msa.nseq(), msa2.nseq());
    try std.testing.expectEqual(msa.alen, msa2.alen);
    for (0..msa.nseq()) |i| {
        try std.testing.expectEqualStrings(msa.names[i], msa2.names[i]);
        try std.testing.expectEqualSlices(u8, msa.seqs[i], msa2.seqs[i]);
    }
}

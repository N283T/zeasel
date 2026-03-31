// PHYLIP interleaved and sequential MSA format parser/writer.
//
// PHYLIP format starts with a header line: "  nseq  alen"
// Then nseq lines of "name       sequence..." (name is fixed-width, typically 10 chars).
// Interleaved format repeats blocks of nseq lines until alen is reached.
// Sequential format has all residues for each sequence contiguously.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Msa = @import("../msa.zig").Msa;
const Alphabet = @import("../alphabet.zig").Alphabet;

/// Parse a PHYLIP format alignment (interleaved or sequential).
/// Detects variant automatically. Returns a new Msa.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    var lines = std.mem.splitScalar(u8, data, '\n');

    // Parse header: nseq  alen
    const header = skipEmpty(&lines) orelse return error.InvalidFormat;
    var header_it = std.mem.tokenizeAny(u8, header, " \t");
    const nseq_str = header_it.next() orelse return error.InvalidFormat;
    const alen_str = header_it.next() orelse return error.InvalidFormat;
    const nseq = std.fmt.parseInt(usize, nseq_str, 10) catch return error.InvalidFormat;
    const alen = std.fmt.parseInt(usize, alen_str, 10) catch return error.InvalidFormat;

    if (nseq == 0 or alen == 0) return error.InvalidFormat;

    // Collect remaining non-empty lines
    var all_lines: std.ArrayList([]const u8) = .empty;
    defer all_lines.deinit(allocator);
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, " \t\r");
        if (trimmed.len > 0) {
            try all_lines.append(allocator, trimmed);
        }
    }

    // Determine name width from first line
    // PHYLIP standard: 10 chars for name, but we detect it
    const name_width = detectNameWidth(all_lines.items[0], alen);

    // Parse names and sequences
    var names = try allocator.alloc([]const u8, nseq);
    var names_done: usize = 0;
    errdefer {
        for (0..names_done) |i| allocator.free(names[i]);
        allocator.free(names);
    }

    var seq_bufs = try allocator.alloc(std.ArrayList(u8), nseq);
    var bufs_done: usize = 0;
    errdefer {
        for (0..bufs_done) |i| seq_bufs[i].deinit(allocator);
        allocator.free(seq_bufs);
    }
    for (0..nseq) |i| {
        seq_bufs[i] = .empty;
        bufs_done = i + 1;
    }

    // First block: has names
    for (0..nseq) |i| {
        if (i >= all_lines.items.len) return error.InvalidFormat;
        const line = all_lines.items[i];
        const raw_name = std.mem.trimRight(u8, line[0..@min(name_width, line.len)], " \t");
        names[i] = try allocator.dupe(u8, raw_name);
        names_done = i + 1;

        const seq_part = if (line.len > name_width) line[name_width..] else "";
        for (seq_part) |ch| {
            if (ch != ' ' and ch != '\t') {
                try seq_bufs[i].append(allocator, ch);
            }
        }
    }

    // Detect interleaved vs sequential format.
    // After the first block of nseq named lines, check the next line (if present).
    // In interleaved format, subsequent blocks have NO names — just raw sequence.
    // In sequential format, the next line is a continuation of the first sequence,
    // OR (once first seq is complete) a new named line for the next sequence.
    // Heuristic: count total non-space residues that would be contributed by
    // interleaved parsing. If (total_lines - nseq) is divisible by nseq AND
    // that gives the right residue count, use interleaved. Otherwise, sequential.
    const is_sequential = blk: {
        if (seq_bufs[0].items.len >= alen) break :blk false;
        if (all_lines.items.len <= nseq) break :blk false;
        // Count residues in subsequent lines (no name field).
        const remaining_lines = all_lines.items.len - nseq;
        if (remaining_lines % nseq != 0) break :blk true;
        // In interleaved, each subsequent block contributes equally to each sequence.
        // Count residues from interleaved blocks for seq 0.
        var interleaved_residues = seq_bufs[0].items.len;
        var li: usize = nseq;
        while (li < all_lines.items.len) : (li += nseq) {
            const line = all_lines.items[li];
            for (line) |ch| {
                if (ch != ' ' and ch != '\t') interleaved_residues += 1;
            }
        }
        if (interleaved_residues == alen) break :blk false;
        break :blk true;
    };

    if (is_sequential) {
        // Sequential format: for each sequence, name+data on first line then
        // continuation lines until alen residues are reached.
        // Re-extract names and sequence data from scratch since the initial
        // pass assumed interleaved layout (first nseq lines as named lines).
        for (0..names_done) |i| allocator.free(names[i]);
        names_done = 0;
        for (0..nseq) |i| {
            seq_bufs[i].clearRetainingCapacity();
        }

        var seq_line_idx: usize = 0;
        for (0..nseq) |i| {
            if (seq_line_idx >= all_lines.items.len) return error.InvalidFormat;
            const first_line = all_lines.items[seq_line_idx];
            seq_line_idx += 1;

            const nw = detectNameWidth(first_line, alen);
            const raw_name = std.mem.trimRight(u8, first_line[0..@min(nw, first_line.len)], " \t");
            names[i] = try allocator.dupe(u8, raw_name);
            names_done = i + 1;

            const seq_part = if (first_line.len > nw) first_line[nw..] else "";
            for (seq_part) |ch| {
                if (ch != ' ' and ch != '\t') {
                    try seq_bufs[i].append(allocator, ch);
                }
            }

            while (seq_bufs[i].items.len < alen) {
                if (seq_line_idx >= all_lines.items.len) break;
                const cont_line = all_lines.items[seq_line_idx];
                seq_line_idx += 1;
                for (cont_line) |ch| {
                    if (ch != ' ' and ch != '\t') {
                        try seq_bufs[i].append(allocator, ch);
                    }
                }
            }
        }
    } else {
        // Interleaved: subsequent blocks have no names.
        var line_idx: usize = nseq;
        while (seq_bufs[0].items.len < alen) {
            for (0..nseq) |i| {
                if (line_idx >= all_lines.items.len) break;
                const line = all_lines.items[line_idx];
                line_idx += 1;
                for (line) |ch| {
                    if (ch != ' ' and ch != '\t') {
                        try seq_bufs[i].append(allocator, ch);
                    }
                }
            }
            if (line_idx >= all_lines.items.len) break;
        }
    }

    // Convert to text sequences
    var text_seqs = try allocator.alloc([]const u8, nseq);
    var text_done: usize = 0;
    errdefer {
        for (0..text_done) |i| allocator.free(text_seqs[i]);
        allocator.free(text_seqs);
    }
    for (0..nseq) |i| {
        text_seqs[i] = try seq_bufs[i].toOwnedSlice(allocator);
        text_done = i + 1;
    }

    // Free seq_bufs (already transferred to text_seqs)
    for (0..nseq) |i| seq_bufs[i].deinit(allocator);
    allocator.free(seq_bufs);

    const result = try Msa.fromText(allocator, abc, names, text_seqs);

    // Free temporary text sequences and names
    for (0..nseq) |i| {
        allocator.free(text_seqs[i]);
        allocator.free(names[i]);
    }
    allocator.free(text_seqs);
    allocator.free(names);

    return result;
}

/// Write an MSA in PHYLIP interleaved format.
pub fn write(msa: Msa, dest: std.io.AnyWriter) !void {
    // Header
    try dest.print(" {d} {d}\n", .{ msa.nseq(), msa.alen });

    const block_size: usize = 60;
    var col: usize = 0;

    while (col < msa.alen) {
        const end = @min(col + block_size, msa.alen);
        for (0..msa.nseq()) |i| {
            if (col == 0) {
                // First block: print name (padded to 10 chars)
                const name = msa.names[i];
                try dest.writeAll(name);
                if (name.len < 10) {
                    for (0..10 - name.len) |_| try dest.writeByte(' ');
                }
            }
            // Print sequence chunk
            for (col..end) |j| {
                try dest.writeByte(msa.abc.decode(msa.seqs[i][j]));
            }
            try dest.writeByte('\n');
        }
        if (end < msa.alen) try dest.writeByte('\n');
        col = end;
    }
}

fn skipEmpty(lines: *std.mem.SplitIterator(u8, .scalar)) ?[]const u8 {
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, " \t\r");
        if (trimmed.len > 0) return trimmed;
    }
    return null;
}

fn detectNameWidth(first_line: []const u8, alen: usize) usize {
    _ = alen;
    // Standard PHYLIP uses 10-char names. Try to detect by finding
    // where sequence data starts (first non-space after initial token).
    // Default to 10 if we can't detect.
    if (first_line.len <= 10) return first_line.len;

    // Check if there's a clear space boundary around position 10
    if (first_line.len > 10 and first_line[9] != ' ' and first_line[10] != ' ') {
        // Strict format: names are exactly 10 chars
        return 10;
    }

    // Relaxed: find first space followed by sequence-like chars
    var i: usize = 1;
    while (i < first_line.len) : (i += 1) {
        if (first_line[i] == ' ' or first_line[i] == '\t') {
            // Skip spaces to find sequence start
            var j = i;
            while (j < first_line.len and (first_line[j] == ' ' or first_line[j] == '\t')) : (j += 1) {}
            if (j < first_line.len) return j;
        }
    }
    return 10;
}

// --- Tests ---

test "parse: simple interleaved PHYLIP" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\ 2 8
        \\seq1      ACGT
        \\seq2      TGCA
        \\
        \\ACGT
        \\TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 8), msa.alen);
    try std.testing.expectEqualStrings("seq1", msa.names[0]);
}

test "parse: sequential PHYLIP (single block)" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\ 2 4
        \\seq1      ACGT
        \\seq2      TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 4), msa.alen);
}

test "parse: sequential PHYLIP (multi-line)" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    // Sequential: each sequence's data is contiguous (name line + continuation lines).
    // 2 sequences, 8 residues each. Total non-blank lines = 4 (not divisible by 2 as blocks).
    const data =
        \\ 2 8
        \\seq1      ACGT
        \\ACGT
        \\seq2      TGCA
        \\TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 8), msa.alen);
    try std.testing.expectEqualStrings("seq1", msa.names[0]);
    try std.testing.expectEqualStrings("seq2", msa.names[1]);
}

test "write: round trip" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try write(msa, buf.writer(allocator).any());

    // Parse back
    var msa2 = try parse(allocator, abc, buf.items);
    defer msa2.deinit();

    try std.testing.expectEqual(msa.nseq(), msa2.nseq());
    try std.testing.expectEqual(msa.alen, msa2.alen);
}

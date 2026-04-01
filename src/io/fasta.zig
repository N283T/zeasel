// FASTA format parser and writer.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const Sequence = @import("../sequence.zig").Sequence;

pub const ParseError = error{
    InvalidFormat,
    InvalidCharacter,
    OutOfMemory,
};

/// Parse all FASTA records from a byte buffer.
/// Returns an owned slice of Sequences. Caller owns the slice and each Sequence.
pub fn parseAll(allocator: Allocator, abc: *const Alphabet, data: []const u8) ![]Sequence {
    var sequences = std.ArrayList(Sequence){};
    errdefer {
        for (sequences.items) |*seq| seq.deinit();
        sequences.deinit(allocator);
    }

    var pos: usize = 0;
    while (pos < data.len) {
        // Skip blank lines and whitespace before '>'
        while (pos < data.len) {
            const c = data[pos];
            if (c == '\n' or c == '\r' or c == ' ' or c == '\t') {
                pos += 1;
            } else {
                break;
            }
        }
        if (pos >= data.len) break;

        if (data[pos] != '>') return error.InvalidFormat;

        var seq = try parseOne(allocator, abc, data, &pos);
        errdefer seq.deinit();
        try sequences.append(allocator, seq);
    }

    return sequences.toOwnedSlice(allocator);
}

/// Parse a single FASTA record starting at data[pos] (which must be '>').
/// Updates pos to point past the parsed record (to the next '>' or EOF).
pub fn parseOne(allocator: Allocator, abc: *const Alphabet, data: []const u8, pos: *usize) !Sequence {
    // Consume '>'
    pos.* += 1;

    // --- Parse header line ---
    const header_start = pos.*;
    while (pos.* < data.len and data[pos.*] != '\n' and data[pos.*] != '\r') {
        pos.* += 1;
    }
    const header = data[header_start..pos.*];

    // Skip newline(s) after header
    if (pos.* < data.len and data[pos.*] == '\r') pos.* += 1;
    if (pos.* < data.len and data[pos.*] == '\n') pos.* += 1;

    // Split header into name, optional accession, and optional description.
    // Format: >name [accession [description...]]
    // The second whitespace-delimited word is treated as accession.
    var name: []const u8 = header;
    var accession: ?[]const u8 = null;
    var description: ?[]const u8 = null;

    for (header, 0..) |c, i| {
        if (c == ' ' or c == '\t') {
            name = header[0..i];
            // Skip whitespace after name.
            var rest_start = i + 1;
            while (rest_start < header.len and (header[rest_start] == ' ' or header[rest_start] == '\t')) {
                rest_start += 1;
            }
            if (rest_start < header.len) {
                // Find end of the second word (accession).
                var acc_end = rest_start;
                while (acc_end < header.len and header[acc_end] != ' ' and header[acc_end] != '\t') {
                    acc_end += 1;
                }
                accession = header[rest_start..acc_end];
                // Skip whitespace after accession to find description.
                var desc_start = acc_end;
                while (desc_start < header.len and (header[desc_start] == ' ' or header[desc_start] == '\t')) {
                    desc_start += 1;
                }
                if (desc_start < header.len) {
                    description = header[desc_start..];
                }
            }
            break;
        }
    }

    // --- Collect sequence lines until next '>' or EOF ---
    var seq_buf = std.ArrayList(u8){};
    defer seq_buf.deinit(allocator);

    while (pos.* < data.len and data[pos.*] != '>') {
        const c = data[pos.*];
        pos.* += 1;
        // Skip whitespace (newlines, spaces, tabs, carriage returns)
        if (c == '\n' or c == '\r' or c == ' ' or c == '\t') continue;
        try seq_buf.append(allocator, c);
    }

    // Digitize accumulated sequence text (empty sequences are valid in FASTA)
    const dsq = try abc.digitize(allocator, seq_buf.items);
    errdefer allocator.free(dsq);

    const name_copy = try allocator.dupe(u8, name);
    errdefer allocator.free(name_copy);

    var acc_copy: ?[]const u8 = null;
    if (accession) |acc| {
        acc_copy = try allocator.dupe(u8, acc);
    }
    errdefer if (acc_copy) |a| allocator.free(a);

    var desc_copy: ?[]const u8 = null;
    if (description) |desc| {
        desc_copy = try allocator.dupe(u8, desc);
    }
    errdefer if (desc_copy) |d| allocator.free(d);

    return Sequence{
        .name = name_copy,
        .accession = acc_copy,
        .description = desc_copy,
        .taxonomy_id = null,
        .dsq = dsq,
        .secondary_structure = null,
        .source = null,
        .abc = abc,
        .allocator = allocator,
    };
}

/// Write a single sequence in FASTA format to dest.
/// line_width controls how many residues per line (0 means no wrapping).
pub fn write(dest: std.io.AnyWriter, seq: Sequence, line_width: usize) !void {
    try dest.writeByte('>');
    try dest.writeAll(seq.name);
    if (seq.description) |desc| {
        try dest.writeByte(' ');
        try dest.writeAll(desc);
    }
    try dest.writeByte('\n');

    const text = try seq.toText();
    defer seq.allocator.free(text);

    if (line_width == 0 or text.len == 0) {
        try dest.writeAll(text);
        if (text.len > 0) try dest.writeByte('\n');
        return;
    }

    var i: usize = 0;
    while (i < text.len) {
        const end = @min(i + line_width, text.len);
        try dest.writeAll(text[i..end]);
        try dest.writeByte('\n');
        i = end;
    }
}

/// Write multiple sequences in FASTA format.
pub fn writeAll(dest: std.io.AnyWriter, seqs: []const Sequence, line_width: usize) !void {
    for (seqs) |seq| {
        try write(dest, seq, line_width);
    }
}

// --- Tests ---

const alphabet_mod = @import("../alphabet.zig");

test "parseAll: single sequence" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, seqs[0].dsq);
    try std.testing.expectEqual(@as(?[]const u8, null), seqs[0].description);
}

test "parseAll: single sequence with description" {
    const allocator = std.testing.allocator;
    const data = ">seq1 A DNA sequence\nACGT\n";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqualStrings("A", seqs[0].accession.?);
    try std.testing.expectEqualStrings("DNA sequence", seqs[0].description.?);
}

test "parseAll: accession parsed from header" {
    const allocator = std.testing.allocator;
    const data = ">seq1 ACC001 A DNA sequence\nACGT\n";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqualStrings("ACC001", seqs[0].accession.?);
    try std.testing.expectEqualStrings("A DNA sequence", seqs[0].description.?);
}

test "parseAll: accession only, no description" {
    const allocator = std.testing.allocator;
    const data = ">seq1 ACC001\nACGT\n";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqualStrings("ACC001", seqs[0].accession.?);
    try std.testing.expectEqual(@as(?[]const u8, null), seqs[0].description);
}

test "parseAll: name only, no accession" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqual(@as(?[]const u8, null), seqs[0].accession);
    try std.testing.expectEqual(@as(?[]const u8, null), seqs[0].description);
}

test "parseAll: multi-line sequence" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\nACGT\n";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqual(@as(usize, 8), seqs[0].dsq.len);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3, 0, 1, 2, 3 }, seqs[0].dsq);
}

test "parseAll: multiple sequences" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 2), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqualStrings("seq2", seqs[1].name);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, seqs[0].dsq);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 2, 2, 2, 2 }, seqs[1].dsq);
}

test "parseAll: lowercase sequence" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nacgt\n";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, seqs[0].dsq);
}

test "parseAll: empty input returns empty slice" {
    const allocator = std.testing.allocator;
    const seqs = try parseAll(allocator, &alphabet_mod.dna, "");
    defer allocator.free(seqs);

    try std.testing.expectEqual(@as(usize, 0), seqs.len);
}

test "parseAll: empty sequence accepted" {
    const allocator = std.testing.allocator;
    const seqs = try parseAll(allocator, &alphabet_mod.dna, ">seq1\n>seq2\nACGT\n");
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 2), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqual(@as(usize, 0), seqs[0].dsq.len);
    try std.testing.expectEqualStrings("seq2", seqs[1].name);
    try std.testing.expectEqual(@as(usize, 4), seqs[1].dsq.len);
}

test "parseAll: empty sequence at end" {
    const allocator = std.testing.allocator;
    const seqs = try parseAll(allocator, &alphabet_mod.dna, ">seq1\nACGT\n>seq2\n");
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 2), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqual(@as(usize, 4), seqs[0].dsq.len);
    try std.testing.expectEqualStrings("seq2", seqs[1].name);
    try std.testing.expectEqual(@as(usize, 0), seqs[1].dsq.len);
}

test "parseAll: invalid format (no '>') returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(
        error.InvalidFormat,
        parseAll(allocator, &alphabet_mod.dna, "ACGT\n"),
    );
}

test "parseAll: no trailing newline" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT";
    const seqs = try parseAll(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, seqs[0].dsq);
}

test "write: single sequence" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), seq, 60);
    try std.testing.expectEqualStrings(">seq1\nACGT\n", buf.items);
}

test "write: sequence with description" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq.deinit();
    seq.description = try allocator.dupe(u8, "my description");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), seq, 60);
    try std.testing.expectEqualStrings(">seq1 my description\nACGT\n", buf.items);
}

test "write: line wrapping at 4 chars" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "s", "ACGTACGT");
    defer seq.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), seq, 4);
    try std.testing.expectEqualStrings(">s\nACGT\nACGT\n", buf.items);
}

test "write and parseAll: round trip" {
    const allocator = std.testing.allocator;

    var seq1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGTACGT");
    defer seq1.deinit();
    seq1.description = try allocator.dupe(u8, "first sequence");

    var seq2 = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq2", "GGGGCCCC");
    defer seq2.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    const originals = [_]Sequence{ seq1, seq2 };
    try writeAll(buf.writer(allocator).any(), &originals, 60);

    const parsed = try parseAll(allocator, &alphabet_mod.dna, buf.items);
    defer {
        for (parsed) |*s| @constCast(s).deinit();
        allocator.free(parsed);
    }

    try std.testing.expectEqual(@as(usize, 2), parsed.len);
    try std.testing.expectEqualStrings("seq1", parsed[0].name);
    // "first sequence" is written as description; when re-parsed, "first" becomes accession.
    try std.testing.expectEqualStrings("first", parsed[0].accession.?);
    try std.testing.expectEqualStrings("sequence", parsed[0].description.?);
    try std.testing.expectEqualSlices(u8, seq1.dsq, parsed[0].dsq);
    try std.testing.expectEqualStrings("seq2", parsed[1].name);
    try std.testing.expectEqualSlices(u8, seq2.dsq, parsed[1].dsq);
}

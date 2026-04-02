// GenBank and EMBL/UniProt format parsers and writers.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const Sequence = @import("../sequence.zig").Sequence;

// ---------------------------------------------------------------------------
// GenBank parser
// ---------------------------------------------------------------------------

/// Parse all GenBank records from a buffer.
/// Returns an owned slice of Sequences. Caller owns the slice and each Sequence.
pub fn parseAllGenBank(allocator: Allocator, abc: *const Alphabet, data: []const u8) ![]Sequence {
    var sequences = std.ArrayList(Sequence){};
    errdefer {
        for (sequences.items) |*seq| seq.deinit();
        sequences.deinit(allocator);
    }

    var pos: usize = 0;
    while (pos < data.len) {
        // Skip blank lines before next record
        while (pos < data.len) {
            const c = data[pos];
            if (c == '\n' or c == '\r' or c == ' ' or c == '\t') {
                pos += 1;
            } else {
                break;
            }
        }
        if (pos >= data.len) break;

        // Expect "LOCUS" at start of a record
        if (!std.mem.startsWith(u8, data[pos..], "LOCUS")) return error.InvalidFormat;

        var seq = try parseOneGenBank(allocator, abc, data, &pos);
        errdefer seq.deinit();
        try sequences.append(allocator, seq);
    }

    return sequences.toOwnedSlice(allocator);
}

/// Parse a single GenBank record starting at data[pos] (which must be at "LOCUS").
/// Updates pos to point past the "//" record terminator.
pub fn parseOneGenBank(allocator: Allocator, abc: *const Alphabet, data: []const u8, pos: *usize) !Sequence {
    var name: ?[]const u8 = null;
    var accession: ?[]const u8 = null;
    // description may span continuation lines
    var desc_buf = std.ArrayList(u8){};
    defer desc_buf.deinit(allocator);

    var in_definition = false;
    var in_origin = false;
    var seq_buf = std.ArrayList(u8){};
    defer seq_buf.deinit(allocator);
    var found_terminator = false;

    while (pos.* < data.len) {
        // Read one line
        const line_start = pos.*;
        while (pos.* < data.len and data[pos.*] != '\n') {
            pos.* += 1;
        }
        const line_end = pos.*;
        // Consume the newline
        if (pos.* < data.len) pos.* += 1;

        var line = data[line_start..line_end];
        // Strip trailing CR
        if (line.len > 0 and line[line.len - 1] == '\r') {
            line = line[0 .. line.len - 1];
        }

        // Record terminator
        if (std.mem.eql(u8, line, "//")) {
            found_terminator = true;
            break;
        }

        if (in_origin) {
            // Sequence data: skip position numbers and spaces, collect letters
            for (line) |c| {
                if (std.ascii.isAlphabetic(c)) {
                    try seq_buf.append(allocator, c);
                }
            }
            continue;
        }

        if (line.len == 0) {
            in_definition = false;
            continue;
        }

        // Check for keyword lines (no leading space)
        if (line[0] != ' ') {
            in_definition = false;

            if (std.mem.startsWith(u8, line, "LOCUS")) {
                // LOCUS <name> ...
                // Name is the token right after "LOCUS" keyword
                const rest = std.mem.trimLeft(u8, line["LOCUS".len..], " \t");
                var end: usize = 0;
                while (end < rest.len and !std.ascii.isWhitespace(rest[end])) {
                    end += 1;
                }
                name = rest[0..end];
            } else if (std.mem.startsWith(u8, line, "DEFINITION")) {
                const rest = std.mem.trimLeft(u8, line["DEFINITION".len..], " \t");
                try desc_buf.appendSlice(allocator, rest);
                in_definition = true;
            } else if (std.mem.startsWith(u8, line, "VERSION")) {
                // VERSION line format: "VERSION   <accession.version>"
                // Like Easel, use VERSION for accession (includes version suffix).
                const rest = std.mem.trimLeft(u8, line["VERSION".len..], " \t");
                // Take first whitespace-delimited token
                var end: usize = 0;
                while (end < rest.len and !std.ascii.isWhitespace(rest[end])) {
                    end += 1;
                }
                if (end > 0) accession = rest[0..end];
            } else if (std.mem.startsWith(u8, line, "ORIGIN")) {
                in_origin = true;
            }
        } else {
            // Continuation line (has leading whitespace)
            if (in_definition) {
                const cont = std.mem.trimLeft(u8, line, " \t");
                if (desc_buf.items.len > 0) {
                    try desc_buf.append(allocator, ' ');
                }
                try desc_buf.appendSlice(allocator, cont);
            }
        }
    }

    if (name == null) return error.InvalidFormat;
    if (!found_terminator) return error.InvalidFormat;

    const dsq = try abc.digitize(allocator, seq_buf.items);
    errdefer allocator.free(dsq);

    const name_copy = try allocator.dupe(u8, name.?);
    errdefer allocator.free(name_copy);

    var acc_copy: ?[]const u8 = null;
    if (accession) |acc| {
        acc_copy = try allocator.dupe(u8, acc);
    }
    errdefer if (acc_copy) |acc| allocator.free(acc);

    var desc_copy: ?[]const u8 = null;
    if (desc_buf.items.len > 0) {
        desc_copy = try allocator.dupe(u8, desc_buf.items);
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

// ---------------------------------------------------------------------------
// EMBL / UniProt parser
// ---------------------------------------------------------------------------

/// Parse all EMBL/UniProt records from a buffer.
/// Returns an owned slice of Sequences. Caller owns the slice and each Sequence.
pub fn parseAllEmbl(allocator: Allocator, abc: *const Alphabet, data: []const u8) ![]Sequence {
    var sequences = std.ArrayList(Sequence){};
    errdefer {
        for (sequences.items) |*seq| seq.deinit();
        sequences.deinit(allocator);
    }

    var pos: usize = 0;
    while (pos < data.len) {
        // Skip blank lines before next record
        while (pos < data.len) {
            const c = data[pos];
            if (c == '\n' or c == '\r' or c == ' ' or c == '\t') {
                pos += 1;
            } else {
                break;
            }
        }
        if (pos >= data.len) break;

        // Expect "ID   " at start of a record (5-char prefix)
        if (!std.mem.startsWith(u8, data[pos..], "ID   ")) return error.InvalidFormat;

        var seq = try parseOneEmbl(allocator, abc, data, &pos);
        errdefer seq.deinit();
        try sequences.append(allocator, seq);
    }

    return sequences.toOwnedSlice(allocator);
}

/// Parse a single EMBL/UniProt record starting at data[pos] (which must be at "ID   ").
/// Updates pos to point past the "//" record terminator.
pub fn parseOneEmbl(allocator: Allocator, abc: *const Alphabet, data: []const u8, pos: *usize) !Sequence {
    var name: ?[]const u8 = null;
    var accession: ?[]const u8 = null;
    // description may span multiple DE continuation lines
    var desc_buf = std.ArrayList(u8){};
    defer desc_buf.deinit(allocator);
    var in_sequence = false;
    var found_terminator = false;

    var seq_buf = std.ArrayList(u8){};
    defer seq_buf.deinit(allocator);

    while (pos.* < data.len) {
        // Read one line
        const line_start = pos.*;
        while (pos.* < data.len and data[pos.*] != '\n') {
            pos.* += 1;
        }
        const line_end = pos.*;
        if (pos.* < data.len) pos.* += 1;

        var line = data[line_start..line_end];
        if (line.len > 0 and line[line.len - 1] == '\r') {
            line = line[0 .. line.len - 1];
        }

        // Record terminator
        if (std.mem.eql(u8, line, "//")) {
            found_terminator = true;
            break;
        }

        if (in_sequence) {
            // Sequence data lines: leading spaces, lowercase letters, spaces, numbers at end
            for (line) |c| {
                if (std.ascii.isAlphabetic(c)) {
                    try seq_buf.append(allocator, c);
                }
            }
            continue;
        }

        if (line.len < 2) continue;

        // Two-letter tag at start
        const tag = line[0..2];

        if (std.mem.eql(u8, tag, "ID")) {
            // ID   <name>; ...
            const rest = std.mem.trimLeft(u8, line[2..], " \t");
            var end: usize = 0;
            while (end < rest.len and rest[end] != ';' and !std.ascii.isWhitespace(rest[end])) {
                end += 1;
            }
            name = rest[0..end];
        } else if (std.mem.eql(u8, tag, "AC")) {
            // AC   <accession>; ...
            const rest = std.mem.trimLeft(u8, line[2..], " \t");
            var end: usize = 0;
            while (end < rest.len and rest[end] != ';' and !std.ascii.isWhitespace(rest[end])) {
                end += 1;
            }
            if (end > 0 and accession == null) {
                accession = rest[0..end];
            }
        } else if (std.mem.eql(u8, tag, "DE")) {
            const rest = std.mem.trimLeft(u8, line[2..], " \t");
            if (desc_buf.items.len > 0) {
                try desc_buf.append(allocator, ' ');
            }
            try desc_buf.appendSlice(allocator, rest);
        } else if (std.mem.eql(u8, tag, "SQ")) {
            in_sequence = true;
        }
    }

    if (name == null) return error.InvalidFormat;
    if (!found_terminator) return error.InvalidFormat;

    const dsq = try abc.digitize(allocator, seq_buf.items);
    errdefer allocator.free(dsq);

    const name_copy = try allocator.dupe(u8, name.?);
    errdefer allocator.free(name_copy);

    var acc_copy: ?[]const u8 = null;
    if (accession) |acc| {
        acc_copy = try allocator.dupe(u8, acc);
    }
    errdefer if (acc_copy) |acc| allocator.free(acc);

    var desc_copy: ?[]const u8 = null;
    if (desc_buf.items.len > 0) {
        desc_copy = try allocator.dupe(u8, desc_buf.items);
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

// ---------------------------------------------------------------------------
// Writers
// ---------------------------------------------------------------------------

/// Write a single sequence in GenBank format.
pub fn writeGenBank(dest: std.io.AnyWriter, seq: Sequence) !void {
    // LOCUS line
    try dest.print("LOCUS       {s}\n", .{seq.name});

    // DEFINITION line
    if (seq.description) |desc| {
        try dest.print("DEFINITION  {s}\n", .{desc});
    } else {
        try dest.writeAll("DEFINITION  .\n");
    }

    // ACCESSION line (base accession without version suffix)
    if (seq.accession) |acc| {
        // Strip version suffix (everything after last '.') for ACCESSION line
        const base_acc = if (std.mem.lastIndexOfScalar(u8, acc, '.')) |dot|
            acc[0..dot]
        else
            acc;
        try dest.print("ACCESSION   {s}\n", .{base_acc});
    } else {
        try dest.print("ACCESSION   {s}\n", .{seq.name});
    }

    // VERSION line (full accession.version)
    if (seq.accession) |acc| {
        try dest.print("VERSION     {s}\n", .{acc});
    }

    // ORIGIN section
    try dest.writeAll("ORIGIN\n");

    const text = try seq.toText();
    defer seq.allocator.free(text);

    // GenBank sequence: 60 residues per line, 10 per chunk, position prefix
    var i: usize = 0;
    while (i < text.len) {
        const line_end = @min(i + 60, text.len);
        // Position (1-indexed, right-justified in 9 chars)
        try dest.print("{d:>9} ", .{i + 1});
        var j = i;
        while (j < line_end) {
            const chunk_end = @min(j + 10, line_end);
            // Lowercase as in canonical GenBank
            var k = j;
            while (k < chunk_end) : (k += 1) {
                try dest.writeByte(std.ascii.toLower(text[k]));
            }
            j = chunk_end;
            if (j < line_end) try dest.writeByte(' ');
        }
        try dest.writeByte('\n');
        i = line_end;
    }

    try dest.writeAll("//\n");
}

/// Write a single sequence in EMBL format.
pub fn writeEmbl(dest: std.io.AnyWriter, seq: Sequence) !void {
    // ID line
    try dest.print("ID   {s};\n", .{seq.name});

    // AC line
    if (seq.accession) |acc| {
        try dest.print("AC   {s};\n", .{acc});
    } else {
        try dest.print("AC   {s};\n", .{seq.name});
    }

    // DE line
    if (seq.description) |desc| {
        try dest.print("DE   {s}\n", .{desc});
    }

    // SQ line and sequence data
    const text = try seq.toText();
    defer seq.allocator.free(text);

    try dest.print("SQ   Sequence {d} BP;\n", .{text.len});

    // EMBL sequence: 60 per line in 10-char chunks, position number right-justified
    // Build each line in a fixed buffer to avoid AnyWriter arithmetic issues
    var i: usize = 0;
    while (i < text.len) {
        const line_end = @min(i + 60, text.len);
        // Build line content: "     " + chunks separated by spaces
        var line_buf: [80]u8 = undefined;
        var pos: usize = 0;
        // 5-space indent
        @memcpy(line_buf[0..5], "     ");
        pos = 5;
        var j = i;
        while (j < line_end) {
            const chunk_end = @min(j + 10, line_end);
            var k = j;
            while (k < chunk_end) : (k += 1) {
                line_buf[pos] = std.ascii.toLower(text[k]);
                pos += 1;
            }
            j = chunk_end;
            if (j < line_end) {
                line_buf[pos] = ' ';
                pos += 1;
            }
        }
        // Pad with spaces to column 70, then write position number
        while (pos < 70) {
            line_buf[pos] = ' ';
            pos += 1;
        }
        // Write the buffered line content
        try dest.writeAll(line_buf[0..pos]);
        // Write position number and newline
        try dest.print("{d:>10}\n", .{line_end});
        i = line_end;
    }

    try dest.writeAll("//\n");
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const alphabet_mod = @import("../alphabet.zig");

test "parseAllGenBank: single record with VERSION" {
    const allocator = std.testing.allocator;
    const data =
        \\LOCUS       SCU49845     5028 bp    DNA
        \\DEFINITION  Saccharomyces cerevisiae TCP1-beta gene, partial cds.
        \\ACCESSION   U49845
        \\VERSION     U49845.1
        \\ORIGIN
        \\        1 acgt acgt
        \\//
        \\
    ;
    const seqs = try parseAllGenBank(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("SCU49845", seqs[0].name);
    // Accession comes from VERSION line, not ACCESSION
    try std.testing.expectEqualStrings("U49845.1", seqs[0].accession.?);
    try std.testing.expectEqualStrings(
        "Saccharomyces cerevisiae TCP1-beta gene, partial cds.",
        seqs[0].description.?,
    );
    try std.testing.expectEqual(@as(usize, 8), seqs[0].dsq.len);
}

test "parseAllGenBank: no VERSION line gives null accession" {
    const allocator = std.testing.allocator;
    const data =
        \\LOCUS       SEQ1
        \\ACCESSION   X12345
        \\ORIGIN
        \\        1 acgt
        \\//
        \\
    ;
    const seqs = try parseAllGenBank(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("SEQ1", seqs[0].name);
    // No VERSION line means no accession
    try std.testing.expectEqual(@as(?[]const u8, null), seqs[0].accession);
}

test "parseAllGenBank: sequence strips numbers and spaces" {
    const allocator = std.testing.allocator;
    const data =
        \\LOCUS       SEQ1
        \\ORIGIN
        \\        1 gatc ctcc
        \\       11 atat acaa
        \\//
        \\
    ;
    const seqs = try parseAllGenBank(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqual(@as(usize, 16), seqs[0].dsq.len);
}

test "parseAllGenBank: multi-record buffer" {
    const allocator = std.testing.allocator;
    const data =
        \\LOCUS       SEQ1
        \\ORIGIN
        \\        1 acgt
        \\//
        \\LOCUS       SEQ2
        \\ORIGIN
        \\        1 gggg
        \\//
        \\
    ;
    const seqs = try parseAllGenBank(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 2), seqs.len);
    try std.testing.expectEqualStrings("SEQ1", seqs[0].name);
    try std.testing.expectEqualStrings("SEQ2", seqs[1].name);
}

test "parseAllGenBank: multi-line DEFINITION" {
    const allocator = std.testing.allocator;
    const data =
        \\LOCUS       SEQ1
        \\DEFINITION  First line of description
        \\            continuation here.
        \\ORIGIN
        \\        1 acgt
        \\//
        \\
    ;
    const seqs = try parseAllGenBank(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings(
        "First line of description continuation here.",
        seqs[0].description.?,
    );
}

test "parseAllEmbl: single record" {
    const allocator = std.testing.allocator;
    const data =
        \\ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
        \\AC   X56734; S46826;
        \\DE   Trifolium repens mRNA for non-cyanogenic beta-glucosidase
        \\SQ   Sequence 8 BP;
        \\     aaacaaac
        \\//
        \\
    ;
    const seqs = try parseAllEmbl(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("X56734", seqs[0].name);
    try std.testing.expectEqualStrings("X56734", seqs[0].accession.?);
    try std.testing.expectEqualStrings(
        "Trifolium repens mRNA for non-cyanogenic beta-glucosidase",
        seqs[0].description.?,
    );
    try std.testing.expectEqual(@as(usize, 8), seqs[0].dsq.len);
}

test "parseAllEmbl: name strips trailing semicolon" {
    const allocator = std.testing.allocator;
    const data =
        \\ID   MYSEQ; SV 1;
        \\AC   ACC001;
        \\SQ   Sequence 4 BP;
        \\     acgt
        \\//
        \\
    ;
    const seqs = try parseAllEmbl(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings("MYSEQ", seqs[0].name);
    try std.testing.expectEqualStrings("ACC001", seqs[0].accession.?);
}

test "parseAllEmbl: sequence strips numbers and spaces" {
    const allocator = std.testing.allocator;
    const data =
        \\ID   SEQ1;
        \\SQ   Sequence 20 BP;
        \\     aaacaaacca aatatggatt        20
        \\//
        \\
    ;
    const seqs = try parseAllEmbl(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqual(@as(usize, 20), seqs[0].dsq.len);
}

test "parseAllEmbl: multi-record buffer" {
    const allocator = std.testing.allocator;
    const data =
        \\ID   SEQ1;
        \\SQ   Sequence 4 BP;
        \\     acgt
        \\//
        \\ID   SEQ2;
        \\SQ   Sequence 4 BP;
        \\     gggg
        \\//
        \\
    ;
    const seqs = try parseAllEmbl(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 2), seqs.len);
    try std.testing.expectEqualStrings("SEQ1", seqs[0].name);
    try std.testing.expectEqualStrings("SEQ2", seqs[1].name);
}

test "parseAllEmbl: multi-line DE description" {
    const allocator = std.testing.allocator;
    const data =
        \\ID   P12345; SV 1;
        \\AC   P12345;
        \\DE   RecName: Full=Glucosidase;
        \\DE   Short=Beta-Glu;
        \\SQ   Sequence 4 BP;
        \\     acgt
        \\//
        \\
    ;
    const seqs = try parseAllEmbl(allocator, &alphabet_mod.dna, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 1), seqs.len);
    try std.testing.expectEqualStrings(
        "RecName: Full=Glucosidase; Short=Beta-Glu;",
        seqs[0].description.?,
    );
}

test "writeGenBank: basic structure" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "SEQ1", "ACGT");
    defer seq.deinit();
    seq.accession = try allocator.dupe(u8, "ACC001");
    seq.description = try allocator.dupe(u8, "A test sequence");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeGenBank(buf.writer(allocator).any(), seq);
    const out = buf.items;

    // Check structural elements
    try std.testing.expect(std.mem.indexOf(u8, out, "LOCUS       SEQ1") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "DEFINITION  A test sequence") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "ACCESSION   ACC001") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "VERSION     ACC001") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "ORIGIN") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "acgt") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "//") != null);
}

test "writeEmbl: basic structure" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "X56734", "ACGT");
    defer seq.deinit();
    seq.accession = try allocator.dupe(u8, "X56734");
    seq.description = try allocator.dupe(u8, "Test sequence");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEmbl(buf.writer(allocator).any(), seq);
    const out = buf.items;

    try std.testing.expect(std.mem.indexOf(u8, out, "ID   X56734;") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "AC   X56734;") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "DE   Test sequence") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "SQ   Sequence 4 BP;") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "acgt") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "//") != null);
}

test "GenBank round-trip: write then parse" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "SEQ1", "ACGTACGT");
    defer seq.deinit();
    seq.accession = try allocator.dupe(u8, "U12345.1");
    seq.description = try allocator.dupe(u8, "Round trip test");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);
    try writeGenBank(buf.writer(allocator).any(), seq);

    // Verify VERSION line is written
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "VERSION     U12345.1") != null);
    // Verify ACCESSION line strips version suffix
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "ACCESSION   U12345\n") != null);

    const parsed = try parseAllGenBank(allocator, &alphabet_mod.dna, buf.items);
    defer {
        for (parsed) |*s| @constCast(s).deinit();
        allocator.free(parsed);
    }

    try std.testing.expectEqual(@as(usize, 1), parsed.len);
    try std.testing.expectEqualStrings("SEQ1", parsed[0].name);
    // Round-trip: accession comes from VERSION line
    try std.testing.expectEqualStrings("U12345.1", parsed[0].accession.?);
    try std.testing.expectEqualStrings("Round trip test", parsed[0].description.?);
    try std.testing.expectEqualSlices(u8, seq.dsq, parsed[0].dsq);
}

test "EMBL round-trip: write then parse" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "X56734", "ACGTACGT");
    defer seq.deinit();
    seq.accession = try allocator.dupe(u8, "X56734");
    seq.description = try allocator.dupe(u8, "EMBL round trip");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);
    try writeEmbl(buf.writer(allocator).any(), seq);

    const parsed = try parseAllEmbl(allocator, &alphabet_mod.dna, buf.items);
    defer {
        for (parsed) |*s| @constCast(s).deinit();
        allocator.free(parsed);
    }

    try std.testing.expectEqual(@as(usize, 1), parsed.len);
    try std.testing.expectEqualStrings("X56734", parsed[0].name);
    try std.testing.expectEqualStrings("X56734", parsed[0].accession.?);
    try std.testing.expectEqualSlices(u8, seq.dsq, parsed[0].dsq);
}

test "writeEmbl: exactly 60 residues" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.amino, "test60", "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY");
    defer seq.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEmbl(buf.writer(allocator).any(), seq);
    const out = buf.items;
    try std.testing.expect(std.mem.indexOf(u8, out, "//") != null);
}

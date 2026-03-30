// Unified sequence reader with format auto-detection.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const Sequence = @import("../sequence.zig").Sequence;
const fasta = @import("fasta.zig");
const stockholm = @import("stockholm.zig");

pub const Format = enum {
    fasta,
    stockholm,

    /// Detect format from the first non-whitespace bytes of data.
    pub fn detect(header: []const u8) ?Format {
        var i: usize = 0;
        while (i < header.len and (header[i] == ' ' or header[i] == '\t' or header[i] == '\n' or header[i] == '\r')) {
            i += 1;
        }
        if (i < header.len and header[i] == '>') return .fasta;
        // Detect Stockholm: first non-blank line starts with "# STOCKHOLM"
        const rest = header[i..];
        if (std.mem.startsWith(u8, rest, "# STOCKHOLM")) return .stockholm;
        return null;
    }
};

pub const Reader = struct {
    format: Format,
    data: []const u8,
    pos: usize,
    abc: *const Alphabet,
    allocator: Allocator,
    owns_data: bool,
    // Stockholm buffering: sequences are extracted from the MSA on first call.
    stk_seqs: ?[]Sequence = null,
    stk_idx: usize = 0,

    /// Open a file and create a reader (reads entire file into memory).
    /// If format is null, the format is auto-detected from the file header.
    pub fn fromFile(allocator: Allocator, abc: *const Alphabet, path: []const u8, format: ?Format) !Reader {
        const file = try std.fs.cwd().openFile(path, .{});
        defer file.close();

        const data = try file.readToEndAlloc(allocator, std.math.maxInt(usize));
        errdefer allocator.free(data);

        const fmt = format orelse Format.detect(data) orelse return error.UnknownFormat;

        return Reader{
            .format = fmt,
            .data = data,
            .pos = 0,
            .abc = abc,
            .allocator = allocator,
            .owns_data = true,
        };
    }

    /// Create a reader from in-memory data (does not take ownership of data).
    /// If format is null, the format is auto-detected from the data header.
    pub fn fromMemory(allocator: Allocator, abc: *const Alphabet, data: []const u8, format: ?Format) !Reader {
        const fmt = format orelse Format.detect(data) orelse return error.UnknownFormat;

        return Reader{
            .format = fmt,
            .data = data,
            .pos = 0,
            .abc = abc,
            .allocator = allocator,
            .owns_data = false,
        };
    }

    /// Read the next sequence record. Returns null at EOF.
    pub fn next(self: *Reader) !?Sequence {
        switch (self.format) {
            .fasta => return try nextFasta(self),
            .stockholm => return try nextStockholm(self),
        }
    }

    fn nextFasta(self: *Reader) !?Sequence {
        // Skip blank lines
        while (self.pos < self.data.len) {
            const c = self.data[self.pos];
            if (c == '\n' or c == '\r' or c == ' ' or c == '\t') {
                self.pos += 1;
            } else {
                break;
            }
        }
        if (self.pos >= self.data.len) return null;
        if (self.data[self.pos] != '>') return error.InvalidFormat;

        return try parseFastaRecord(self.allocator, self.abc, self.data, &self.pos);
    }

    /// Read all remaining sequences.
    pub fn readAll(self: *Reader) ![]Sequence {
        var list = std.ArrayList(Sequence){};
        errdefer {
            for (list.items) |*seq| seq.deinit();
            list.deinit(self.allocator);
        }

        while (try self.next()) |seq| {
            var s = seq;
            errdefer s.deinit();
            try list.append(self.allocator, s);
        }

        return list.toOwnedSlice(self.allocator);
    }

    fn nextStockholm(self: *Reader) !?Sequence {
        // On first call, parse the entire MSA and extract all sequences.
        if (self.stk_seqs == null) {
            var msa = try stockholm.parse(self.allocator, self.abc, self.data);
            errdefer msa.deinit();

            const n = msa.nseq();
            const seqs = try self.allocator.alloc(Sequence, n);
            errdefer self.allocator.free(seqs);

            var done: usize = 0;
            errdefer for (0..done) |i| seqs[i].deinit();

            for (0..n) |i| {
                seqs[i] = try msa.extractSeq(i);
                done += 1;
            }

            msa.deinit();
            self.stk_seqs = seqs;
            self.stk_idx = 0;
            // Mark data as consumed so pos-based calls are harmless.
            self.pos = self.data.len;
        }

        const seqs = self.stk_seqs.?;
        if (self.stk_idx >= seqs.len) return null;

        const seq = seqs[self.stk_idx];
        self.stk_idx += 1;
        return seq;
    }

    pub fn deinit(self: *Reader) void {
        if (self.stk_seqs) |seqs| {
            // Sequences already returned to caller are not freed here;
            // we only free any that were never consumed.
            for (seqs[self.stk_idx..]) |*s| @constCast(s).deinit();
            self.allocator.free(seqs);
        }
        if (self.owns_data) {
            self.allocator.free(self.data);
        }
    }
};

/// Parse one FASTA record from data starting at pos (must point to '>').
/// Advances pos past the record.
fn parseFastaRecord(allocator: Allocator, abc: *const Alphabet, data: []const u8, pos: *usize) !Sequence {
    // Consume '>'
    pos.* += 1;

    // Parse header line
    const header_start = pos.*;
    while (pos.* < data.len and data[pos.*] != '\n' and data[pos.*] != '\r') {
        pos.* += 1;
    }
    const header = data[header_start..pos.*];

    if (pos.* < data.len and data[pos.*] == '\r') pos.* += 1;
    if (pos.* < data.len and data[pos.*] == '\n') pos.* += 1;

    var name: []const u8 = header;
    var description: ?[]const u8 = null;
    for (header, 0..) |c, i| {
        if (c == ' ' or c == '\t') {
            name = header[0..i];
            var ds = i + 1;
            while (ds < header.len and (header[ds] == ' ' or header[ds] == '\t')) ds += 1;
            if (ds < header.len) description = header[ds..];
            break;
        }
    }

    // Collect sequence characters
    var seq_buf = std.ArrayList(u8){};
    defer seq_buf.deinit(allocator);

    while (pos.* < data.len and data[pos.*] != '>') {
        const c = data[pos.*];
        pos.* += 1;
        if (c == '\n' or c == '\r' or c == ' ' or c == '\t') continue;
        try seq_buf.append(allocator, c);
    }

    const dsq = try abc.digitize(allocator, seq_buf.items);
    errdefer allocator.free(dsq);

    const name_copy = try allocator.dupe(u8, name);
    errdefer allocator.free(name_copy);

    var desc_copy: ?[]const u8 = null;
    if (description) |desc| {
        desc_copy = try allocator.dupe(u8, desc);
    }
    errdefer if (desc_copy) |d| allocator.free(d);

    return Sequence{
        .name = name_copy,
        .accession = null,
        .description = desc_copy,
        .taxonomy_id = null,
        .dsq = dsq,
        .secondary_structure = null,
        .source = null,
        .abc = abc,
        .allocator = allocator,
    };
}

// --- Tests ---

const alphabet_mod = @import("../alphabet.zig");

test "Format.detect: fasta" {
    try std.testing.expectEqual(@as(?Format, .fasta), Format.detect(">seq1\nACGT\n"));
}

test "Format.detect: unknown returns null" {
    try std.testing.expectEqual(@as(?Format, null), Format.detect("ACGT\n"));
}

test "Format.detect: stockholm" {
    try std.testing.expectEqual(@as(?Format, .stockholm), Format.detect("# STOCKHOLM 1.0\n"));
}

test "Format.detect: stockholm with leading newline" {
    try std.testing.expectEqual(@as(?Format, .stockholm), Format.detect("\n# STOCKHOLM 1.0\n"));
}

test "Format.detect: fasta with leading whitespace" {
    try std.testing.expectEqual(@as(?Format, .fasta), Format.detect("\n>seq1\nACGT\n"));
}

test "Reader.fromMemory: auto-detect fasta" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n";

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, data, null);
    defer reader.deinit();

    try std.testing.expectEqual(Format.fasta, reader.format);
}

test "Reader.fromMemory: unknown format returns error" {
    const allocator = std.testing.allocator;
    const data = "ACGT\n";

    try std.testing.expectError(
        error.UnknownFormat,
        Reader.fromMemory(allocator, &alphabet_mod.dna, data, null),
    );
}

test "Reader.next: iterator pattern" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, data, .fasta);
    defer reader.deinit();

    var seq1 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq1.deinit();
    try std.testing.expectEqualStrings("seq1", seq1.name);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, seq1.dsq);

    var seq2 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq2.deinit();
    try std.testing.expectEqualStrings("seq2", seq2.name);

    const eof = try reader.next();
    try std.testing.expectEqual(@as(?Sequence, null), eof);
}

test "Reader.readAll" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, data, .fasta);
    defer reader.deinit();

    const seqs = try reader.readAll();
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 2), seqs.len);
    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqualStrings("seq2", seqs[1].name);
}

test "Reader.next: stockholm yields ungapped sequences" {
    const allocator = std.testing.allocator;
    const data = "# STOCKHOLM 1.0\n\nseq1  AC-GT\nseq2  ACGGT\n//\n";

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, data, null);
    defer reader.deinit();

    try std.testing.expectEqual(Format.stockholm, reader.format);

    var seq1 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq1.deinit();
    try std.testing.expectEqualStrings("seq1", seq1.name);
    // Gaps removed: AC-GT -> ACGT (4 residues)
    try std.testing.expectEqual(@as(usize, 4), seq1.dsq.len);

    var seq2 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq2.deinit();
    try std.testing.expectEqualStrings("seq2", seq2.name);
    try std.testing.expectEqual(@as(usize, 5), seq2.dsq.len);

    const eof = try reader.next();
    try std.testing.expectEqual(@as(?Sequence, null), eof);
}

test "Reader.readAll: second call on exhausted reader returns empty" {
    const allocator = std.testing.allocator;

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, ">s\nACGT\n", .fasta);
    defer reader.deinit();

    const seqs = try reader.readAll();
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }
    try std.testing.expectEqual(@as(usize, 1), seqs.len);

    const more = try reader.readAll();
    defer allocator.free(more);
    try std.testing.expectEqual(@as(usize, 0), more.len);
}

// Unified sequence writer.

const std = @import("std");
const Sequence = @import("../sequence.zig").Sequence;
const fasta = @import("fasta.zig");
const Format = @import("reader.zig").Format;

pub const Writer = struct {
    format: Format,
    dest: std.io.AnyWriter,
    line_width: usize,

    /// Create a Writer for the given format.
    /// line_width controls residues per line when writing FASTA (60 is standard).
    pub fn init(dest: std.io.AnyWriter, format: Format, line_width: usize) Writer {
        return Writer{
            .format = format,
            .dest = dest,
            .line_width = line_width,
        };
    }

    /// Write a single sequence record.
    /// Note: for stockholm format, use stockholm.write() directly with an Msa —
    /// writing a single ungapped Sequence in Stockholm does not make semantic sense.
    /// This path writes the sequence in FASTA as a fallback for stockholm format.
    pub fn write(self: *Writer, seq: Sequence) !void {
        switch (self.format) {
            .fasta => try fasta.write(self.dest, seq, self.line_width),
            .stockholm => try fasta.write(self.dest, seq, self.line_width),
        }
    }

    /// Write multiple sequence records.
    pub fn writeAll(self: *Writer, seqs: []const Sequence) !void {
        for (seqs) |seq| {
            try self.write(seq);
        }
    }
};

// --- Tests ---

const alphabet_mod = @import("../alphabet.zig");

test "Writer.write: basic FASTA output" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "myseq", "ACGTACGT");
    defer seq.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    var w = Writer.init(buf.writer(allocator).any(), .fasta, 60);
    try w.write(seq);

    try std.testing.expectEqualStrings(">myseq\nACGTACGT\n", buf.items);
}

test "Writer.writeAll: two sequences" {
    const allocator = std.testing.allocator;
    var seq1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "s1", "ACGT");
    defer seq1.deinit();
    var seq2 = try Sequence.fromText(allocator, &alphabet_mod.dna, "s2", "TTTT");
    defer seq2.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    var w = Writer.init(buf.writer(allocator).any(), .fasta, 60);
    const seqs = [_]Sequence{ seq1, seq2 };
    try w.writeAll(&seqs);

    try std.testing.expectEqualStrings(">s1\nACGT\n>s2\nTTTT\n", buf.items);
}

test "Writer.write: line wrapping" {
    const allocator = std.testing.allocator;
    var seq = try Sequence.fromText(allocator, &alphabet_mod.dna, "s", "ACGTACGT");
    defer seq.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    var w = Writer.init(buf.writer(allocator).any(), .fasta, 4);
    try w.write(seq);

    try std.testing.expectEqualStrings(">s\nACGT\nACGT\n", buf.items);
}

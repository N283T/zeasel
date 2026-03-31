// Unified sequence reader with format auto-detection.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const Sequence = @import("../sequence.zig").Sequence;
const fasta = @import("fasta.zig");
const stockholm = @import("stockholm.zig");
const genbank = @import("genbank.zig");
const clustal = @import("clustal.zig");
const phylip = @import("phylip.zig");
const a2m = @import("a2m.zig");
const psiblast = @import("psiblast.zig");
const selex = @import("selex.zig");

pub const Format = enum {
    fasta,
    stockholm,
    genbank,
    embl,
    /// Clustal alignment format (CLUSTAL W / CLUSTAL OMEGA).
    clustal,
    /// Aligned FASTA — syntactically identical to FASTA but all sequences
    /// have the same length including gap columns. Cannot be auto-detected;
    /// must be specified explicitly.
    afa,
    /// PHYLIP interleaved format.
    phylip,
    /// A2M (aligned FASTA with insert annotation via case).
    /// Cannot be auto-detected (looks like FASTA); must be specified explicitly.
    a2m,
    /// PSI-BLAST flat text alignment format.
    /// Cannot be auto-detected; must be specified explicitly.
    psiblast,
    /// SELEX (old alignment format from Sean Eddy).
    /// Cannot be auto-detected; must be specified explicitly.
    selex,

    /// Detect format from the first non-whitespace bytes of data.
    pub fn detect(header: []const u8) ?Format {
        var i: usize = 0;
        while (i < header.len and (header[i] == ' ' or header[i] == '\t' or header[i] == '\n' or header[i] == '\r')) {
            i += 1;
        }
        if (i >= header.len) return null;
        if (header[i] == '>') return .fasta;
        const rest = header[i..];
        if (std.mem.startsWith(u8, rest, "# STOCKHOLM")) return .stockholm;
        if (std.mem.startsWith(u8, rest, "LOCUS")) return .genbank;
        // EMBL/UniProt: two-letter tag "ID" followed by exactly 3 spaces
        if (std.mem.startsWith(u8, rest, "ID   ")) return .embl;
        if (std.mem.startsWith(u8, rest, "CLUSTAL")) return .clustal;
        // CLUSTAL-like: MUSCLE, PROBCONS etc. output with "multiple sequence alignment"
        if (std.mem.indexOf(u8, rest[0..@min(rest.len, 80)], "multiple sequence alignment") != null) return .clustal;
        // PHYLIP: first non-whitespace content is two integers (nseq alen)
        if (detectPhylip(rest)) return .phylip;
        return null;
    }

    fn detectPhylip(data: []const u8) bool {
        // PHYLIP starts with two whitespace-separated integers on the first line
        var it = std.mem.tokenizeAny(u8, data, " \t");
        const tok1 = it.next() orelse return false;
        // First token must end before newline and be all digits
        for (tok1) |c| {
            if (c < '0' or c > '9') return false;
        }
        const tok2 = it.next() orelse return false;
        for (tok2) |c| {
            if (c == '\n' or c == '\r') break;
            if (c < '0' or c > '9') return false;
        }
        return true;
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
    // GenBank/EMBL buffering: all records parsed on first call, then iterated.
    gb_seqs: ?[]Sequence = null,
    gb_idx: usize = 0,
    // Clustal/AFA buffering: sequences extracted from the MSA on first call.
    msa_seqs: ?[]Sequence = null,
    msa_idx: usize = 0,

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
            .genbank => return try nextGenBank(self),
            .embl => return try nextEmbl(self),
            .clustal => return try nextClustal(self),
            .afa => return try nextAfa(self),
            .phylip => return try nextPhylip(self),
            .a2m => return try nextA2m(self),
            .psiblast => return try nextPsiblast(self),
            .selex => return try nextSelex(self),
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

        return try fasta.parseOne(self.allocator, self.abc, self.data, &self.pos);
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

    fn nextGenBank(self: *Reader) !?Sequence {
        if (self.gb_seqs == null) {
            const seqs = try genbank.parseAllGenBank(self.allocator, self.abc, self.data);
            self.gb_seqs = seqs;
            self.gb_idx = 0;
            self.pos = self.data.len;
        }
        const seqs = self.gb_seqs.?;
        if (self.gb_idx >= seqs.len) return null;
        const seq = seqs[self.gb_idx];
        self.gb_idx += 1;
        return seq;
    }

    fn nextEmbl(self: *Reader) !?Sequence {
        if (self.gb_seqs == null) {
            const seqs = try genbank.parseAllEmbl(self.allocator, self.abc, self.data);
            self.gb_seqs = seqs;
            self.gb_idx = 0;
            self.pos = self.data.len;
        }
        const seqs = self.gb_seqs.?;
        if (self.gb_idx >= seqs.len) return null;
        const seq = seqs[self.gb_idx];
        self.gb_idx += 1;
        return seq;
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

    fn nextMsaFormat(self: *Reader, comptime parseFn: anytype) !?Sequence {
        if (self.msa_seqs == null) {
            var msa = try parseFn(self.allocator, self.abc, self.data);
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
            self.msa_seqs = seqs;
            self.msa_idx = 0;
            self.pos = self.data.len;
        }

        const seqs = self.msa_seqs.?;
        if (self.msa_idx >= seqs.len) return null;

        const seq = seqs[self.msa_idx];
        self.msa_idx += 1;
        return seq;
    }

    fn nextClustal(self: *Reader) !?Sequence {
        return self.nextMsaFormat(clustal.parse);
    }

    fn nextAfa(self: *Reader) !?Sequence {
        const afa_mod = @import("afa.zig");
        return self.nextMsaFormat(afa_mod.parse);
    }

    fn nextPhylip(self: *Reader) !?Sequence {
        return self.nextMsaFormat(phylip.parse);
    }

    fn nextA2m(self: *Reader) !?Sequence {
        return self.nextMsaFormat(a2m.parse);
    }

    fn nextPsiblast(self: *Reader) !?Sequence {
        return self.nextMsaFormat(psiblast.parse);
    }

    fn nextSelex(self: *Reader) !?Sequence {
        return self.nextMsaFormat(selex.parse);
    }

    pub fn deinit(self: *Reader) void {
        if (self.stk_seqs) |seqs| {
            // Sequences already returned to caller are not freed here;
            // we only free any that were never consumed.
            for (seqs[self.stk_idx..]) |*s| @constCast(s).deinit();
            self.allocator.free(seqs);
        }
        if (self.gb_seqs) |seqs| {
            for (seqs[self.gb_idx..]) |*s| @constCast(s).deinit();
            self.allocator.free(seqs);
        }
        if (self.msa_seqs) |seqs| {
            for (seqs[self.msa_idx..]) |*s| @constCast(s).deinit();
            self.allocator.free(seqs);
        }
        if (self.owns_data) {
            self.allocator.free(self.data);
        }
    }
};

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

test "Format.detect: genbank" {
    try std.testing.expectEqual(@as(?Format, .genbank), Format.detect("LOCUS       SEQ1\n"));
}

test "Format.detect: embl" {
    try std.testing.expectEqual(@as(?Format, .embl), Format.detect("ID   X56734; SV 1;\n"));
}

test "Reader.next: genbank auto-detect and iterate" {
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

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, data, null);
    defer reader.deinit();

    try std.testing.expectEqual(Format.genbank, reader.format);

    var seq1 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq1.deinit();
    try std.testing.expectEqualStrings("SEQ1", seq1.name);
    try std.testing.expectEqual(@as(usize, 4), seq1.dsq.len);

    var seq2 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq2.deinit();
    try std.testing.expectEqualStrings("SEQ2", seq2.name);

    const eof = try reader.next();
    try std.testing.expectEqual(@as(?Sequence, null), eof);
}

test "Reader.next: embl auto-detect and iterate" {
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

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, data, null);
    defer reader.deinit();

    try std.testing.expectEqual(Format.embl, reader.format);

    var seq1 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq1.deinit();
    try std.testing.expectEqualStrings("SEQ1", seq1.name);

    var seq2 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq2.deinit();
    try std.testing.expectEqualStrings("SEQ2", seq2.name);

    const eof = try reader.next();
    try std.testing.expectEqual(@as(?Sequence, null), eof);
}

test "Format.detect: clustal" {
    try std.testing.expectEqual(@as(?Format, .clustal), Format.detect("CLUSTAL W (1.83) multiple sequence alignment\n"));
}

test "Reader.next: clustal auto-detect yields ungapped sequences" {
    const allocator = std.testing.allocator;
    const data =
        \\CLUSTAL W (1.83) multiple sequence alignment
        \\
        \\seq1      AC-GT
        \\seq2      ACGGT
        \\
        \\
    ;

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, data, null);
    defer reader.deinit();

    try std.testing.expectEqual(Format.clustal, reader.format);

    var seq1 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq1.deinit();
    try std.testing.expectEqualStrings("seq1", seq1.name);
    // Gaps removed: AC-GT -> ACGT (4 residues).
    try std.testing.expectEqual(@as(usize, 4), seq1.dsq.len);

    var seq2 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq2.deinit();
    try std.testing.expectEqualStrings("seq2", seq2.name);
    try std.testing.expectEqual(@as(usize, 5), seq2.dsq.len);

    const eof = try reader.next();
    try std.testing.expectEqual(@as(?Sequence, null), eof);
}

test "Reader.next: afa explicit format yields ungapped sequences" {
    const allocator = std.testing.allocator;
    const data =
        \\>seq1
        \\AC-GT
        \\>seq2
        \\ACGGT
        \\
    ;

    var reader = try Reader.fromMemory(allocator, &alphabet_mod.dna, data, .afa);
    defer reader.deinit();

    try std.testing.expectEqual(Format.afa, reader.format);

    var seq1 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq1.deinit();
    try std.testing.expectEqualStrings("seq1", seq1.name);
    // Gaps removed by extractSeq: AC-GT -> ACGT (4 residues).
    try std.testing.expectEqual(@as(usize, 4), seq1.dsq.len);

    var seq2 = (try reader.next()) orelse return error.ExpectedSequence;
    defer seq2.deinit();
    try std.testing.expectEqualStrings("seq2", seq2.name);
    try std.testing.expectEqual(@as(usize, 5), seq2.dsq.len);

    const eof = try reader.next();
    try std.testing.expectEqual(@as(?Sequence, null), eof);
}

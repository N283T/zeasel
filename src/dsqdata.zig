// Binary sequence database format for fast sequential reading.
// Layout:
//   Header (32 bytes): magic, version, alphabet_type, padding, num_sequences, total_residues, reserved
//   Per-sequence records: name_len, seq_len, name bytes, dsq bytes

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("alphabet.zig").Alphabet;
const alphabet_mod = @import("alphabet.zig");
const Sequence = @import("sequence.zig").Sequence;

pub const MAGIC = "ZSQD".*;
pub const VERSION: u32 = 1;

pub const Header = struct {
    alphabet_type: alphabet_mod.AlphabetType,
    num_sequences: u64,
    total_residues: u64,
};

/// Write sequences to a binary dsqdata stream.
pub fn write(dest: std.io.AnyWriter, abc: *const Alphabet, seqs: []const Sequence) !void {
    // Header: magic (4) + version (4) + alphabet_type (1) + padding (3) +
    //         num_sequences (8) + total_residues (8) + reserved (4) = 32 bytes
    try dest.writeAll(&MAGIC);
    try dest.writeInt(u32, VERSION, .little);
    try dest.writeByte(@intFromEnum(abc.kind));
    try dest.writeAll(&[_]u8{ 0, 0, 0 }); // padding

    var total_res: u64 = 0;
    for (seqs) |s| total_res += @intCast(s.len());
    try dest.writeInt(u64, @intCast(seqs.len), .little);
    try dest.writeInt(u64, total_res, .little);
    try dest.writeAll(&[_]u8{ 0, 0, 0, 0 }); // reserved

    // Sequence records
    for (seqs) |seq| {
        try dest.writeInt(u32, @intCast(seq.name.len), .little);
        try dest.writeInt(u64, @intCast(seq.dsq.len), .little);
        try dest.writeAll(seq.name);
        try dest.writeAll(seq.dsq);
    }
}

/// Read the 32-byte header from a dsqdata stream.
pub fn readHeader(src: std.io.AnyReader) !Header {
    var magic: [4]u8 = undefined;
    _ = try src.readAll(&magic);
    if (!std.mem.eql(u8, &magic, &MAGIC)) return error.InvalidFormat;

    const version = try src.readInt(u32, .little);
    if (version != VERSION) return error.UnsupportedVersion;

    const abc_byte = try src.readByte();

    var padding: [3]u8 = undefined;
    _ = try src.readAll(&padding);

    const num_seqs = try src.readInt(u64, .little);
    const total_res = try src.readInt(u64, .little);

    var reserved: [4]u8 = undefined;
    _ = try src.readAll(&reserved);

    return Header{
        .alphabet_type = @enumFromInt(abc_byte),
        .num_sequences = num_seqs,
        .total_residues = total_res,
    };
}

/// Read the next sequence from a dsqdata stream.
/// Returns null at end of stream. Caller owns the returned Sequence.
pub fn readNext(allocator: Allocator, src: std.io.AnyReader, abc: *const Alphabet) !?Sequence {
    const name_len = src.readInt(u32, .little) catch |e| {
        if (e == error.EndOfStream) return null;
        return e;
    };
    const seq_len = try src.readInt(u64, .little);

    const name = try allocator.alloc(u8, name_len);
    errdefer allocator.free(name);
    _ = try src.readAll(name);

    const dsq = try allocator.alloc(u8, seq_len);
    errdefer allocator.free(dsq);
    _ = try src.readAll(dsq);

    return Sequence{
        .name = name,
        .accession = null,
        .description = null,
        .taxonomy_id = null,
        .dsq = dsq,
        .secondary_structure = null,
        .source = null,
        .abc = abc,
        .allocator = allocator,
    };
}

/// Read all sequences from a dsqdata stream (after the header has been read).
/// Caller owns the returned slice and each Sequence in it.
pub fn readAll(allocator: Allocator, src: std.io.AnyReader, abc: *const Alphabet) ![]Sequence {
    var seqs = std.ArrayList(Sequence){};
    errdefer {
        for (seqs.items) |*s| s.deinit();
        seqs.deinit(allocator);
    }
    while (try readNext(allocator, src, abc)) |seq| {
        try seqs.append(allocator, seq);
    }
    return seqs.toOwnedSlice(allocator);
}

// --- Threaded prefetch reader ---

const SequenceBlock = @import("sequence.zig").SequenceBlock;
const WorkQueue = @import("threads.zig").WorkQueue;

/// A threaded dsqdata reader that prefetches sequence blocks in a background
/// loader thread, handing them to the main thread via a work queue.
/// Must be heap-allocated (use create/destroy) because the loader thread
/// holds a pointer to the queue.
pub const PrefetchReader = struct {
    queue: WorkQueue(*SequenceBlock),
    loader: ?std.Thread,
    allocator: Allocator,

    /// Create a PrefetchReader on the heap and spawn the loader thread.
    /// The header must already have been read from src.
    /// Caller must call destroy() when done.
    pub fn create(allocator: Allocator, src: std.io.AnyReader, abc: *const Alphabet, block_size: usize) !*PrefetchReader {
        const self = try allocator.create(PrefetchReader);
        errdefer allocator.destroy(self);

        self.* = PrefetchReader{
            .queue = try WorkQueue(*SequenceBlock).init(allocator, 64),
            .loader = null,
            .allocator = allocator,
        };

        const ctx = try allocator.create(LoaderCtx);
        errdefer allocator.destroy(ctx);
        ctx.* = LoaderCtx{
            .reader = src,
            .queue = &self.queue,
            .allocator = allocator,
            .abc = abc,
            .block_size = block_size,
        };

        self.loader = try std.Thread.spawn(.{}, loaderThread, .{ctx});
        return self;
    }

    /// Read the next block of sequences. Returns null when all sequences
    /// have been read. Caller should call recycle() on each returned block.
    pub fn read(self: *PrefetchReader) ?*SequenceBlock {
        return self.queue.receive();
    }

    /// Recycle a used block by freeing its sequences and the block itself.
    pub fn recycle(self: *PrefetchReader, block: *SequenceBlock) void {
        block.deinit();
        self.allocator.destroy(block);
    }

    /// Destroy the reader, waiting for the loader thread to finish.
    pub fn destroy(self: *PrefetchReader) void {
        if (self.loader) |t| t.join();
        self.queue.deinit();
        self.allocator.destroy(self);
    }

    const LoaderCtx = struct {
        reader: std.io.AnyReader,
        queue: *WorkQueue(*SequenceBlock),
        allocator: Allocator,
        abc: *const Alphabet,
        block_size: usize,
    };

    fn loaderThread(ctx: *LoaderCtx) void {
        const alloc = ctx.allocator;
        const q = ctx.queue;
        const abc = ctx.abc;
        const bs = ctx.block_size;
        const rdr = ctx.reader;
        alloc.destroy(ctx);

        while (true) {
            const block = alloc.create(SequenceBlock) catch {
                q.close();
                return;
            };
            block.* = SequenceBlock.initCapacity(alloc, bs) catch {
                alloc.destroy(block);
                q.close();
                return;
            };

            var count: usize = 0;
            while (count < bs) {
                const seq = readNext(alloc, rdr, abc) catch break;
                if (seq) |s| {
                    block.add(s) catch break;
                    count += 1;
                } else break;
            }

            if (count == 0) {
                block.deinit();
                alloc.destroy(block);
                q.close();
                return;
            }

            q.send(block) catch {
                block.deinit();
                alloc.destroy(block);
                q.close();
                return;
            };

            if (count < bs) {
                q.close();
                return;
            }
        }
    }
};

// --- Tests ---

test "write and readAll: round trip with 3 DNA sequences" {
    const allocator = std.testing.allocator;

    var seq1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq1", "ACGT");
    defer seq1.deinit();
    var seq2 = try Sequence.fromText(allocator, &alphabet_mod.dna, "seq2", "TTTTAAAA");
    defer seq2.deinit();
    var seq3 = try Sequence.fromText(allocator, &alphabet_mod.dna, "longname", "GCGCGC");
    defer seq3.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    const originals = [_]Sequence{ seq1, seq2, seq3 };
    try write(buf.writer(allocator).any(), &alphabet_mod.dna, &originals);

    var fbs = std.io.fixedBufferStream(buf.items);
    const hdr = try readHeader(fbs.reader().any());
    const seqs = try readAll(allocator, fbs.reader().any(), &alphabet_mod.dna);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    try std.testing.expectEqual(@as(usize, 3), seqs.len);
    try std.testing.expectEqual(hdr.num_sequences, 3);
    try std.testing.expectEqual(hdr.total_residues, 18); // 4 + 8 + 6

    try std.testing.expectEqualStrings("seq1", seqs[0].name);
    try std.testing.expectEqualSlices(u8, seq1.dsq, seqs[0].dsq);

    try std.testing.expectEqualStrings("seq2", seqs[1].name);
    try std.testing.expectEqualSlices(u8, seq2.dsq, seqs[1].dsq);

    try std.testing.expectEqualStrings("longname", seqs[2].name);
    try std.testing.expectEqualSlices(u8, seq3.dsq, seqs[2].dsq);
}

test "header: num_sequences, total_residues, alphabet_type" {
    const allocator = std.testing.allocator;

    var seq1 = try Sequence.fromText(allocator, &alphabet_mod.amino, "prot1", "ACDEF");
    defer seq1.deinit();
    var seq2 = try Sequence.fromText(allocator, &alphabet_mod.amino, "prot2", "GHIKLM");
    defer seq2.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    const seqs = [_]Sequence{ seq1, seq2 };
    try write(buf.writer(allocator).any(), &alphabet_mod.amino, &seqs);

    var fbs = std.io.fixedBufferStream(buf.items);
    const hdr = try readHeader(fbs.reader().any());

    try std.testing.expectEqual(alphabet_mod.AlphabetType.amino, hdr.alphabet_type);
    try std.testing.expectEqual(@as(u64, 2), hdr.num_sequences);
    try std.testing.expectEqual(@as(u64, 11), hdr.total_residues); // 5 + 6
}

test "empty file: write header only, readNext returns null" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), &alphabet_mod.dna, &[_]Sequence{});

    var fbs = std.io.fixedBufferStream(buf.items);
    const hdr = try readHeader(fbs.reader().any());

    try std.testing.expectEqual(@as(u64, 0), hdr.num_sequences);
    try std.testing.expectEqual(@as(u64, 0), hdr.total_residues);

    const next = try readNext(allocator, fbs.reader().any(), &alphabet_mod.dna);
    try std.testing.expect(next == null);
}

test "invalid magic returns error.InvalidFormat" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    // Write a valid file first then corrupt the magic
    try write(buf.writer(allocator).any(), &alphabet_mod.dna, &[_]Sequence{});
    buf.items[0] = 'X'; // corrupt first byte of magic

    var fbs = std.io.fixedBufferStream(buf.items);
    try std.testing.expectError(error.InvalidFormat, readHeader(fbs.reader().any()));
}

test "round trip with amino acid sequences" {
    const allocator = std.testing.allocator;

    var seq = try Sequence.fromText(allocator, &alphabet_mod.amino, "myprot", "ACDEFGHIKLM");
    defer seq.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    const seqs = [_]Sequence{seq};
    try write(buf.writer(allocator).any(), &alphabet_mod.amino, &seqs);

    var fbs = std.io.fixedBufferStream(buf.items);
    _ = try readHeader(fbs.reader().any());
    const result = try readAll(allocator, fbs.reader().any(), &alphabet_mod.amino);
    defer {
        for (result) |*s| @constCast(s).deinit();
        allocator.free(result);
    }

    try std.testing.expectEqual(@as(usize, 1), result.len);
    try std.testing.expectEqualStrings("myprot", result[0].name);
    try std.testing.expectEqualSlices(u8, seq.dsq, result[0].dsq);
}

test "PrefetchReader: threaded read of 5 sequences in blocks of 2" {
    const allocator = std.testing.allocator;

    // Write 5 sequences
    var seqs_arr: [5]Sequence = undefined;
    var to_free: usize = 0;
    defer for (0..to_free) |i| seqs_arr[i].deinit();

    for (0..5) |i| {
        var name_buf: [8]u8 = undefined;
        const name = std.fmt.bufPrint(&name_buf, "s{d}", .{i}) catch unreachable;
        seqs_arr[i] = try Sequence.fromText(allocator, &alphabet_mod.dna, name, "ACGTACGT");
        to_free = i + 1;
    }

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);
    try write(buf.writer(allocator).any(), &alphabet_mod.dna, &seqs_arr);

    // Read with prefetch (heap-allocated PrefetchReader)
    var fbs = std.io.fixedBufferStream(buf.items);
    _ = try readHeader(fbs.reader().any());

    const reader = try PrefetchReader.create(allocator, fbs.reader().any(), &alphabet_mod.dna, 2);
    defer reader.destroy();

    var total_seqs: usize = 0;
    while (reader.read()) |block| {
        defer reader.recycle(block);
        total_seqs += block.count;
    }

    try std.testing.expectEqual(@as(usize, 5), total_seqs);
}

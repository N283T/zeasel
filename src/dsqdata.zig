// dsqdata — binary sequence database format with dual-format support.
//
// DsqData is a tagged union that dispatches to either the native zeasel
// format (single-file ZSQD) or the Easel format (4-file .dsqi/.dsqm/.dsqs).
// The open() constructor auto-detects the format. Writing always produces
// zeasel ZSQD format.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("alphabet.zig").Alphabet;
const alphabet_mod = @import("alphabet.zig");
const Sequence = @import("sequence.zig").Sequence;

// --- Sub-module re-exports ---

pub const zeasel_dsqdata = @import("dsqdata/zeasel.zig");
pub const easel_dsqdata = @import("dsqdata/easel.zig");
pub const pack = @import("dsqdata/pack.zig");

// --- Backwards-compatible API (delegates to zeasel format) ---

pub const write = zeasel_dsqdata.write;
pub const readHeader = zeasel_dsqdata.readHeader;
pub const readNext = zeasel_dsqdata.readNext;
pub const readAll = zeasel_dsqdata.readAll;
pub const PrefetchReader = zeasel_dsqdata.PrefetchReader;
pub const MAGIC = zeasel_dsqdata.MAGIC;
pub const VERSION = zeasel_dsqdata.VERSION;

// --- Shared Header type ---

pub const Header = struct {
    alphabet_type: alphabet_mod.AlphabetType,
    num_sequences: u64,
    total_residues: u64,
};

// --- ZeaselStream: internal helper for streaming zeasel format through the union ---

const ZeaselStream = struct {
    file: std.fs.File,
    hdr: Header,

    fn open(file: std.fs.File) !ZeaselStream {
        // File is already positioned past the 4-byte magic. Read the rest of
        // the 32-byte header (28 remaining bytes).
        const reader = file.deprecatedReader();

        const version = try reader.readInt(u32, .little);
        if (version != zeasel_dsqdata.VERSION) return error.UnsupportedVersion;

        const abc_byte = try reader.readByte();

        var padding: [3]u8 = undefined;
        _ = try reader.readAll(&padding);

        const num_seqs = try reader.readInt(u64, .little);
        const total_res = try reader.readInt(u64, .little);

        var reserved: [4]u8 = undefined;
        _ = try reader.readAll(&reserved);

        return ZeaselStream{
            .file = file,
            .hdr = Header{
                .alphabet_type = @enumFromInt(abc_byte),
                .num_sequences = num_seqs,
                .total_residues = total_res,
            },
        };
    }

    fn next(self: *ZeaselStream, allocator: Allocator, abc: *const Alphabet) !?Sequence {
        return zeasel_dsqdata.readNext(allocator, self.file.deprecatedReader().any(), abc);
    }

    fn deinit(self: *ZeaselStream) void {
        self.file.close();
    }
};

// --- DsqData tagged union ---

pub const DsqData = union(enum) {
    zeasel: ZeaselStream,
    easel: easel_dsqdata.EaselDsqData,

    /// Open a dsqdata database, auto-detecting format from the file contents.
    ///
    /// Detection logic:
    ///   1. Try opening `path` as a file and read the first 4 bytes.
    ///   2. If magic == "ZSQD" -> zeasel format (single file, keep open for streaming).
    ///   3. If magic == 0xc4d3d1b1 (Easel .dsqi magic) -> easel format; strip ".dsqi"
    ///      suffix if present to get the basename.
    ///   4. If file doesn't exist, try path ++ ".dsqi" -> easel basename.
    ///   5. Otherwise -> error.InvalidFormat.
    pub fn open(allocator: Allocator, path: []const u8) !DsqData {
        // Step 1: Try opening path directly.
        const file = std.fs.cwd().openFile(path, .{}) catch |err| {
            if (err == error.FileNotFound) {
                // Step 4: Try path ++ ".dsqi" to detect easel basename.
                const dsqi_path = try std.fmt.allocPrint(allocator, "{s}.dsqi", .{path});
                defer allocator.free(dsqi_path);
                std.fs.cwd().access(dsqi_path, .{}) catch return error.FileNotFound;
                // .dsqi exists, so path is the easel basename.
                return DsqData{ .easel = try easel_dsqdata.EaselDsqData.open(allocator, path) };
            }
            return err;
        };

        // Step 2/3: Read first 4 bytes to detect format.
        var magic: [4]u8 = undefined;
        const n = file.preadAll(&magic, 0) catch {
            file.close();
            return error.InvalidFormat;
        };
        if (n != 4) {
            // File too small to have a valid magic. Could be a stub file for
            // an Easel database. Try path ++ ".dsqi" as fallback.
            file.close();
            const dsqi_path = try std.fmt.allocPrint(allocator, "{s}.dsqi", .{path});
            defer allocator.free(dsqi_path);
            std.fs.cwd().access(dsqi_path, .{}) catch return error.InvalidFormat;
            return DsqData{ .easel = try easel_dsqdata.EaselDsqData.open(allocator, path) };
        }

        // Check for zeasel ZSQD magic.
        if (std.mem.eql(u8, &magic, &zeasel_dsqdata.MAGIC)) {
            // Seek past the 4-byte magic that we already read via pread.
            file.seekTo(4) catch {
                file.close();
                return error.InvalidFormat;
            };
            return DsqData{ .zeasel = ZeaselStream.open(file) catch |err| {
                file.close();
                return err;
            } };
        }

        // Check for Easel .dsqi magic.
        const magic_u32 = std.mem.readInt(u32, &magic, .little);
        if (magic_u32 == 0xc4d3d1b1) {
            file.close();
            // Derive the basename by stripping ".dsqi" suffix if present.
            const basename = if (std.mem.endsWith(u8, path, ".dsqi"))
                path[0 .. path.len - 5]
            else
                path;
            // Need a mutable copy for EaselDsqData.open since it allocPrints.
            return DsqData{ .easel = try easel_dsqdata.EaselDsqData.open(allocator, basename) };
        }

        file.close();
        return error.InvalidFormat;
    }

    /// Read the next sequence from the database.
    /// Returns null at end of stream. Caller owns the returned Sequence.
    pub fn next(self: *DsqData, allocator: Allocator, abc: *const Alphabet) !?Sequence {
        switch (self.*) {
            .zeasel => |*z| return z.next(allocator, abc),
            .easel => |*e| return e.readNext(allocator, abc),
        }
    }

    /// Return the header metadata for this database.
    pub fn header(self: DsqData) Header {
        switch (self) {
            .zeasel => |z| return z.hdr,
            .easel => |e| return Header{
                .alphabet_type = e.alphabet_type,
                .num_sequences = e.num_sequences,
                .total_residues = e.total_residues,
            },
        }
    }

    /// Seek to and read a specific sequence by index (0-based).
    /// Only supported for Easel format (which has a .dsqi index file).
    /// Returns null if seq_idx is out of range.
    pub fn readSequence(self: *DsqData, allocator: Allocator, abc: *const Alphabet, seq_idx: u64) !?Sequence {
        switch (self.*) {
            .zeasel => return error.RandomAccessNotSupported,
            .easel => |*e| return e.readSequence(allocator, abc, seq_idx),
        }
    }

    /// Read a chunk of sequences starting at seq_idx.
    /// Only supported for Easel format (which has a .dsqi index file).
    /// Returns up to max_seqs sequences. Caller owns the returned slice and each Sequence.
    pub fn readChunk(self: *DsqData, allocator: Allocator, abc: *const Alphabet, seq_idx: u64, max_seqs: usize) ![]Sequence {
        switch (self.*) {
            .zeasel => return error.RandomAccessNotSupported,
            .easel => |*e| return e.readChunk(allocator, abc, seq_idx, max_seqs),
        }
    }

    /// Close all file handles and release resources.
    pub fn deinit(self: *DsqData) void {
        switch (self.*) {
            .zeasel => |*z| z.deinit(),
            .easel => |*e| e.deinit(),
        }
    }
};

// --- Sub-module test references ---

test {
    _ = zeasel_dsqdata;
    _ = easel_dsqdata;
    _ = pack;
}

// --- Hub tests ---

test "DsqData.open: auto-detects zeasel format" {
    const allocator = std.testing.allocator;

    // Create a zeasel ZSQD file via write().
    var seq1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "hub_seq", "ACGT");
    defer seq1.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    const originals = [_]Sequence{seq1};
    try write(buf.writer(allocator).any(), &alphabet_mod.dna, &originals);

    // Write to a temp file.
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(.{ .sub_path = "test.zsqd", .data = buf.items });

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test.zsqd", .{dir_path});
    defer allocator.free(full_path);

    var db = try DsqData.open(allocator, full_path);
    defer db.deinit();

    // Verify it detected zeasel format.
    try std.testing.expect(db == .zeasel);

    // Verify header.
    const hdr = db.header();
    try std.testing.expectEqual(@as(u64, 1), hdr.num_sequences);
    try std.testing.expectEqual(@as(u64, 4), hdr.total_residues);

    // Read back the sequence.
    var seq = (try db.next(allocator, &alphabet_mod.dna)).?;
    defer seq.deinit();
    try std.testing.expectEqualStrings("hub_seq", seq.name);
    try std.testing.expectEqualSlices(u8, seq1.dsq, seq.dsq);

    // No more sequences.
    const end = try db.next(allocator, &alphabet_mod.dna);
    try std.testing.expect(end == null);
}

test "DsqData.open: auto-detects easel format by basename" {
    const allocator = std.testing.allocator;

    // Create a test Easel database.
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // 5-bit amino packet: A(0), C(1), D(2), E(3), F(4) + sentinel
    const pkt: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 0) << 25) |
        (@as(u32, 1) << 20) |
        (@as(u32, 2) << 15) |
        (@as(u32, 3) << 10) |
        (@as(u32, 4) << 5) |
        (@as(u32, 31) << 0);

    const packets = [_]u32{pkt};
    const names = [_][]const u8{"easel_seq"};
    const pkt_seqs = [_][]const u32{&packets};
    const lens = [_]u64{5};

    try easel_dsqdata.writeTestEaselDb(tmp_dir.dir, "mydb", 0xAABBCCDD, 3, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/mydb", .{dir_path});
    defer allocator.free(full_path);

    var db = try DsqData.open(allocator, full_path);
    defer db.deinit();

    // Verify it detected easel format.
    try std.testing.expect(db == .easel);

    // Verify header.
    const hdr = db.header();
    try std.testing.expectEqual(@as(u64, 1), hdr.num_sequences);
    try std.testing.expectEqual(alphabet_mod.AlphabetType.amino, hdr.alphabet_type);

    // Read sequence.
    var seq = (try db.next(allocator, &alphabet_mod.amino)).?;
    defer seq.deinit();
    try std.testing.expectEqualStrings("easel_seq", seq.name);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3, 4 }, seq.dsq);
}

test "DsqData.open: auto-detects easel format from .dsqi path" {
    const allocator = std.testing.allocator;

    // Create a test Easel database.
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const pkt: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 0) << 25) |
        (@as(u32, 1) << 20) |
        (@as(u32, 2) << 15) |
        (@as(u32, 3) << 10) |
        (@as(u32, 4) << 5) |
        (@as(u32, 31) << 0);

    const packets = [_]u32{pkt};
    const names = [_][]const u8{"dsqi_seq"};
    const pkt_seqs = [_][]const u32{&packets};
    const lens = [_]u64{5};

    try easel_dsqdata.writeTestEaselDb(tmp_dir.dir, "pathtest", 0x11223344, 3, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp_dir.dir.realpath(".", &path_buf);
    // Pass the .dsqi path explicitly.
    const full_path = try std.fmt.allocPrint(allocator, "{s}/pathtest.dsqi", .{dir_path});
    defer allocator.free(full_path);

    var db = try DsqData.open(allocator, full_path);
    defer db.deinit();

    // Verify it detected easel format.
    try std.testing.expect(db == .easel);

    // Read sequence to verify it works end-to-end.
    var seq = (try db.next(allocator, &alphabet_mod.amino)).?;
    defer seq.deinit();
    try std.testing.expectEqualStrings("dsqi_seq", seq.name);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3, 4 }, seq.dsq);
}

test "DsqData.readSequence: dispatches to easel format" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // Create 2-sequence amino database
    const pkt0: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 0) << 25) |
        (@as(u32, 1) << 20) |
        (@as(u32, 31) << 15) |
        (@as(u32, 31) << 10) |
        (@as(u32, 31) << 5) |
        (@as(u32, 31) << 0);

    const pkt1: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 2) << 25) |
        (@as(u32, 3) << 20) |
        (@as(u32, 4) << 15) |
        (@as(u32, 31) << 10) |
        (@as(u32, 31) << 5) |
        (@as(u32, 31) << 0);

    const packets0 = [_]u32{pkt0};
    const packets1 = [_]u32{pkt1};
    const names = [_][]const u8{ "first", "second" };
    const pkt_seqs = [_][]const u32{ &packets0, &packets1 };
    const lens = [_]u64{ 2, 3 };

    try easel_dsqdata.writeTestEaselDb(tmp_dir.dir, "hubrand", 0x12345678, 3, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/hubrand", .{dir_path});
    defer allocator.free(full_path);

    var db = try DsqData.open(allocator, full_path);
    defer db.deinit();

    try std.testing.expect(db == .easel);

    // Random access: read sequence 1 first
    {
        var seq = (try db.readSequence(allocator, &alphabet_mod.amino, 1)).?;
        defer seq.deinit();
        try std.testing.expectEqualStrings("second", seq.name);
        try std.testing.expectEqualSlices(u8, &[_]u8{ 2, 3, 4 }, seq.dsq);
    }

    // Random access: read sequence 0
    {
        var seq = (try db.readSequence(allocator, &alphabet_mod.amino, 0)).?;
        defer seq.deinit();
        try std.testing.expectEqualStrings("first", seq.name);
        try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1 }, seq.dsq);
    }

    // Out of range
    const out = try db.readSequence(allocator, &alphabet_mod.amino, 99);
    try std.testing.expect(out == null);
}

test "DsqData.readSequence: returns error for zeasel format" {
    const allocator = std.testing.allocator;

    // Create a zeasel ZSQD file.
    var seq1 = try Sequence.fromText(allocator, &alphabet_mod.dna, "test", "ACGT");
    defer seq1.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    const originals = [_]Sequence{seq1};
    try write(buf.writer(allocator).any(), &alphabet_mod.dna, &originals);

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(.{ .sub_path = "test.zsqd", .data = buf.items });

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test.zsqd", .{dir_path});
    defer allocator.free(full_path);

    var db = try DsqData.open(allocator, full_path);
    defer db.deinit();

    try std.testing.expect(db == .zeasel);

    // readSequence should return error for zeasel format
    try std.testing.expectError(error.RandomAccessNotSupported, db.readSequence(allocator, &alphabet_mod.dna, 0));
}

test "DsqData.readChunk: dispatches to easel format" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const make_pkt = struct {
        fn f(code: u32) u32 {
            return pack.EOD | pack.FIVEBIT |
                (code << 25) |
                (@as(u32, 31) << 20) |
                (@as(u32, 31) << 15) |
                (@as(u32, 31) << 10) |
                (@as(u32, 31) << 5) |
                (@as(u32, 31) << 0);
        }
    }.f;

    const pkts = [3][1]u32{
        .{make_pkt(0)},
        .{make_pkt(1)},
        .{make_pkt(2)},
    };

    const names = [_][]const u8{ "a", "b", "c" };
    const pkt_seqs = [_][]const u32{ &pkts[0], &pkts[1], &pkts[2] };
    const lens = [_]u64{ 1, 1, 1 };

    try easel_dsqdata.writeTestEaselDb(tmp_dir.dir, "hubchunk", 0xABCD1234, 3, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp_dir.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/hubchunk", .{dir_path});
    defer allocator.free(full_path);

    var db = try DsqData.open(allocator, full_path);
    defer db.deinit();

    const chunk = try db.readChunk(allocator, &alphabet_mod.amino, 1, 2);
    defer {
        for (chunk) |*s| {
            var seq = s.*;
            seq.deinit();
        }
        allocator.free(chunk);
    }

    try std.testing.expectEqual(@as(usize, 2), chunk.len);
    try std.testing.expectEqualStrings("b", chunk[0].name);
    try std.testing.expectEqualStrings("c", chunk[1].name);
}

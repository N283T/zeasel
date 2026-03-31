// SSI (Simple Sequence Index) — module hub for sequence index formats.
//
// SsiIndex is a tagged union that dispatches to either the native zeasel
// format (in-memory hash map) or the Easel v3.0 format (disk-based binary
// search). The open() constructor auto-detects the format from magic bytes.

const std = @import("std");
const Allocator = std.mem.Allocator;

const zeasel_ssi = @import("ssi/zeasel.zig");
const easel_ssi = @import("ssi/easel.zig");

pub const ZeaselIndex = zeasel_ssi.ZeaselIndex;
pub const EaselIndex = easel_ssi.EaselIndex;
pub const FileInfo = easel_ssi.FileInfo;

pub const SsiEntry = struct {
    name: []const u8,
    offset: u64,
    data_offset: u64,
    seq_len: u64,
    file_id: u16 = 0, // Index into file_names array for multi-file support
};

pub const SsiIndex = union(enum) {
    zeasel: ZeaselIndex,
    easel: EaselIndex,

    /// Open an SSI index file, auto-detecting format by magic bytes.
    /// For Easel: file handle stays open for disk-based lookup.
    /// For zeasel: file is read into memory and closed.
    pub fn open(allocator: Allocator, path: []const u8) !SsiIndex {
        const file = try std.fs.cwd().openFile(path, .{});

        // Read first 4 bytes to detect format.
        var magic: [4]u8 = undefined;
        const n = file.preadAll(&magic, 0) catch {
            file.close();
            return error.InvalidFormat;
        };
        if (n != 4) {
            file.close();
            return error.InvalidFormat;
        }

        const magic_u32 = std.mem.readInt(u32, &magic, .big);

        if (magic_u32 == easel_ssi.EASEL_MAGIC or magic_u32 == easel_ssi.EASEL_MAGIC_SWAPPED) {
            // Easel format — EaselIndex.read keeps the file handle open.
            // On error, close file since EaselIndex won't own it.
            return SsiIndex{ .easel = EaselIndex.read(allocator, file) catch |err| {
                file.close();
                return err;
            } };
        } else if (std.mem.eql(u8, &magic, "ZSSI")) {
            // Zeasel format — read entire file into memory, close, parse.
            // preadAll does not change the file position, so readToEndAlloc
            // reads from position 0 (the file was just opened).
            const data = file.readToEndAlloc(allocator, std.math.maxInt(usize)) catch |err| {
                file.close();
                return err;
            };
            defer allocator.free(data);
            file.close();

            var stream = std.io.fixedBufferStream(data);
            const zidx = try ZeaselIndex.read(allocator, stream.reader().any());
            return SsiIndex{ .zeasel = zidx };
        } else {
            file.close();
            return error.InvalidFormat;
        }
    }

    /// Look up a sequence by name.
    /// Caller owns the returned SsiEntry.name and must free it with `allocator`.
    pub fn lookup(self: *SsiIndex, allocator: Allocator, key: []const u8) !?SsiEntry {
        switch (self.*) {
            .zeasel => |zidx| {
                const raw = zidx.lookup(key) orelse return null;
                return SsiEntry{
                    .name = try allocator.dupe(u8, raw.name),
                    .offset = raw.offset,
                    .data_offset = raw.data_offset,
                    .seq_len = raw.seq_len,
                    .file_id = raw.file_id,
                };
            },
            .easel => |*eidx| return try eidx.lookup(allocator, key),
        }
    }

    /// Look up by ordinal position (0-indexed).
    /// Caller owns the returned SsiEntry.name and must free it with `allocator`.
    pub fn findNumber(self: *SsiIndex, allocator: Allocator, number: usize) !?SsiEntry {
        switch (self.*) {
            .zeasel => |zidx| {
                const raw = zidx.findNumber(number) orelse return null;
                return SsiEntry{
                    .name = try allocator.dupe(u8, raw.name),
                    .offset = raw.offset,
                    .data_offset = raw.data_offset,
                    .seq_len = raw.seq_len,
                    .file_id = raw.file_id,
                };
            },
            .easel => |*eidx| return try eidx.findNumber(allocator, @intCast(number)),
        }
    }

    /// Compute subsequence positioning info.
    /// Only supported for zeasel variant; returns null for easel.
    pub fn findSubseq(self: SsiIndex, name: []const u8, start: u64, end: u64) ?struct { entry: SsiEntry, start: u64, end: u64 } {
        switch (self) {
            .zeasel => |zidx| return zidx.findSubseq(name, start, end),
            .easel => return null,
        }
    }

    /// Get the file name for a given file_id.
    pub fn fileInfo(self: SsiIndex, file_id: u16) ?[]const u8 {
        switch (self) {
            .zeasel => |zidx| return zidx.fileInfo(file_id),
            .easel => |eidx| {
                if (file_id >= eidx.files.len) return null;
                return eidx.files[file_id].name;
            },
        }
    }

    /// Write the index to a binary stream.
    /// Only supported for zeasel variant.
    pub fn write(self: SsiIndex, dest: std.io.AnyWriter) !void {
        switch (self) {
            .zeasel => |zidx| return zidx.write(dest),
            .easel => return error.ReadOnly,
        }
    }

    pub fn deinit(self: *SsiIndex) void {
        switch (self.*) {
            .zeasel => |*zidx| zidx.deinit(),
            .easel => |*eidx| eidx.deinit(),
        }
    }
};

test {
    _ = zeasel_ssi;
    _ = easel_ssi;
}

test "SsiIndex.open: auto-detects zeasel format" {
    const allocator = std.testing.allocator;

    const fasta_data = ">seq1\nACGT\n";
    var zidx = try ZeaselIndex.buildFromFasta(allocator, fasta_data);
    defer zidx.deinit();

    var ssi_buf = std.ArrayList(u8){};
    defer ssi_buf.deinit(allocator);
    try zidx.write(ssi_buf.writer(allocator).any());

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(.{ .sub_path = "test.ssi", .data = ssi_buf.items });

    const dir_path = try tmp_dir.dir.realpathAlloc(allocator, ".");
    defer allocator.free(dir_path);
    const ssi_path = try std.fmt.allocPrint(allocator, "{s}/test.ssi", .{dir_path});
    defer allocator.free(ssi_path);

    var idx = try SsiIndex.open(allocator, ssi_path);
    defer idx.deinit();

    try std.testing.expect(idx == .zeasel);

    // SsiIndex.lookup always returns an owned name — caller must free.
    const e = (try idx.lookup(allocator, "seq1")) orelse return error.ExpectedEntry;
    defer allocator.free(e.name);
    try std.testing.expectEqual(@as(u64, 4), e.seq_len);
}

test "SsiIndex.open: auto-detects easel format" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try easel_ssi.writeEaselHeader(&buf, allocator, .{ .nprimary = 1 });
    try easel_ssi.writePrimaryKey(&buf, allocator, .{ .key = "test_seq", .r_off = 0, .d_off = 10, .len = 42 });

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(.{ .sub_path = "test.ssi", .data = buf.items });

    const dir_path = try tmp_dir.dir.realpathAlloc(allocator, ".");
    defer allocator.free(dir_path);
    const ssi_path = try std.fmt.allocPrint(allocator, "{s}/test.ssi", .{dir_path});
    defer allocator.free(ssi_path);

    var idx = try SsiIndex.open(allocator, ssi_path);
    defer idx.deinit();

    try std.testing.expect(idx == .easel);

    // EaselIndex.lookup returns owned name, must free.
    const e = (try idx.lookup(allocator, "test_seq")) orelse return error.ExpectedEntry;
    defer allocator.free(e.name);
    try std.testing.expectEqual(@as(u64, 42), e.seq_len);
}

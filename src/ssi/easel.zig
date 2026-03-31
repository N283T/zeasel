// EaselIndex — parser for Easel SSI v3.0 binary index files.
//
// The Easel SSI format is big-endian by default but may be byteswapped
// if the index was created on a little-endian machine. The magic value
// determines endianness: 0xd3d3c9b3 (native/big-endian) or 0xb3c9d3d3
// (byteswapped/little-endian).

const std = @import("std");
const Allocator = std.mem.Allocator;

const ssi = @import("../ssi.zig");
const SsiEntry = ssi.SsiEntry;

pub const EASEL_MAGIC: u32 = 0xd3d3c9b3;
pub const EASEL_MAGIC_SWAPPED: u32 = 0xb3c9d3d3;

pub const FileInfo = struct {
    name: []const u8,
    format: u32,
    flags: u32,
    bpl: u32, // bytes per line
    rpl: u32, // residues per line
};

pub const EaselIndex = struct {
    file: std.fs.File,
    byteswap: bool,
    offsz: u8,
    nfiles: u16,
    nprimary: u64,
    nsecondary: u64,
    plen: u32,
    slen: u32,
    precsize: u32,
    srecsize: u32,
    poffset: u64,
    soffset: u64,
    files: []FileInfo,
    allocator: Allocator,

    // Maximum header size: fixed (54 bytes) + 3 * 8-byte offsets = 78 bytes.
    const max_header_size = 78;

    /// Parse an Easel SSI v3.0 index from an open file handle.
    /// The file handle is kept open for later disk-based lookups.
    /// Caller owns the returned EaselIndex and must call deinit().
    pub fn read(allocator: Allocator, file: std.fs.File) !EaselIndex {
        // Read the maximum possible header into a stack buffer.
        var hdr_buf: [max_header_size]u8 = undefined;
        const n_read = try file.preadAll(&hdr_buf, 0);
        if (n_read < 14) return error.InvalidFormat; // Need at least magic + flags + offsz + nfiles

        var pos: usize = 0;

        // Read magic (always big-endian on disk).
        const raw_magic = readIntFromBuf(u32, hdr_buf[pos..], .big);
        pos += 4;

        const byteswap = if (raw_magic == EASEL_MAGIC)
            false
        else if (raw_magic == EASEL_MAGIC_SWAPPED)
            true
        else
            return error.InvalidFormat;

        const endian: std.builtin.Endian = if (byteswap) .little else .big;

        // flags (skip)
        pos += 4;

        const offsz_raw = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        if (offsz_raw != 4 and offsz_raw != 8) return error.InvalidFormat;
        const offsz: u8 = @intCast(offsz_raw);

        // Verify we read enough header bytes for the variable-size offsets.
        const actual_header_size: usize = 54 + 3 * @as(usize, offsz);
        if (n_read < actual_header_size) return error.InvalidFormat;

        const nfiles = readIntFromBuf(u16, hdr_buf[pos..], endian);
        pos += 2;
        const nprimary = readIntFromBuf(u64, hdr_buf[pos..], endian);
        pos += 8;
        const nsecondary = readIntFromBuf(u64, hdr_buf[pos..], endian);
        pos += 8;
        const flen = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const plen = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const slen = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const frecsize = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const precsize = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const srecsize = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;

        const foffset = readOffsetFromBuf(hdr_buf[pos..], offsz, endian);
        pos += offsz;
        const poffset = readOffsetFromBuf(hdr_buf[pos..], offsz, endian);
        pos += offsz;
        const soffset = readOffsetFromBuf(hdr_buf[pos..], offsz, endian);

        // Read file info section from disk.
        const file_section_size = @as(usize, frecsize) * @as(usize, nfiles);
        const file_section_buf = try allocator.alloc(u8, file_section_size);
        defer allocator.free(file_section_buf);

        const file_n_read = try file.preadAll(file_section_buf, foffset);
        if (file_n_read < file_section_size) return error.InvalidFormat;

        // Parse file info entries.
        const files = try allocator.alloc(FileInfo, nfiles);
        var files_done: usize = 0;
        errdefer {
            for (files[0..files_done]) |f| allocator.free(f.name);
            allocator.free(files);
        }

        var fpos: usize = 0;
        for (0..nfiles) |i| {
            const rec_start = fpos;
            const name_bytes = file_section_buf[fpos..][0..flen];
            fpos += flen;

            // Find actual string length (up to first null).
            var name_len: usize = 0;
            for (name_bytes) |c| {
                if (c == 0) break;
                name_len += 1;
            }
            const name = try allocator.dupe(u8, name_bytes[0..name_len]);
            errdefer allocator.free(name);

            const format = readIntFromBuf(u32, file_section_buf[fpos..], endian);
            fpos += 4;
            const file_flags = readIntFromBuf(u32, file_section_buf[fpos..], endian);
            fpos += 4;
            const bpl = readIntFromBuf(u32, file_section_buf[fpos..], endian);
            fpos += 4;
            const rpl = readIntFromBuf(u32, file_section_buf[fpos..], endian);
            fpos += 4;

            // Advance to next record using frecsize, which may be larger than
            // flen + 16 if a future Easel version adds extra fields.
            fpos = rec_start + @as(usize, frecsize);

            files[i] = FileInfo{
                .name = name,
                .format = format,
                .flags = file_flags,
                .bpl = bpl,
                .rpl = rpl,
            };
            files_done += 1;
        }

        return EaselIndex{
            .file = file,
            .byteswap = byteswap,
            .offsz = offsz,
            .nfiles = nfiles,
            .nprimary = nprimary,
            .nsecondary = nsecondary,
            .plen = plen,
            .slen = slen,
            .precsize = precsize,
            .srecsize = srecsize,
            .poffset = poffset,
            .soffset = soffset,
            .files = files,
            .allocator = allocator,
        };
    }

    /// Look up a key in the index, trying primary keys first, then secondary.
    /// Returns an SsiEntry on success (caller owns entry.name) or null if not found.
    pub fn lookup(self: *EaselIndex, allocator: Allocator, key: []const u8) !?SsiEntry {
        // Try primary keys first.
        if (try self.binarySearchPrimary(allocator, key)) |entry| {
            return entry;
        }

        // Try secondary keys if present.
        if (self.nsecondary > 0) {
            if (try self.binarySearchSecondary(allocator, key)) |primary_key| {
                defer allocator.free(primary_key);

                // Look up the primary key to get the full entry.
                if (try self.binarySearchPrimary(allocator, primary_key)) |entry| {
                    // Replace the primary key name with the search key
                    // so the caller sees what they asked for.
                    allocator.free(entry.name);
                    const name_copy = try allocator.dupe(u8, key);
                    return SsiEntry{
                        .name = name_copy,
                        .offset = entry.offset,
                        .data_offset = entry.data_offset,
                        .seq_len = entry.seq_len,
                        .file_id = entry.file_id,
                    };
                }
            }
        }

        return null;
    }

    /// Look up by ordinal position (0-indexed).
    /// Caller owns the returned SsiEntry.name and must free it.
    pub fn findNumber(self: *EaselIndex, allocator: Allocator, number: u64) !?SsiEntry {
        if (number >= self.nprimary) return null;

        const endian: std.builtin.Endian = if (self.byteswap) .little else .big;

        // Read the full primary key record.
        const rec_buf = try allocator.alloc(u8, self.precsize);
        defer allocator.free(rec_buf);

        const offset = self.poffset + @as(u64, self.precsize) * number;
        const n = try self.file.preadAll(rec_buf, offset);
        if (n != self.precsize) return error.InvalidFormat;

        // Parse key name (null-padded C string).
        var name_len: usize = 0;
        for (rec_buf[0..self.plen]) |c| {
            if (c == 0) break;
            name_len += 1;
        }
        const name = try allocator.dupe(u8, rec_buf[0..name_len]);
        errdefer allocator.free(name);

        // Parse remaining fields after the key.
        var pos: usize = self.plen;
        const fnum = readIntFromBuf(u16, rec_buf[pos..], endian);
        pos += 2;
        const r_off = readOffsetFromBuf(rec_buf[pos..], self.offsz, endian);
        pos += self.offsz;
        const d_off = readOffsetFromBuf(rec_buf[pos..], self.offsz, endian);
        pos += self.offsz;
        const len_raw = readIntFromBuf(i64, rec_buf[pos..], endian);

        return SsiEntry{
            .name = name,
            .offset = r_off,
            .data_offset = d_off,
            .seq_len = if (len_raw >= 0) @intCast(len_raw) else 0,
            .file_id = fnum,
        };
    }

    /// Binary search the primary key section for the given key.
    /// Returns an SsiEntry (caller owns entry.name) or null.
    fn binarySearchPrimary(self: *EaselIndex, allocator: Allocator, key: []const u8) !?SsiEntry {
        if (self.nprimary == 0) return null;

        const endian: std.builtin.Endian = if (self.byteswap) .little else .big;
        const rec_buf = try allocator.alloc(u8, self.precsize);
        defer allocator.free(rec_buf);

        var low: u64 = 0;
        var high: u64 = self.nprimary;

        while (low < high) {
            const mid = low + (high - low) / 2;
            const offset = self.poffset + @as(u64, self.precsize) * mid;

            const n_read = try self.file.preadAll(rec_buf, offset);
            if (n_read < self.precsize) return error.InvalidFormat;

            // Extract the null-padded key.
            const rec_key_bytes = rec_buf[0..self.plen];
            var rec_key_len: usize = 0;
            for (rec_key_bytes) |c| {
                if (c == 0) break;
                rec_key_len += 1;
            }
            const rec_key = rec_key_bytes[0..rec_key_len];

            const order = std.mem.order(u8, rec_key, key);
            switch (order) {
                .eq => {
                    // Parse the rest of the record after the key.
                    var pos: usize = self.plen;
                    const fnum = readIntFromBuf(u16, rec_buf[pos..], endian);
                    pos += 2;
                    const r_off = readOffsetFromBuf(rec_buf[pos..], self.offsz, endian);
                    pos += self.offsz;
                    const d_off = readOffsetFromBuf(rec_buf[pos..], self.offsz, endian);
                    pos += self.offsz;
                    const len_raw = readIntFromBuf(i64, rec_buf[pos..], endian);
                    const seq_len: u64 = if (len_raw < 0) 0 else @intCast(len_raw);

                    const name_copy = try allocator.dupe(u8, key);
                    return SsiEntry{
                        .name = name_copy,
                        .offset = r_off,
                        .data_offset = d_off,
                        .seq_len = seq_len,
                        .file_id = fnum,
                    };
                },
                .lt => low = mid + 1,
                .gt => high = mid,
            }
        }

        return null;
    }

    /// Binary search the secondary key section for the given key.
    /// Returns the associated primary key as an owned slice, or null.
    fn binarySearchSecondary(self: *EaselIndex, allocator: Allocator, key: []const u8) !?[]const u8 {
        if (self.nsecondary == 0) return null;

        const rec_buf = try allocator.alloc(u8, self.srecsize);
        defer allocator.free(rec_buf);

        var low: u64 = 0;
        var high: u64 = self.nsecondary;

        while (low < high) {
            const mid = low + (high - low) / 2;
            const offset = self.soffset + @as(u64, self.srecsize) * mid;

            const n_read = try self.file.preadAll(rec_buf, offset);
            if (n_read < self.srecsize) return error.InvalidFormat;

            // Extract the null-padded secondary key.
            const skey_bytes = rec_buf[0..self.slen];
            var skey_len: usize = 0;
            for (skey_bytes) |c| {
                if (c == 0) break;
                skey_len += 1;
            }
            const skey = skey_bytes[0..skey_len];

            const order = std.mem.order(u8, skey, key);
            switch (order) {
                .eq => {
                    // Extract the null-padded primary key.
                    const pkey_bytes = rec_buf[self.slen..][0..self.plen];
                    var pkey_len: usize = 0;
                    for (pkey_bytes) |c| {
                        if (c == 0) break;
                        pkey_len += 1;
                    }
                    return try allocator.dupe(u8, pkey_bytes[0..pkey_len]);
                },
                .lt => low = mid + 1,
                .gt => high = mid,
            }
        }

        return null;
    }

    /// Free all allocated memory and close the file handle.
    pub fn deinit(self: *EaselIndex) void {
        for (self.files) |f| self.allocator.free(f.name);
        self.allocator.free(self.files);
        self.file.close();
    }
};

/// Read an integer of type T from a byte slice at the given endianness.
fn readIntFromBuf(comptime T: type, buf: []const u8, endian: std.builtin.Endian) T {
    const n = @divExact(@typeInfo(T).int.bits, 8);
    return std.mem.readInt(T, buf[0..n], endian);
}

/// Read a 4-byte or 8-byte offset from a byte slice.
fn readOffsetFromBuf(buf: []const u8, offsz: u8, endian: std.builtin.Endian) u64 {
    return switch (offsz) {
        4 => @as(u64, readIntFromBuf(u32, buf, endian)),
        8 => readIntFromBuf(u64, buf, endian),
        else => unreachable,
    };
}

/// Write a valid Easel SSI v3.0 header and file info section for testing.
/// Calculates section offsets automatically based on the provided parameters.
pub fn writeEaselHeader(
    buf: *std.ArrayList(u8),
    allocator: Allocator,
    opts: struct {
        nfiles: u16 = 1,
        nprimary: u64 = 0,
        nsecondary: u64 = 0,
        flen: u32 = 32,
        plen: u32 = 32,
        slen: u32 = 32,
        offsz: u8 = 8,
        byteswap: bool = false,
        file_names: []const []const u8 = &.{"testfile.fa"},
        file_formats: []const u32 = &.{0},
        file_flags: []const u32 = &.{0},
        file_bpls: []const u32 = &.{80},
        file_rpls: []const u32 = &.{60},
    },
) !void {
    const writer = buf.writer(allocator);

    const endian: std.builtin.Endian = if (opts.byteswap) .little else .big;
    const offsz: u8 = opts.offsz;

    // Calculate sizes.
    const frecsize: u32 = opts.flen + 16; // flen + 4*u32
    const precsize: u32 = opts.plen + 2 + 2 * @as(u32, offsz) + 8;
    const srecsize: u32 = opts.slen + opts.plen;

    // Header size: fixed fields + 3 offsets.
    const header_size: u64 = 4 + 4 + 4 + 2 + 8 + 8 + 4 + 4 + 4 + 4 + 4 + 4 + 3 * @as(u64, offsz);
    const foffset: u64 = header_size;
    const poffset: u64 = foffset + @as(u64, frecsize) * @as(u64, opts.nfiles);
    const soffset: u64 = poffset + @as(u64, precsize) * opts.nprimary;

    // Write magic.
    const magic: u32 = if (opts.byteswap) EASEL_MAGIC_SWAPPED else EASEL_MAGIC;
    try writer.writeInt(u32, magic, .big);

    // flags
    try writer.writeInt(u32, 0, endian);
    // offsz
    try writer.writeInt(u32, @as(u32, offsz), endian);
    // nfiles
    try writer.writeInt(u16, opts.nfiles, endian);
    // nprimary
    try writer.writeInt(u64, opts.nprimary, endian);
    // nsecondary
    try writer.writeInt(u64, opts.nsecondary, endian);
    // flen, plen, slen
    try writer.writeInt(u32, opts.flen, endian);
    try writer.writeInt(u32, opts.plen, endian);
    try writer.writeInt(u32, opts.slen, endian);
    // frecsize, precsize, srecsize
    try writer.writeInt(u32, frecsize, endian);
    try writer.writeInt(u32, precsize, endian);
    try writer.writeInt(u32, srecsize, endian);

    // offsets (foffset, poffset, soffset)
    try writeOffset(writer, foffset, offsz, endian);
    try writeOffset(writer, poffset, offsz, endian);
    try writeOffset(writer, soffset, offsz, endian);

    // File info section.
    for (0..opts.nfiles) |i| {
        const name = if (i < opts.file_names.len) opts.file_names[i] else "unknown";
        const fmt = if (i < opts.file_formats.len) opts.file_formats[i] else 0;
        const flg = if (i < opts.file_flags.len) opts.file_flags[i] else 0;
        const bpl = if (i < opts.file_bpls.len) opts.file_bpls[i] else 80;
        const rpl = if (i < opts.file_rpls.len) opts.file_rpls[i] else 60;

        // Write null-padded filename.
        try writer.writeAll(name);
        const pad_len = opts.flen - @as(u32, @intCast(name.len));
        for (0..pad_len) |_| try writer.writeByte(0);

        try writer.writeInt(u32, fmt, endian);
        try writer.writeInt(u32, flg, endian);
        try writer.writeInt(u32, bpl, endian);
        try writer.writeInt(u32, rpl, endian);
    }
}

pub fn writeOffset(writer: anytype, value: u64, offsz: u8, endian: std.builtin.Endian) !void {
    switch (offsz) {
        4 => try writer.writeInt(u32, @intCast(value), endian),
        8 => try writer.writeInt(u64, value, endian),
        else => return error.InvalidFormat,
    }
}

// --- Tests ---

test "EaselIndex.read: parses header and file info" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{
        .nfiles = 1,
        .nprimary = 3,
        .nsecondary = 0,
        .flen = 32,
        .plen = 32,
        .slen = 32,
        .offsz = 8,
        .file_names = &.{"sequences.fasta"},
        .file_bpls = &.{80},
        .file_rpls = &.{60},
    });

    // Write to a temp file since EaselIndex needs a seekable file handle.
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(.{ .sub_path = "test.ssi", .data = buf.items });
    const file = try tmp_dir.dir.openFile("test.ssi", .{});
    // File will be closed by deinit.

    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    try std.testing.expectEqual(false, idx.byteswap);
    try std.testing.expectEqual(@as(u8, 8), idx.offsz);
    try std.testing.expectEqual(@as(u16, 1), idx.nfiles);
    try std.testing.expectEqual(@as(u64, 3), idx.nprimary);
    try std.testing.expectEqual(@as(u64, 0), idx.nsecondary);
    try std.testing.expectEqual(@as(u32, 32), idx.plen);
    try std.testing.expectEqual(@as(u32, 32), idx.slen);

    // File info.
    try std.testing.expectEqual(@as(usize, 1), idx.files.len);
    try std.testing.expectEqualStrings("sequences.fasta", idx.files[0].name);
    try std.testing.expectEqual(@as(u32, 80), idx.files[0].bpl);
    try std.testing.expectEqual(@as(u32, 60), idx.files[0].rpl);
}

test "EaselIndex.read: rejects invalid magic" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // Write "ZSSI" followed by padding — not a valid Easel magic.
    var bad_data: [64]u8 = undefined;
    @memset(&bad_data, 0);
    @memcpy(bad_data[0..4], "ZSSI");

    try tmp_dir.dir.writeFile(.{ .sub_path = "bad.ssi", .data = &bad_data });
    const file = try tmp_dir.dir.openFile("bad.ssi", .{});
    defer file.close();

    try std.testing.expectError(error.InvalidFormat, EaselIndex.read(allocator, file));
}

test "EaselIndex.read: parses byteswapped header" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{
        .nfiles = 1,
        .nprimary = 5,
        .nsecondary = 2,
        .byteswap = true,
        .file_names = &.{"swapped.fa"},
        .file_bpls = &.{100},
        .file_rpls = &.{80},
    });

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(.{ .sub_path = "swap.ssi", .data = buf.items });
    const file = try tmp_dir.dir.openFile("swap.ssi", .{});

    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    try std.testing.expectEqual(true, idx.byteswap);
    try std.testing.expectEqual(@as(u64, 5), idx.nprimary);
    try std.testing.expectEqual(@as(u64, 2), idx.nsecondary);
    try std.testing.expectEqualStrings("swapped.fa", idx.files[0].name);
    try std.testing.expectEqual(@as(u32, 100), idx.files[0].bpl);
    try std.testing.expectEqual(@as(u32, 80), idx.files[0].rpl);
}

/// Write a single primary key record to a buffer.
pub fn writePrimaryKey(
    buf: *std.ArrayList(u8),
    allocator: Allocator,
    opts: struct {
        key: []const u8,
        plen: u32 = 32,
        fnum: u16 = 0,
        r_off: u64 = 0,
        d_off: u64 = 0,
        len: i64 = 0,
        offsz: u8 = 8,
        byteswap: bool = false,
    },
) !void {
    const writer = buf.writer(allocator);
    const endian: std.builtin.Endian = if (opts.byteswap) .little else .big;

    // Write null-padded key.
    try writer.writeAll(opts.key);
    const pad_len = opts.plen - @as(u32, @intCast(opts.key.len));
    for (0..pad_len) |_| try writer.writeByte(0);

    // fnum
    try writer.writeInt(u16, opts.fnum, endian);
    // r_off
    try writeOffset(writer, opts.r_off, opts.offsz, endian);
    // d_off
    try writeOffset(writer, opts.d_off, opts.offsz, endian);
    // len
    try writer.writeInt(i64, opts.len, endian);
}

/// Write a single secondary key record to a buffer.
pub fn writeSecondaryKey(
    buf: *std.ArrayList(u8),
    allocator: Allocator,
    opts: struct {
        skey: []const u8,
        pkey: []const u8,
        slen: u32 = 32,
        plen: u32 = 32,
    },
) !void {
    const writer = buf.writer(allocator);

    // Write null-padded secondary key.
    try writer.writeAll(opts.skey);
    const spad = opts.slen - @as(u32, @intCast(opts.skey.len));
    for (0..spad) |_| try writer.writeByte(0);

    // Write null-padded primary key.
    try writer.writeAll(opts.pkey);
    const ppad = opts.plen - @as(u32, @intCast(opts.pkey.len));
    for (0..ppad) |_| try writer.writeByte(0);
}

test "EaselIndex.read: parses header with 4-byte offsets" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{
        .nfiles = 2,
        .nprimary = 10,
        .nsecondary = 3,
        .flen = 16,
        .plen = 16,
        .slen = 16,
        .offsz = 4,
        .file_names = &.{ "small1.fa", "small2.fa" },
        .file_formats = &.{ 1, 2 },
        .file_flags = &.{ 0, 0 },
        .file_bpls = &.{ 60, 70 },
        .file_rpls = &.{ 50, 55 },
    });

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(.{ .sub_path = "off4.ssi", .data = buf.items });
    const file = try tmp_dir.dir.openFile("off4.ssi", .{});

    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    try std.testing.expectEqual(false, idx.byteswap);
    try std.testing.expectEqual(@as(u8, 4), idx.offsz);
    try std.testing.expectEqual(@as(u16, 2), idx.nfiles);
    try std.testing.expectEqual(@as(u64, 10), idx.nprimary);
    try std.testing.expectEqual(@as(u64, 3), idx.nsecondary);
    try std.testing.expectEqual(@as(u32, 16), idx.plen);
    try std.testing.expectEqual(@as(u32, 16), idx.slen);

    // File info.
    try std.testing.expectEqual(@as(usize, 2), idx.files.len);
    try std.testing.expectEqualStrings("small1.fa", idx.files[0].name);
    try std.testing.expectEqual(@as(u32, 1), idx.files[0].format);
    try std.testing.expectEqual(@as(u32, 60), idx.files[0].bpl);
    try std.testing.expectEqual(@as(u32, 50), idx.files[0].rpl);
    try std.testing.expectEqualStrings("small2.fa", idx.files[1].name);
    try std.testing.expectEqual(@as(u32, 2), idx.files[1].format);
    try std.testing.expectEqual(@as(u32, 70), idx.files[1].bpl);
    try std.testing.expectEqual(@as(u32, 55), idx.files[1].rpl);
}

test "EaselIndex.lookup: finds primary key by binary search" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{
        .nfiles = 1,
        .nprimary = 3,
        .nsecondary = 0,
        .plen = 32,
        .slen = 32,
        .offsz = 8,
        .file_names = &.{"test.fa"},
    });

    // Write 3 sorted primary keys: alpha, beta, gamma.
    try writePrimaryKey(&buf, allocator, .{ .key = "alpha", .fnum = 0, .r_off = 100, .d_off = 110, .len = 50 });
    try writePrimaryKey(&buf, allocator, .{ .key = "beta", .fnum = 0, .r_off = 200, .d_off = 220, .len = 75 });
    try writePrimaryKey(&buf, allocator, .{ .key = "gamma", .fnum = 0, .r_off = 300, .d_off = 330, .len = 120 });

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(.{ .sub_path = "prim.ssi", .data = buf.items });
    const file = try tmp_dir.dir.openFile("prim.ssi", .{});

    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    // Find "alpha".
    {
        const entry = (try idx.lookup(allocator, "alpha")).?;
        defer allocator.free(entry.name);
        try std.testing.expectEqualStrings("alpha", entry.name);
        try std.testing.expectEqual(@as(u64, 100), entry.offset);
        try std.testing.expectEqual(@as(u64, 110), entry.data_offset);
        try std.testing.expectEqual(@as(u64, 50), entry.seq_len);
        try std.testing.expectEqual(@as(u16, 0), entry.file_id);
    }

    // Find "beta".
    {
        const entry = (try idx.lookup(allocator, "beta")).?;
        defer allocator.free(entry.name);
        try std.testing.expectEqualStrings("beta", entry.name);
        try std.testing.expectEqual(@as(u64, 200), entry.offset);
        try std.testing.expectEqual(@as(u64, 220), entry.data_offset);
        try std.testing.expectEqual(@as(u64, 75), entry.seq_len);
    }

    // Find "gamma".
    {
        const entry = (try idx.lookup(allocator, "gamma")).?;
        defer allocator.free(entry.name);
        try std.testing.expectEqualStrings("gamma", entry.name);
        try std.testing.expectEqual(@as(u64, 300), entry.offset);
        try std.testing.expectEqual(@as(u64, 330), entry.data_offset);
        try std.testing.expectEqual(@as(u64, 120), entry.seq_len);
    }

    // Missing key returns null.
    {
        const result = try idx.lookup(allocator, "delta");
        try std.testing.expect(result == null);
    }
}

test "EaselIndex.lookup: finds secondary key and resolves to primary" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{
        .nfiles = 1,
        .nprimary = 2,
        .nsecondary = 2,
        .plen = 32,
        .slen = 32,
        .offsz = 8,
        .file_names = &.{"test.fa"},
    });

    // Write 2 sorted primary keys.
    try writePrimaryKey(&buf, allocator, .{ .key = "seq1", .fnum = 0, .r_off = 100, .d_off = 110, .len = 50 });
    try writePrimaryKey(&buf, allocator, .{ .key = "seq2", .fnum = 0, .r_off = 200, .d_off = 220, .len = 75 });

    // Write 2 sorted secondary keys: acc1 -> seq1, acc2 -> seq2.
    try writeSecondaryKey(&buf, allocator, .{ .skey = "acc1", .pkey = "seq1" });
    try writeSecondaryKey(&buf, allocator, .{ .skey = "acc2", .pkey = "seq2" });

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(.{ .sub_path = "sec.ssi", .data = buf.items });
    const file = try tmp_dir.dir.openFile("sec.ssi", .{});

    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    // Look up via secondary key "acc2".
    {
        const entry = (try idx.lookup(allocator, "acc2")).?;
        defer allocator.free(entry.name);
        // The returned name should be the SEARCH key, not the primary key.
        try std.testing.expectEqualStrings("acc2", entry.name);
        try std.testing.expectEqual(@as(u64, 200), entry.offset);
        try std.testing.expectEqual(@as(u64, 220), entry.data_offset);
        try std.testing.expectEqual(@as(u64, 75), entry.seq_len);
    }

    // Look up via secondary key "acc1".
    {
        const entry = (try idx.lookup(allocator, "acc1")).?;
        defer allocator.free(entry.name);
        try std.testing.expectEqualStrings("acc1", entry.name);
        try std.testing.expectEqual(@as(u64, 100), entry.offset);
        try std.testing.expectEqual(@as(u64, 110), entry.data_offset);
        try std.testing.expectEqual(@as(u64, 50), entry.seq_len);
    }

    // Missing secondary key returns null.
    {
        const result = try idx.lookup(allocator, "acc3");
        try std.testing.expect(result == null);
    }
}

test "EaselIndex.findNumber: ordinal access" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{ .nprimary = 3 });

    try writePrimaryKey(&buf, allocator, .{ .key = "alpha", .r_off = 0, .d_off = 10, .len = 100 });
    try writePrimaryKey(&buf, allocator, .{ .key = "beta", .r_off = 200, .d_off = 210, .len = 50 });
    try writePrimaryKey(&buf, allocator, .{ .key = "gamma", .r_off = 400, .d_off = 410, .len = 75 });

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(.{ .sub_path = "test.ssi", .data = buf.items });

    const file = try tmp_dir.dir.openFile("test.ssi", .{});
    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    const e0 = (try idx.findNumber(allocator, 0)) orelse return error.ExpectedEntry;
    defer allocator.free(e0.name);
    try std.testing.expectEqualStrings("alpha", e0.name);
    try std.testing.expectEqual(@as(u64, 0), e0.offset);
    try std.testing.expectEqual(@as(u64, 10), e0.data_offset);
    try std.testing.expectEqual(@as(u64, 100), e0.seq_len);

    const e2 = (try idx.findNumber(allocator, 2)) orelse return error.ExpectedEntry;
    defer allocator.free(e2.name);
    try std.testing.expectEqualStrings("gamma", e2.name);
    try std.testing.expectEqual(@as(u64, 400), e2.offset);

    // Out of range
    const missing = try idx.findNumber(allocator, 3);
    try std.testing.expect(missing == null);
}

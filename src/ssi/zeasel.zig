// ZeaselIndex — zeasel-native SSI (Simple Sequence Index) format.
//
// Binary format:
//   Header:
//     magic:       [4]u8  = "ZSSI"
//     version:     u32    = 1   (little-endian)
//     num_entries: u64          (little-endian)
//   Per entry:
//     name_len:    u32          (little-endian)
//     name:        [name_len]u8
//     offset:      u64          (little-endian) — byte offset of '>' header
//     data_offset: u64          (little-endian) — byte offset of first residue
//     seq_len:     u64          (little-endian) — number of residues (excl. gaps/whitespace)

const std = @import("std");
const Allocator = std.mem.Allocator;

const ssi = @import("../ssi.zig");
const SsiEntry = ssi.SsiEntry;

const MAGIC = "ZSSI";
const VERSION: u32 = 2; // v2 adds file_id, aliases, file_names

pub const ZeaselIndex = struct {
    entries: []SsiEntry,
    name_map: std.StringHashMap(usize),
    alias_map: std.StringHashMap(usize), // Secondary key -> primary entry index
    file_names: [][]const u8, // Multi-file support: file_names[file_id]
    allocator: Allocator,

    /// Build an index from a FASTA byte buffer.
    /// Scans the data recording byte offsets of each '>' header.
    pub fn buildFromFasta(allocator: Allocator, data: []const u8) !ZeaselIndex {
        var entry_list = std.ArrayList(SsiEntry){};
        errdefer {
            for (entry_list.items) |e| allocator.free(e.name);
            entry_list.deinit(allocator);
        }

        var pos: usize = 0;
        while (pos < data.len) {
            // Skip whitespace / blank lines between records.
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

            const record_offset = pos;

            // Advance past '>'.
            pos += 1;

            // Parse name (up to first whitespace or newline).
            const name_start = pos;
            while (pos < data.len and data[pos] != ' ' and data[pos] != '\t' and
                data[pos] != '\n' and data[pos] != '\r')
            {
                pos += 1;
            }
            const name = try allocator.dupe(u8, data[name_start..pos]);
            errdefer allocator.free(name);

            // Skip rest of header line.
            while (pos < data.len and data[pos] != '\n' and data[pos] != '\r') {
                pos += 1;
            }
            if (pos < data.len and data[pos] == '\r') pos += 1;
            if (pos < data.len and data[pos] == '\n') pos += 1;

            const data_start = pos;

            // Count residues and advance to next '>'.
            var seq_len: u64 = 0;
            while (pos < data.len and data[pos] != '>') {
                const c = data[pos];
                if (c != '\n' and c != '\r' and c != ' ' and c != '\t') {
                    seq_len += 1;
                }
                pos += 1;
            }

            try entry_list.append(allocator, SsiEntry{
                .name = name,
                .offset = record_offset,
                .data_offset = data_start,
                .seq_len = seq_len,
            });
        }

        const entries = try entry_list.toOwnedSlice(allocator);
        errdefer {
            for (entries) |e| allocator.free(e.name);
            allocator.free(entries);
        }

        var name_map = std.StringHashMap(usize).init(allocator);
        errdefer name_map.deinit();

        for (entries, 0..) |e, i| {
            try name_map.put(e.name, i);
        }

        return ZeaselIndex{
            .entries = entries,
            .name_map = name_map,
            .alias_map = std.StringHashMap(usize).init(allocator),
            .file_names = &.{},
            .allocator = allocator,
        };
    }

    /// Look up a sequence by primary name or alias. Returns the entry or null.
    pub fn lookup(self: ZeaselIndex, name: []const u8) ?SsiEntry {
        const idx = self.name_map.get(name) orelse
            self.alias_map.get(name) orelse return null;
        return self.entries[idx];
    }

    /// Look up a sequence by ordinal number (0-indexed).
    pub fn findNumber(self: ZeaselIndex, number: usize) ?SsiEntry {
        if (number >= self.entries.len) return null;
        return self.entries[number];
    }

    /// Add an alias (secondary key) mapping to an existing primary entry.
    pub fn addAlias(self: *ZeaselIndex, alias: []const u8, primary_name: []const u8) !void {
        const idx = self.name_map.get(primary_name) orelse return error.KeyNotFound;
        const alias_copy = try self.allocator.dupe(u8, alias);
        errdefer self.allocator.free(alias_copy);
        try self.alias_map.put(alias_copy, idx);
    }

    /// Add a file to the multi-file index. Returns the file_id.
    pub fn addFile(self: *ZeaselIndex, file_name: []const u8) !u16 {
        const new_id: u16 = @intCast(self.file_names.len);
        const name_copy = try self.allocator.dupe(u8, file_name);
        errdefer self.allocator.free(name_copy);

        const new_list = try self.allocator.alloc([]const u8, self.file_names.len + 1);
        @memcpy(new_list[0..self.file_names.len], self.file_names);
        new_list[self.file_names.len] = name_copy;
        if (self.file_names.len > 0) self.allocator.free(self.file_names);
        self.file_names = new_list;

        return new_id;
    }

    /// Get the file name for a given file_id.
    pub fn fileInfo(self: ZeaselIndex, file_id: u16) ?[]const u8 {
        if (file_id >= self.file_names.len) return null;
        return self.file_names[file_id];
    }

    /// Compute subsequence position: given a sequence name, start, and end
    /// residue coordinates (1-indexed), returns the entry with positioning info.
    /// The caller can use data_offset + computed byte offset for seeking.
    pub fn findSubseq(self: ZeaselIndex, name: []const u8, start: u64, end: u64) ?struct { entry: SsiEntry, start: u64, end: u64 } {
        const entry = self.lookup(name) orelse return null;
        if (start == 0 or end == 0 or start > entry.seq_len or end > entry.seq_len) return null;
        return .{ .entry = entry, .start = start, .end = end };
    }

    /// Write the index to a binary stream (little-endian).
    /// v2 format: entries with file_id, file_names, alias_map.
    pub fn write(self: ZeaselIndex, dest: std.io.AnyWriter) !void {
        // Header: magic + version + num_entries + num_files + num_aliases
        try dest.writeAll(MAGIC);
        try dest.writeInt(u32, VERSION, .little);
        try dest.writeInt(u64, @intCast(self.entries.len), .little);
        try dest.writeInt(u32, @intCast(self.file_names.len), .little);
        try dest.writeInt(u32, @intCast(self.alias_map.count()), .little);

        // Entries (with file_id)
        for (self.entries) |e| {
            try dest.writeInt(u32, @intCast(e.name.len), .little);
            try dest.writeAll(e.name);
            try dest.writeInt(u64, e.offset, .little);
            try dest.writeInt(u64, e.data_offset, .little);
            try dest.writeInt(u64, e.seq_len, .little);
            try dest.writeInt(u16, e.file_id, .little);
        }

        // File names
        for (self.file_names) |fname| {
            try dest.writeInt(u32, @intCast(fname.len), .little);
            try dest.writeAll(fname);
        }

        // Aliases: each is (alias_name_len, alias_name, entry_index)
        var alias_it = self.alias_map.iterator();
        while (alias_it.next()) |kv| {
            try dest.writeInt(u32, @intCast(kv.key_ptr.len), .little);
            try dest.writeAll(kv.key_ptr.*);
            try dest.writeInt(u64, @intCast(kv.value_ptr.*), .little);
        }
    }

    /// Read an index from a binary stream.
    /// Supports both v1 (no file_id/aliases/file_names) and v2 formats.
    /// Caller owns the returned ZeaselIndex and must call deinit().
    pub fn read(allocator: Allocator, src: std.io.AnyReader) !ZeaselIndex {
        var magic: [4]u8 = undefined;
        _ = try src.readAll(&magic);
        if (!std.mem.eql(u8, &magic, MAGIC)) return error.InvalidFormat;

        const version = try src.readInt(u32, .little);
        if (version != 1 and version != VERSION) return error.UnsupportedVersion;

        const num_entries = try src.readInt(u64, .little);

        // v2 has additional header fields
        const num_files: u32 = if (version >= 2) try src.readInt(u32, .little) else 0;
        const num_aliases: u32 = if (version >= 2) try src.readInt(u32, .little) else 0;

        const entries = try allocator.alloc(SsiEntry, num_entries);
        var entries_done: usize = 0;
        errdefer {
            for (0..entries_done) |i| allocator.free(entries[i].name);
            allocator.free(entries);
        }

        for (0..num_entries) |i| {
            const name_len = try src.readInt(u32, .little);
            const name = try allocator.alloc(u8, name_len);
            errdefer allocator.free(name);
            _ = try src.readAll(name);

            const offset = try src.readInt(u64, .little);
            const data_offset = try src.readInt(u64, .little);
            const seq_len = try src.readInt(u64, .little);
            const file_id: u16 = if (version >= 2) try src.readInt(u16, .little) else 0;

            entries[i] = SsiEntry{
                .name = name,
                .offset = offset,
                .data_offset = data_offset,
                .seq_len = seq_len,
                .file_id = file_id,
            };
            entries_done += 1;
        }

        // Read file names
        var file_names: [][]const u8 = &.{};
        if (num_files > 0) {
            file_names = try allocator.alloc([]const u8, num_files);
            var files_done: usize = 0;
            errdefer {
                for (0..files_done) |fi| allocator.free(file_names[fi]);
                allocator.free(file_names);
            }
            for (0..num_files) |fi| {
                const fname_len = try src.readInt(u32, .little);
                const fname = try allocator.alloc(u8, fname_len);
                errdefer allocator.free(fname);
                _ = try src.readAll(fname);
                file_names[fi] = fname;
                files_done += 1;
            }
        }

        var name_map = std.StringHashMap(usize).init(allocator);
        errdefer name_map.deinit();
        for (entries, 0..) |e, i| {
            try name_map.put(e.name, i);
        }

        // Read aliases
        var alias_map = std.StringHashMap(usize).init(allocator);
        errdefer {
            var ait = alias_map.keyIterator();
            while (ait.next()) |k| allocator.free(k.*);
            alias_map.deinit();
        }
        for (0..num_aliases) |_| {
            const alias_len = try src.readInt(u32, .little);
            const alias_name = try allocator.alloc(u8, alias_len);
            errdefer allocator.free(alias_name);
            _ = try src.readAll(alias_name);
            const entry_idx = try src.readInt(u64, .little);
            try alias_map.put(alias_name, entry_idx);
        }

        return ZeaselIndex{
            .entries = entries,
            .name_map = name_map,
            .alias_map = alias_map,
            .file_names = file_names,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *ZeaselIndex) void {
        for (self.entries) |e| self.allocator.free(e.name);
        self.allocator.free(self.entries);
        self.name_map.deinit();
        // Free alias keys (they were duped)
        var alias_it = self.alias_map.keyIterator();
        while (alias_it.next()) |key| self.allocator.free(key.*);
        self.alias_map.deinit();
        for (self.file_names) |f| self.allocator.free(f);
        if (self.file_names.len > 0) self.allocator.free(self.file_names);
    }
};

/// Maximum number of entries before returning error.TooManyKeys.
/// External sort for larger indexes can be implemented later.
pub const MAX_BUILDER_ENTRIES: usize = 10_000_000;

/// Incremental SSI index builder.
///
/// Allows adding entries from arbitrary formats one at a time,
/// registering multiple source files, and adding aliases.
/// Call `build()` to produce a `ZeaselIndex` or `write()` to
/// serialise directly without an intermediate in-memory index.
pub const IndexBuilder = struct {
    allocator: Allocator,
    entries: std.ArrayList(Entry),
    aliases: std.ArrayList(Alias),
    file_names: std.ArrayList([]const u8),

    pub const Entry = struct {
        name: []const u8,
        offset: u64,
        data_offset: u64,
        seq_len: u64,
        file_id: u16,
    };

    pub const Alias = struct {
        alias: []const u8,
        target: []const u8,
    };

    pub fn init(allocator: Allocator) IndexBuilder {
        return .{
            .allocator = allocator,
            .entries = std.ArrayList(Entry){},
            .aliases = std.ArrayList(Alias){},
            .file_names = std.ArrayList([]const u8){},
        };
    }

    /// Register a source file. Returns the file_id (0-indexed).
    pub fn addFile(self: *IndexBuilder, filename: []const u8) !u16 {
        const name_copy = try self.allocator.dupe(u8, filename);
        errdefer self.allocator.free(name_copy);
        try self.file_names.append(self.allocator, name_copy);
        return @intCast(self.file_names.items.len - 1);
    }

    /// Add a primary key entry.
    pub fn addKey(self: *IndexBuilder, name: []const u8, file_id: u16, offset: u64, data_offset: u64, seq_len: u64) !void {
        if (self.entries.items.len >= MAX_BUILDER_ENTRIES) return error.TooManyKeys;
        const name_copy = try self.allocator.dupe(u8, name);
        errdefer self.allocator.free(name_copy);
        try self.entries.append(self.allocator, .{
            .name = name_copy,
            .offset = offset,
            .data_offset = data_offset,
            .seq_len = seq_len,
            .file_id = file_id,
        });
    }

    /// Add an alias (secondary key pointing to a primary key name).
    pub fn addAlias(self: *IndexBuilder, alias: []const u8, primary: []const u8) !void {
        const alias_copy = try self.allocator.dupe(u8, alias);
        errdefer self.allocator.free(alias_copy);
        const target_copy = try self.allocator.dupe(u8, primary);
        errdefer self.allocator.free(target_copy);
        try self.aliases.append(self.allocator, .{
            .alias = alias_copy,
            .target = target_copy,
        });
    }

    /// Sort entries by name and check for duplicate primary keys.
    /// Returns error.DuplicateKey if any two entries share a name.
    fn validateAndSort(self: *IndexBuilder) !void {
        const items = self.entries.items;
        std.mem.sort(Entry, items, {}, struct {
            fn lessThan(_: void, a: Entry, b: Entry) bool {
                return std.mem.order(u8, a.name, b.name) == .lt;
            }
        }.lessThan);
        // Check adjacent pairs for duplicates.
        if (items.len > 1) {
            for (0..items.len - 1) |i| {
                if (std.mem.eql(u8, items[i].name, items[i + 1].name)) {
                    return error.DuplicateKey;
                }
            }
        }
    }

    /// Build a `ZeaselIndex` from the accumulated entries.
    /// Validates uniqueness of primary keys. Caller owns the returned index.
    pub fn build(self: *IndexBuilder) !ZeaselIndex {
        try self.validateAndSort();

        // Transfer entries into SsiEntry slice.
        const entries = try self.allocator.alloc(SsiEntry, self.entries.items.len);
        errdefer {
            for (entries) |e| self.allocator.free(e.name);
            self.allocator.free(entries);
        }

        for (self.entries.items, 0..) |be, i| {
            entries[i] = SsiEntry{
                .name = try self.allocator.dupe(u8, be.name),
                .offset = be.offset,
                .data_offset = be.data_offset,
                .seq_len = be.seq_len,
                .file_id = be.file_id,
            };
        }

        // Build name_map.
        var name_map = std.StringHashMap(usize).init(self.allocator);
        errdefer name_map.deinit();
        for (entries, 0..) |e, i| {
            try name_map.put(e.name, i);
        }

        // Build alias_map: resolve alias -> primary entry index.
        var alias_map = std.StringHashMap(usize).init(self.allocator);
        errdefer {
            var ait = alias_map.keyIterator();
            while (ait.next()) |k| self.allocator.free(k.*);
            alias_map.deinit();
        }
        for (self.aliases.items) |a| {
            const idx = name_map.get(a.target) orelse return error.KeyNotFound;
            const alias_copy = try self.allocator.dupe(u8, a.alias);
            errdefer self.allocator.free(alias_copy);
            try alias_map.put(alias_copy, idx);
        }

        // Copy file names.
        var file_names: [][]const u8 = &.{};
        if (self.file_names.items.len > 0) {
            file_names = try self.allocator.alloc([]const u8, self.file_names.items.len);
            var files_done: usize = 0;
            errdefer {
                for (0..files_done) |fi| self.allocator.free(file_names[fi]);
                self.allocator.free(file_names);
            }
            for (self.file_names.items, 0..) |fname, fi| {
                file_names[fi] = try self.allocator.dupe(u8, fname);
                files_done += 1;
            }
        }

        return ZeaselIndex{
            .entries = entries,
            .name_map = name_map,
            .alias_map = alias_map,
            .file_names = file_names,
            .allocator = self.allocator,
        };
    }

    /// Write the index directly to a binary stream in zeasel SSI v2 format,
    /// without building an intermediate ZeaselIndex.
    /// Validates uniqueness of primary keys.
    pub fn write(self: *IndexBuilder, dest: std.io.AnyWriter) !void {
        try self.validateAndSort();

        // To write aliases we need to resolve target names to sorted entry indices.
        // Build a temporary name->index map.
        var name_map = std.StringHashMap(usize).init(self.allocator);
        defer name_map.deinit();
        for (self.entries.items, 0..) |e, i| {
            try name_map.put(e.name, i);
        }

        // Header
        try dest.writeAll(MAGIC);
        try dest.writeInt(u32, VERSION, .little);
        try dest.writeInt(u64, @intCast(self.entries.items.len), .little);
        try dest.writeInt(u32, @intCast(self.file_names.items.len), .little);
        try dest.writeInt(u32, @intCast(self.aliases.items.len), .little);

        // Entries (sorted)
        for (self.entries.items) |e| {
            try dest.writeInt(u32, @intCast(e.name.len), .little);
            try dest.writeAll(e.name);
            try dest.writeInt(u64, e.offset, .little);
            try dest.writeInt(u64, e.data_offset, .little);
            try dest.writeInt(u64, e.seq_len, .little);
            try dest.writeInt(u16, e.file_id, .little);
        }

        // File names
        for (self.file_names.items) |fname| {
            try dest.writeInt(u32, @intCast(fname.len), .little);
            try dest.writeAll(fname);
        }

        // Aliases
        for (self.aliases.items) |a| {
            const idx = name_map.get(a.target) orelse return error.KeyNotFound;
            try dest.writeInt(u32, @intCast(a.alias.len), .little);
            try dest.writeAll(a.alias);
            try dest.writeInt(u64, @intCast(idx), .little);
        }
    }

    pub fn deinit(self: *IndexBuilder) void {
        for (self.entries.items) |e| self.allocator.free(e.name);
        self.entries.deinit(self.allocator);
        for (self.aliases.items) |a| {
            self.allocator.free(a.alias);
            self.allocator.free(a.target);
        }
        self.aliases.deinit(self.allocator);
        for (self.file_names.items) |f| self.allocator.free(f);
        self.file_names.deinit(self.allocator);
    }
};

// --- Tests ---

// Use ZeaselIndex directly for tests (SsiIndex is now a tagged union).
const SsiIndex = ZeaselIndex;

test "buildFromFasta: correct offsets" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectEqual(@as(usize, 2), idx.entries.len);
    try std.testing.expectEqualStrings("seq1", idx.entries[0].name);
    try std.testing.expectEqualStrings("seq2", idx.entries[1].name);

    // ">seq1" starts at byte 0.
    try std.testing.expectEqual(@as(u64, 0), idx.entries[0].offset);
    // Data starts after ">seq1\n" = 6 bytes.
    try std.testing.expectEqual(@as(u64, 6), idx.entries[0].data_offset);
    // seq_len is 4 (ACGT).
    try std.testing.expectEqual(@as(u64, 4), idx.entries[0].seq_len);

    // ">seq2" starts at byte 11 (">seq1\nACGT\n" = 11 bytes).
    try std.testing.expectEqual(@as(u64, 11), idx.entries[1].offset);
    try std.testing.expectEqual(@as(u64, 4), idx.entries[1].seq_len);
}

test "lookup: finds entry by name" {
    const allocator = std.testing.allocator;
    const data = ">alpha\nACGT\n>beta\nGGGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    const e = idx.lookup("alpha") orelse return error.TestEntryNotFound;
    try std.testing.expectEqualStrings("alpha", e.name);
    try std.testing.expectEqual(@as(u64, 4), e.seq_len);
}

test "lookup: missing name returns null" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectEqual(@as(?SsiEntry, null), idx.lookup("missing"));
}

test "write and read: round trip" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try idx.write(buf.writer(allocator).any());

    var stream = std.io.fixedBufferStream(buf.items);
    var idx2 = try SsiIndex.read(allocator, stream.reader().any());
    defer idx2.deinit();

    try std.testing.expectEqual(idx.entries.len, idx2.entries.len);
    for (0..idx.entries.len) |i| {
        try std.testing.expectEqualStrings(idx.entries[i].name, idx2.entries[i].name);
        try std.testing.expectEqual(idx.entries[i].offset, idx2.entries[i].offset);
        try std.testing.expectEqual(idx.entries[i].data_offset, idx2.entries[i].data_offset);
        try std.testing.expectEqual(idx.entries[i].seq_len, idx2.entries[i].seq_len);
    }

    // lookup works after round trip.
    const e = idx2.lookup("seq1") orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u64, 0), e.offset);
    try std.testing.expectEqual(@as(u64, 4), e.seq_len);
}

test "buildFromFasta: multi-line sequences" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\nACGT\n>seq2\nGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectEqual(@as(u64, 8), idx.entries[0].seq_len);
    try std.testing.expectEqual(@as(u64, 2), idx.entries[1].seq_len);
}

test "addAlias: lookup by alias" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try idx.addAlias("alias1", "seq1");
    const e = idx.lookup("alias1") orelse return error.TestEntryNotFound;
    try std.testing.expectEqualStrings("seq1", e.name);
    try std.testing.expectEqual(@as(u64, 4), e.seq_len);
}

test "addAlias: unknown primary returns error" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectError(error.KeyNotFound, idx.addAlias("alias1", "nonexistent"));
}

test "findNumber: ordinal lookup" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    const e0 = idx.findNumber(0) orelse return error.TestEntryNotFound;
    try std.testing.expectEqualStrings("seq1", e0.name);

    const e1 = idx.findNumber(1) orelse return error.TestEntryNotFound;
    try std.testing.expectEqualStrings("seq2", e1.name);

    try std.testing.expect(idx.findNumber(2) == null);
}

test "addFile and fileInfo" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    const id0 = try idx.addFile("data/file1.fasta");
    const id1 = try idx.addFile("data/file2.fasta");

    try std.testing.expectEqual(@as(u16, 0), id0);
    try std.testing.expectEqual(@as(u16, 1), id1);
    try std.testing.expectEqualStrings("data/file1.fasta", idx.fileInfo(0).?);
    try std.testing.expectEqualStrings("data/file2.fasta", idx.fileInfo(1).?);
    try std.testing.expect(idx.fileInfo(2) == null);
}

test "findSubseq: valid range" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGTACGT\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    const result = idx.findSubseq("seq1", 3, 6) orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u64, 3), result.start);
    try std.testing.expectEqual(@as(u64, 6), result.end);
}

test "findSubseq: out of range returns null" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expect(idx.findSubseq("seq1", 1, 10) == null);
    try std.testing.expect(idx.findSubseq("seq1", 0, 2) == null); // 0 is invalid (1-indexed)
}

test "buildFromFasta: sequence with description" {
    const allocator = std.testing.allocator;
    const data = ">seq1 my description\nACGT\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectEqualStrings("seq1", idx.entries[0].name);
    try std.testing.expectEqual(@as(u64, 4), idx.entries[0].seq_len);
}

// --- IndexBuilder tests ---

test "IndexBuilder: build from incremental additions" {
    const allocator = std.testing.allocator;

    var builder = IndexBuilder.init(allocator);
    defer builder.deinit();

    _ = try builder.addFile("test.fasta");
    try builder.addKey("seq1", 0, 0, 5, 100);
    try builder.addKey("seq2", 0, 200, 205, 150);

    var idx = try builder.build();
    defer idx.deinit();

    try std.testing.expectEqual(@as(usize, 2), idx.entries.len);

    const e1 = idx.lookup("seq1") orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u64, 100), e1.seq_len);
    try std.testing.expectEqual(@as(u64, 5), e1.data_offset);

    const e2 = idx.lookup("seq2") orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u64, 150), e2.seq_len);
}

test "IndexBuilder: duplicate key detection" {
    const allocator = std.testing.allocator;

    var builder = IndexBuilder.init(allocator);
    defer builder.deinit();

    _ = try builder.addFile("test.fasta");
    try builder.addKey("seq1", 0, 0, 5, 100);
    try builder.addKey("seq1", 0, 200, 205, 150);

    try std.testing.expectError(error.DuplicateKey, builder.build());
}

test "IndexBuilder: write and read round-trip" {
    const allocator = std.testing.allocator;

    var builder = IndexBuilder.init(allocator);
    defer builder.deinit();

    _ = try builder.addFile("data/seqs.fasta");
    try builder.addKey("alpha", 0, 0, 10, 50);
    try builder.addKey("beta", 0, 100, 110, 75);

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);
    try builder.write(buf.writer(allocator).any());

    // Read back and verify.
    var stream = std.io.fixedBufferStream(buf.items);
    var idx = try ZeaselIndex.read(allocator, stream.reader().any());
    defer idx.deinit();

    try std.testing.expectEqual(@as(usize, 2), idx.entries.len);
    try std.testing.expectEqual(@as(usize, 1), idx.file_names.len);
    try std.testing.expectEqualStrings("data/seqs.fasta", idx.file_names[0]);

    const e = idx.lookup("alpha") orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u64, 50), e.seq_len);

    const e2 = idx.lookup("beta") orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u64, 75), e2.seq_len);
}

test "IndexBuilder: aliases" {
    const allocator = std.testing.allocator;

    var builder = IndexBuilder.init(allocator);
    defer builder.deinit();

    _ = try builder.addFile("test.fasta");
    try builder.addKey("seq1", 0, 0, 5, 100);
    try builder.addAlias("accession1", "seq1");

    var idx = try builder.build();
    defer idx.deinit();

    // Lookup by alias returns the primary entry.
    const e = idx.lookup("accession1") orelse return error.TestEntryNotFound;
    try std.testing.expectEqualStrings("seq1", e.name);
    try std.testing.expectEqual(@as(u64, 100), e.seq_len);
}

test "IndexBuilder: alias write round-trip" {
    const allocator = std.testing.allocator;

    var builder = IndexBuilder.init(allocator);
    defer builder.deinit();

    _ = try builder.addFile("test.fasta");
    try builder.addKey("seq1", 0, 0, 5, 100);
    try builder.addAlias("acc1", "seq1");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);
    try builder.write(buf.writer(allocator).any());

    var stream = std.io.fixedBufferStream(buf.items);
    var idx = try ZeaselIndex.read(allocator, stream.reader().any());
    defer idx.deinit();

    // Alias should survive serialization.
    const e = idx.lookup("acc1") orelse return error.TestEntryNotFound;
    try std.testing.expectEqualStrings("seq1", e.name);
}

test "IndexBuilder: multiple files" {
    const allocator = std.testing.allocator;

    var builder = IndexBuilder.init(allocator);
    defer builder.deinit();

    const id0 = try builder.addFile("file1.fasta");
    const id1 = try builder.addFile("file2.fasta");

    try std.testing.expectEqual(@as(u16, 0), id0);
    try std.testing.expectEqual(@as(u16, 1), id1);

    try builder.addKey("seq_a", 0, 0, 5, 100);
    try builder.addKey("seq_b", 1, 0, 10, 200);

    var idx = try builder.build();
    defer idx.deinit();

    const ea = idx.lookup("seq_a") orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u16, 0), ea.file_id);

    const eb = idx.lookup("seq_b") orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u16, 1), eb.file_id);

    try std.testing.expectEqualStrings("file1.fasta", idx.fileInfo(0).?);
    try std.testing.expectEqualStrings("file2.fasta", idx.fileInfo(1).?);
}

test "IndexBuilder: alias to nonexistent primary fails on build" {
    const allocator = std.testing.allocator;

    var builder = IndexBuilder.init(allocator);
    defer builder.deinit();

    try builder.addKey("seq1", 0, 0, 5, 100);
    try builder.addAlias("alias1", "nonexistent");

    try std.testing.expectError(error.KeyNotFound, builder.build());
}

test "IndexBuilder: entries sorted by name after build" {
    const allocator = std.testing.allocator;

    var builder = IndexBuilder.init(allocator);
    defer builder.deinit();

    // Add in reverse order.
    try builder.addKey("zebra", 0, 0, 5, 10);
    try builder.addKey("alpha", 0, 100, 105, 20);
    try builder.addKey("middle", 0, 200, 205, 30);

    var idx = try builder.build();
    defer idx.deinit();

    // Entries should be sorted by name.
    try std.testing.expectEqualStrings("alpha", idx.entries[0].name);
    try std.testing.expectEqualStrings("middle", idx.entries[1].name);
    try std.testing.expectEqualStrings("zebra", idx.entries[2].name);
}

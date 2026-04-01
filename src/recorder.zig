// Line-based I/O recorder with rewind/replay capability.
//
// Saves a history of lines read from an input stream, allowing
// position marking and rollback for look-ahead parsing.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const Recorder = struct {
    lines: std.ArrayList([]const u8),
    position: usize,
    allocator: Allocator,
    max_lines: ?usize,
    /// Total number of lines ever written (may exceed lines.items.len for bounded recorders).
    total_written: usize,

    /// Create a new unbounded recorder.
    pub fn init(allocator: Allocator) Recorder {
        return .{
            .lines = .empty,
            .position = 0,
            .allocator = allocator,
            .max_lines = null,
            .total_written = 0,
        };
    }

    /// Create a bounded recorder that keeps only the most recent `max` entries.
    /// Uses a circular buffer approach: once full, the oldest entry is freed
    /// and overwritten.
    pub fn initBounded(allocator: Allocator, max: usize) Recorder {
        return .{
            .lines = .empty,
            .position = 0,
            .allocator = allocator,
            .max_lines = max,
            .total_written = 0,
        };
    }

    /// Record a line (takes ownership of a copy).
    pub fn recordLine(self: *Recorder, line: []const u8) !void {
        const copy = try self.allocator.dupe(u8, line);
        try self.appendInternal(self.allocator, copy);
    }

    /// Read the next line from the recording, advancing position.
    /// Returns null at end of recording.
    pub fn readLine(self: *Recorder) ?[]const u8 {
        if (self.position >= self.lines.items.len) return null;
        const line = self.lines.items[self.position];
        self.position += 1;
        return line;
    }

    /// Mark the current position for later rollback.
    pub fn mark(self: Recorder) usize {
        return self.position;
    }

    /// Rollback to a previously marked position.
    pub fn rollback(self: *Recorder, saved_position: usize) void {
        self.position = saved_position;
    }

    /// Rewind to the beginning.
    pub fn rewind(self: *Recorder) void {
        self.position = 0;
    }

    /// Number of recorded lines currently held.
    pub fn lineCount(self: Recorder) usize {
        return self.lines.items.len;
    }

    /// Append a line to the recorder (stream integration for look-ahead buffering).
    /// The line is copied; the caller retains ownership of the input slice.
    /// Modelled after Easel's esl_recorder_Read + manual append pattern.
    pub fn addLine(self: *Recorder, allocator: Allocator, line: []const u8) !void {
        const copy = try allocator.dupe(u8, line);
        try self.appendInternal(allocator, copy);
    }

    /// Get a line by logical index. For bounded recorders the index is
    /// relative to the current buffer contents (0 = oldest retained line).
    pub fn getLine(self: Recorder, index: usize) ?[]const u8 {
        if (index >= self.lines.items.len) return null;
        return self.lines.items[index];
    }

    /// Internal helper: append an already-duped line, evicting the oldest
    /// entry if the bounded capacity is reached.
    fn appendInternal(self: *Recorder, allocator: Allocator, copy: []const u8) !void {
        if (self.max_lines) |max| {
            if (max == 0) {
                // Zero capacity: discard immediately.
                self.allocator.free(copy);
                self.total_written += 1;
                return;
            }
            if (self.lines.items.len >= max) {
                // Evict the oldest entry (index 0) and shift left.
                self.allocator.free(self.lines.items[0]);
                std.mem.copyForwards([]const u8, self.lines.items[0 .. self.lines.items.len - 1], self.lines.items[1..self.lines.items.len]);
                self.lines.items[self.lines.items.len - 1] = copy;
                // Adjust position to stay valid after the shift.
                if (self.position > 0) self.position -= 1;
                self.total_written += 1;
                return;
            }
        }
        try self.lines.append(allocator, copy);
        self.total_written += 1;
    }

    /// Mark the current end-of-recording position as the start of a block.
    /// Returns the position marker that can later be passed to getBlock().
    /// Modelled after Easel's esl_recorder_MarkBlock().
    pub fn markBlock(self: Recorder) usize {
        return self.lines.items.len;
    }

    /// Return all lines recorded from a previous block marker to the current
    /// end of the recording.  The returned slice is newly allocated and owned
    /// by the caller; the individual line slices point into recorder memory
    /// and must not be freed separately.
    /// Modelled after Easel's esl_recorder_GetBlock().
    pub fn getBlock(self: Recorder, allocator: Allocator, block_mark: usize) ![][]const u8 {
        if (block_mark > self.lines.items.len) return error.InvalidBlockMark;
        const block_lines = self.lines.items[block_mark..];
        const result = try allocator.alloc([]const u8, block_lines.len);
        @memcpy(result, block_lines);
        return result;
    }

    /// Free all recorded lines and the recorder.
    pub fn deinit(self: *Recorder) void {
        for (self.lines.items) |line| self.allocator.free(line);
        self.lines.deinit(self.allocator);
    }
};

// --- Tests ---

test "record and replay" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    try rec.recordLine("first line");
    try rec.recordLine("second line");
    try rec.recordLine("third line");

    try std.testing.expectEqual(@as(usize, 3), rec.lineCount());

    try std.testing.expectEqualStrings("first line", rec.readLine().?);
    try std.testing.expectEqualStrings("second line", rec.readLine().?);
    try std.testing.expectEqualStrings("third line", rec.readLine().?);
    try std.testing.expectEqual(@as(?[]const u8, null), rec.readLine());
}

test "rewind replays from beginning" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    try rec.recordLine("line1");
    try rec.recordLine("line2");

    _ = rec.readLine();
    _ = rec.readLine();
    rec.rewind();

    try std.testing.expectEqualStrings("line1", rec.readLine().?);
}

test "mark and rollback" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    try rec.recordLine("A");
    try rec.recordLine("B");
    try rec.recordLine("C");

    _ = rec.readLine(); // A
    const saved = rec.mark(); // position = 1
    _ = rec.readLine(); // B
    _ = rec.readLine(); // C
    rec.rollback(saved); // back to after A

    try std.testing.expectEqualStrings("B", rec.readLine().?);
}

test "addLine appends to recording" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    try rec.addLine(allocator, "hello");
    try rec.addLine(allocator, "world");

    try std.testing.expectEqual(@as(usize, 2), rec.lineCount());
    try std.testing.expectEqualStrings("hello", rec.readLine().?);
    try std.testing.expectEqualStrings("world", rec.readLine().?);
}

test "markBlock and getBlock: basic block capture" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    // Simulate reading a header then marking a block of data lines
    try rec.addLine(allocator, "# header");
    const block_start = rec.markBlock(); // mark after header
    try rec.addLine(allocator, "data1");
    try rec.addLine(allocator, "data2");
    try rec.addLine(allocator, "data3");

    const block = try rec.getBlock(allocator, block_start);
    defer allocator.free(block);

    try std.testing.expectEqual(@as(usize, 3), block.len);
    try std.testing.expectEqualStrings("data1", block[0]);
    try std.testing.expectEqualStrings("data2", block[1]);
    try std.testing.expectEqualStrings("data3", block[2]);
}

test "markBlock and getBlock: empty block" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    const block_start = rec.markBlock();
    const block = try rec.getBlock(allocator, block_start);
    defer allocator.free(block);

    try std.testing.expectEqual(@as(usize, 0), block.len);
}

test "markBlock and getBlock: multiple blocks" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    try rec.addLine(allocator, "A");
    const mark1 = rec.markBlock();
    try rec.addLine(allocator, "B");
    try rec.addLine(allocator, "C");
    const mark2 = rec.markBlock();
    try rec.addLine(allocator, "D");

    const block1 = try rec.getBlock(allocator, mark1);
    defer allocator.free(block1);
    const block2 = try rec.getBlock(allocator, mark2);
    defer allocator.free(block2);

    try std.testing.expectEqual(@as(usize, 3), block1.len);
    try std.testing.expectEqualStrings("B", block1[0]);
    try std.testing.expectEqualStrings("D", block1[2]);

    try std.testing.expectEqual(@as(usize, 1), block2.len);
    try std.testing.expectEqualStrings("D", block2[0]);
}

test "getBlock: invalid mark returns error" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    try rec.addLine(allocator, "only line");
    try std.testing.expectError(error.InvalidBlockMark, rec.getBlock(allocator, 999));
}

test "initBounded: keeps only most recent lines" {
    const allocator = std.testing.allocator;
    var rec = Recorder.initBounded(allocator, 3);
    defer rec.deinit();

    try rec.addLine(allocator, "A");
    try rec.addLine(allocator, "B");
    try rec.addLine(allocator, "C");
    try std.testing.expectEqual(@as(usize, 3), rec.lineCount());

    // Adding a 4th line should evict "A".
    try rec.addLine(allocator, "D");
    try std.testing.expectEqual(@as(usize, 3), rec.lineCount());
    try std.testing.expectEqualStrings("B", rec.getLine(0).?);
    try std.testing.expectEqualStrings("C", rec.getLine(1).?);
    try std.testing.expectEqualStrings("D", rec.getLine(2).?);

    // Adding another should evict "B".
    try rec.addLine(allocator, "E");
    try std.testing.expectEqual(@as(usize, 3), rec.lineCount());
    try std.testing.expectEqualStrings("C", rec.getLine(0).?);
    try std.testing.expectEqualStrings("D", rec.getLine(1).?);
    try std.testing.expectEqualStrings("E", rec.getLine(2).?);
}

test "initBounded: total_written tracks all lines" {
    const allocator = std.testing.allocator;
    var rec = Recorder.initBounded(allocator, 2);
    defer rec.deinit();

    try rec.addLine(allocator, "1");
    try rec.addLine(allocator, "2");
    try rec.addLine(allocator, "3");
    try rec.addLine(allocator, "4");

    try std.testing.expectEqual(@as(usize, 4), rec.total_written);
    try std.testing.expectEqual(@as(usize, 2), rec.lineCount());
}

test "initBounded: readLine works with eviction" {
    const allocator = std.testing.allocator;
    var rec = Recorder.initBounded(allocator, 2);
    defer rec.deinit();

    try rec.recordLine("first");
    try rec.recordLine("second");
    try rec.recordLine("third"); // evicts "first"

    try std.testing.expectEqualStrings("second", rec.readLine().?);
    try std.testing.expectEqualStrings("third", rec.readLine().?);
    try std.testing.expectEqual(@as(?[]const u8, null), rec.readLine());
}

test "initBounded: zero capacity discards everything" {
    const allocator = std.testing.allocator;
    var rec = Recorder.initBounded(allocator, 0);
    defer rec.deinit();

    try rec.addLine(allocator, "ignored");
    try rec.addLine(allocator, "also ignored");

    try std.testing.expectEqual(@as(usize, 0), rec.lineCount());
    try std.testing.expectEqual(@as(usize, 2), rec.total_written);
}

test "getLine: returns correct line by index" {
    const allocator = std.testing.allocator;
    var rec = Recorder.init(allocator);
    defer rec.deinit();

    try rec.addLine(allocator, "zero");
    try rec.addLine(allocator, "one");
    try rec.addLine(allocator, "two");

    try std.testing.expectEqualStrings("zero", rec.getLine(0).?);
    try std.testing.expectEqualStrings("one", rec.getLine(1).?);
    try std.testing.expectEqualStrings("two", rec.getLine(2).?);
    try std.testing.expect(rec.getLine(3) == null);
}

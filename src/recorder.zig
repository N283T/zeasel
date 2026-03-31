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

    /// Create a new recorder.
    pub fn init(allocator: Allocator) Recorder {
        return .{
            .lines = .empty,
            .position = 0,
            .allocator = allocator,
        };
    }

    /// Record a line (takes ownership of a copy).
    pub fn recordLine(self: *Recorder, line: []const u8) !void {
        const copy = try self.allocator.dupe(u8, line);
        try self.lines.append(self.allocator, copy);
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

    /// Number of recorded lines.
    pub fn lineCount(self: Recorder) usize {
        return self.lines.items.len;
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

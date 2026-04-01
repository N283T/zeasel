// zeasel-seqstat: report statistics about sequences in a file.
// Equivalent to Easel's esl-seqstat.
//
// Usage: zeasel-seqstat [-a] <seqfile>

const std = @import("std");
const zeasel = @import("zeasel");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 2) {
        std.debug.print("Usage: zeasel-seqstat [-a] <seqfile>\n", .{});
        std.process.exit(1);
    }

    // Parse optional flags.
    var show_per_seq = false;
    var positional_start: usize = 1;
    for (args[1..]) |arg| {
        if (std.mem.eql(u8, arg, "-a")) {
            show_per_seq = true;
            positional_start += 1;
        } else {
            break;
        }
    }

    if (positional_start >= args.len) {
        std.debug.print("Usage: zeasel-seqstat [-a] <seqfile>\n", .{});
        std.process.exit(1);
    }

    const path = args[positional_start];

    // Read entire file into memory.
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    const data = try file.readToEndAlloc(allocator, 512 * 1024 * 1024);
    defer allocator.free(data);

    // Detect format.
    const format = zeasel.io.Format.detect(data) orelse {
        std.debug.print("Error: cannot detect file format\n", .{});
        std.process.exit(1);
    };

    // Guess alphabet from file content.
    const abc_type = zeasel.alphabet.guessType(data) orelse blk: {
        std.debug.print("Warning: could not detect alphabet, defaulting to DNA\n", .{});
        break :blk .dna;
    };
    const abc: *const zeasel.alphabet.Alphabet = switch (abc_type) {
        .dna => &zeasel.alphabet.dna,
        .rna => &zeasel.alphabet.rna,
        .amino => &zeasel.alphabet.amino,
    };

    // Parse all sequences using the unified Reader.
    var reader = try zeasel.io.Reader.fromMemory(allocator, abc, data, format);
    defer reader.deinit();

    const sequences = try reader.readAll();
    defer {
        for (sequences) |*seq| @constCast(seq).deinit();
        allocator.free(sequences);
    }

    const stdout_file = std.fs.File.stdout();
    // TODO(M1): deprecatedWriter() is deprecated; migrate to stdout_file.writer(buffer)
    // once the codebase adopts the new std.Io.Writer API (requires an explicit buffer).
    const stdout = stdout_file.deprecatedWriter();

    if (sequences.len == 0) {
        try stdout.print("No sequences found.\n", .{});
        return;
    }

    // Per-sequence info (when -a is set).
    if (show_per_seq) {
        for (sequences) |seq| {
            const desc = seq.description orelse "";
            try stdout.print("= {s} {d} {s}\n", .{ seq.name, seq.len(), desc });
        }
        try stdout.print("\n", .{});
    }

    // Accumulate statistics.
    var total: u64 = 0;
    var shortest: u64 = std.math.maxInt(u64);
    var longest: u64 = 0;

    for (sequences) |seq| {
        const l: u64 = @intCast(seq.len());
        total += l;
        if (l < shortest) shortest = l;
        if (l > longest) longest = l;
    }

    const avg: f64 = @as(f64, @floatFromInt(total)) / @as(f64, @floatFromInt(sequences.len));

    try stdout.print("Format:          {s}\n", .{@tagName(format)});
    try stdout.print("Alphabet:        {s}\n", .{@tagName(abc_type)});
    try stdout.print("Sequences:       {d}\n", .{sequences.len});
    try stdout.print("Total residues:  {d}\n", .{total});
    try stdout.print("Shortest:        {d}\n", .{shortest});
    try stdout.print("Longest:         {d}\n", .{longest});
    try stdout.print("Average length:  {d:.1}\n", .{avg});
}

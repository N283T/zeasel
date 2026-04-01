// zeasel-seqfetch: fetch specific sequences by name from a file.
//
// Usage: zeasel-seqfetch <seqfile> <name1> [name2] ...
//        zeasel-seqfetch --index <seqfile>
// Outputs matching sequences in FASTA format to stdout.
// With --index, builds a .ssi index file for the given FASTA file.

const std = @import("std");
const zeasel = @import("zeasel");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 2) {
        std.debug.print("Usage: zeasel-seqfetch <seqfile> <name1> [name2] ...\n", .{});
        std.debug.print("       zeasel-seqfetch --index <seqfile>\n", .{});
        std.process.exit(1);
    }

    // Handle --index mode.
    if (std.mem.eql(u8, args[1], "--index")) {
        if (args.len < 3) {
            std.debug.print("Usage: zeasel-seqfetch --index <seqfile>\n", .{});
            std.process.exit(1);
        }
        try buildIndex(allocator, args[2]);
        return;
    }

    if (args.len < 3) {
        std.debug.print("Usage: zeasel-seqfetch <seqfile> <name1> [name2] ...\n", .{});
        std.debug.print("       zeasel-seqfetch --index <seqfile>\n", .{});
        std.process.exit(1);
    }

    const path = args[1];
    const names = args[2..];

    // Read entire file.
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

    // Parse all sequences.
    var reader = try zeasel.io.Reader.fromMemory(allocator, abc, data, format);
    defer reader.deinit();

    const sequences = try reader.readAll();
    defer {
        for (sequences) |*seq| @constCast(seq).deinit();
        allocator.free(sequences);
    }

    // Write matching sequences in FASTA format.
    const stdout_file = std.fs.File.stdout();
    // TODO(M1): deprecatedWriter() is deprecated; migrate to stdout_file.writer(buffer)
    // once the codebase adopts the new std.Io.Writer API (requires an explicit buffer).
    const stdout = stdout_file.deprecatedWriter();
    var writer = zeasel.io.Writer.init(stdout.any(), .fasta, 60);

    var found_count: usize = 0;
    for (names) |wanted| {
        var matched = false;
        for (sequences) |seq| {
            if (std.mem.eql(u8, seq.name, wanted)) {
                try writer.write(seq);
                matched = true;
                found_count += 1;
                break;
            }
        }
        if (!matched) {
            std.debug.print("Warning: sequence '{s}' not found\n", .{wanted});
        }
    }

    if (found_count == 0) {
        std.debug.print("No sequences matched.\n", .{});
        std.process.exit(1);
    }
}

/// Build an SSI index for a FASTA file and write it to <path>.ssi.
fn buildIndex(allocator: std.mem.Allocator, path: []const u8) !void {
    // Read the FASTA file.
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    const data = try file.readToEndAlloc(allocator, 512 * 1024 * 1024);
    defer allocator.free(data);

    // Build the index.
    var index = try zeasel.ssi.ZeaselIndex.buildFromFasta(allocator, data);
    defer index.deinit();

    // Write to <path>.ssi.
    const ssi_path = try std.fmt.allocPrint(allocator, "{s}.ssi", .{path});
    defer allocator.free(ssi_path);

    const out_file = try std.fs.cwd().createFile(ssi_path, .{});
    defer out_file.close();

    // TODO(M1): deprecatedWriter() is deprecated; migrate once the codebase
    // adopts the new std.Io.Writer API.
    const writer = out_file.deprecatedWriter();
    try index.write(writer.any());

    std.debug.print("SSI index written to {s} ({d} entries)\n", .{ ssi_path, index.entries.len });
}

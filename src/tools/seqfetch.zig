// zeacel-seqfetch: fetch specific sequences by name from a file.
//
// Usage: zeacel-seqfetch <seqfile> <name1> [name2] ...
// Outputs matching sequences in FASTA format to stdout.

const std = @import("std");
const zeacel = @import("zeacel");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 3) {
        std.debug.print("Usage: zeacel-seqfetch <seqfile> <name1> [name2] ...\n", .{});
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
    const format = zeacel.io.Format.detect(data) orelse {
        std.debug.print("Error: cannot detect file format\n", .{});
        std.process.exit(1);
    };

    // Guess alphabet from file content.
    const abc_type = zeacel.alphabet.guessType(data) orelse .dna;
    const abc: *const zeacel.alphabet.Alphabet = switch (abc_type) {
        .dna => &zeacel.alphabet.dna,
        .rna => &zeacel.alphabet.rna,
        .amino => &zeacel.alphabet.amino,
    };

    // Parse all sequences.
    var reader = try zeacel.io.Reader.fromMemory(allocator, abc, data, format);
    defer reader.deinit();

    const sequences = try reader.readAll();
    defer {
        for (sequences) |*seq| @constCast(seq).deinit();
        allocator.free(sequences);
    }

    // Write matching sequences in FASTA format.
    const stdout_file = std.fs.File.stdout();
    const stdout = stdout_file.deprecatedWriter();
    var writer = zeacel.io.Writer.init(stdout.any(), .fasta, 60);

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

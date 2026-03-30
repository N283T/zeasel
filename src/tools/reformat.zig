// zeacel-reformat: convert between sequence file formats.
// Equivalent to Easel's esl-reformat.
//
// Usage: zeacel-reformat <output-format> <input-file>
// Where output-format is one of: fasta, genbank, embl

const std = @import("std");
const zeacel = @import("zeacel");

fn parseOutputFormat(name: []const u8) ?zeacel.io.Format {
    if (std.mem.eql(u8, name, "fasta")) return .fasta;
    if (std.mem.eql(u8, name, "genbank")) return .genbank;
    if (std.mem.eql(u8, name, "embl")) return .embl;
    return null;
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 3) {
        std.debug.print("Usage: zeacel-reformat <output-format> <input-file>\n", .{});
        std.debug.print("Output formats: fasta, genbank, embl\n", .{});
        std.process.exit(1);
    }

    const out_fmt_name = args[1];
    const path = args[2];

    const out_format = parseOutputFormat(out_fmt_name) orelse {
        std.debug.print("Error: unknown output format '{s}'\n", .{out_fmt_name});
        std.debug.print("Output formats: fasta, genbank, embl\n", .{});
        std.process.exit(1);
    };

    // Read entire file.
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    const data = try file.readToEndAlloc(allocator, 512 * 1024 * 1024);
    defer allocator.free(data);

    // Detect input format.
    const in_format = zeacel.io.Format.detect(data) orelse {
        std.debug.print("Error: cannot detect input file format\n", .{});
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
    var reader = try zeacel.io.Reader.fromMemory(allocator, abc, data, in_format);
    defer reader.deinit();

    const sequences = try reader.readAll();
    defer {
        for (sequences) |*seq| @constCast(seq).deinit();
        allocator.free(sequences);
    }

    // Write sequences to stdout in the requested format.
    const stdout_file = std.fs.File.stdout();
    const stdout = stdout_file.deprecatedWriter();
    var writer = zeacel.io.Writer.init(stdout.any(), out_format, 60);
    try writer.writeAll(sequences);
}

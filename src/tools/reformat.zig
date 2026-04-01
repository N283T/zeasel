// zeasel-reformat: convert between sequence file formats.
// Equivalent to Easel's esl-reformat.
//
// Usage: zeasel-reformat <output-format> <input-file>
// Where output-format is one of: fasta, genbank, embl, stockholm, afa

const std = @import("std");
const zeasel = @import("zeasel");

const OutputKind = enum { sequence, msa };

fn parseOutputFormat(name: []const u8) ?struct { format: zeasel.io.Format, kind: OutputKind } {
    if (std.mem.eql(u8, name, "fasta")) return .{ .format = .fasta, .kind = .sequence };
    if (std.mem.eql(u8, name, "genbank")) return .{ .format = .genbank, .kind = .sequence };
    if (std.mem.eql(u8, name, "embl")) return .{ .format = .embl, .kind = .sequence };
    if (std.mem.eql(u8, name, "stockholm")) return .{ .format = .stockholm, .kind = .msa };
    if (std.mem.eql(u8, name, "afa")) return .{ .format = .afa, .kind = .msa };
    return null;
}

fn isMsaInputFormat(format: zeasel.io.Format) bool {
    return switch (format) {
        .stockholm, .clustal, .afa, .phylip, .a2m, .psiblast, .selex, .pfam => true,
        .fasta, .genbank, .embl, .ddbj => false,
    };
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 3) {
        std.debug.print("Usage: zeasel-reformat <output-format> <input-file>\n", .{});
        std.debug.print("Output formats: fasta, genbank, embl, stockholm, afa\n", .{});
        std.process.exit(1);
    }

    const out_fmt_name = args[1];
    const path = args[2];

    const out_info = parseOutputFormat(out_fmt_name) orelse {
        std.debug.print("Error: unknown output format '{s}'\n", .{out_fmt_name});
        std.debug.print("Output formats: fasta, genbank, embl, stockholm, afa\n", .{});
        std.process.exit(1);
    };

    // Read entire file.
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    const data = try file.readToEndAlloc(allocator, 512 * 1024 * 1024);
    defer allocator.free(data);

    // Detect input format.
    const in_format = zeasel.io.Format.detect(data) orelse {
        std.debug.print("Error: cannot detect input file format\n", .{});
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

    // Write output to stdout.
    const stdout_file = std.fs.File.stdout();
    // TODO(M1): deprecatedWriter() is deprecated; migrate to stdout_file.writer(buffer)
    // once the codebase adopts the new std.Io.Writer API (requires an explicit buffer).
    const stdout = stdout_file.deprecatedWriter();

    // MSA output formats require an Msa object. When both input and output are
    // MSA-capable, parse as MSA directly so alignment columns are preserved.
    if (out_info.kind == .msa and isMsaInputFormat(in_format)) {
        var msa = switch (in_format) {
            .stockholm, .pfam => try zeasel.io.stockholm.parse(allocator, abc, data),
            .clustal => try zeasel.io.clustal.parse(allocator, abc, data),
            .afa => try zeasel.io.afa.parse(allocator, abc, data),
            .phylip => try zeasel.io.phylip.parse(allocator, abc, data),
            .a2m => try zeasel.io.a2m.parse(allocator, abc, data),
            .psiblast => try zeasel.io.psiblast.parse(allocator, abc, data),
            .selex => try zeasel.io.selex.parse(allocator, abc, data),
            else => unreachable,
        };
        defer msa.deinit();

        switch (out_info.format) {
            .stockholm => try zeasel.io.stockholm.write(stdout.any(), msa),
            .afa => try zeasel.io.afa.write(stdout.any(), msa, 60),
            else => unreachable,
        }
    } else if (out_info.kind == .msa) {
        // MSA output requested but input is not an alignment format.
        std.debug.print("Error: output format '{s}' requires an alignment input (stockholm, clustal, afa, phylip, a2m, psiblast, selex)\n", .{out_fmt_name});
        std.process.exit(1);
    } else {
        // Sequence-based path: read individual sequences and write them out.
        var reader = try zeasel.io.Reader.fromMemory(allocator, abc, data, in_format);
        defer reader.deinit();

        const sequences = try reader.readAll();
        defer {
            for (sequences) |*seq| @constCast(seq).deinit();
            allocator.free(sequences);
        }

        var writer = zeasel.io.Writer.init(stdout.any(), out_info.format, 60);
        try writer.writeAll(sequences);
    }
}

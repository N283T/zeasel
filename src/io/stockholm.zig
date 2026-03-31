// Stockholm format parser and writer.
//
// Format specification:
//   - First line: "# STOCKHOLM 1.0"
//   - Last line: "//"
//   - #=GF <tag> <value>  : per-file annotation
//   - #=GS <name> <tag> <value> : per-sequence annotation
//   - #=GC <tag> <annotation>   : per-column annotation
//   - #=GR <name> <tag> <annotation> : per-residue annotation
//   - <name>  <sequence> : aligned sequence data
//   - Interleaved (multi-block) alignments are supported.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const msa_mod = @import("../msa.zig");
const Msa = msa_mod.Msa;
const GfEntry = msa_mod.GfEntry;
const GcEntry = msa_mod.GcEntry;
const GsEntry = msa_mod.GsEntry;
const GrEntry = msa_mod.GrEntry;

/// Parse a Stockholm format alignment from a byte buffer.
/// All markup lines (#=GF, #=GC, #=GS, #=GR) are stored in the Msa.
/// #=GC SS_cons is additionally stored in msa.consensus_ss.
/// #=GC RF is additionally stored in msa.reference.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    var lines = std.mem.splitScalar(u8, data, '\n');

    // Verify first non-blank line is the Stockholm header.
    var found_header = false;
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, "\r");
        if (trimmed.len == 0) continue;
        if (!std.mem.eql(u8, trimmed, "# STOCKHOLM 1.0")) return error.InvalidFormat;
        found_header = true;
        break;
    }
    if (!found_header) return error.InvalidFormat;

    // Metadata fields.
    var msa_name: ?[]const u8 = null;
    var msa_accession: ?[]const u8 = null;
    var msa_description: ?[]const u8 = null;
    errdefer {
        if (msa_name) |n| allocator.free(n);
        if (msa_accession) |a| allocator.free(a);
        if (msa_description) |d| allocator.free(d);
    }

    // Ordered list of names (first-seen order).
    var name_list = std.ArrayList([]const u8){};
    defer {
        for (name_list.items) |n| allocator.free(n);
        name_list.deinit(allocator);
    }

    // name_index maps sequence name slice (owned by name_list) to its index.
    var name_index = std.StringHashMap(usize).init(allocator);
    defer name_index.deinit();

    // seq_parts[i] collects all text fragments for sequence i across blocks.
    var seq_parts = std.ArrayList(std.ArrayList(u8)){};
    defer {
        for (seq_parts.items) |*parts| parts.deinit(allocator);
        seq_parts.deinit(allocator);
    }

    // Markup accumulators (unmanaged, allocator passed to each call).
    var gf_list = std.ArrayList(GfEntry){};
    errdefer {
        for (gf_list.items) |e| {
            allocator.free(e.tag);
            allocator.free(e.value);
        }
        gf_list.deinit(allocator);
    }

    var gc_list = std.ArrayList(GcEntry){};
    errdefer {
        for (gc_list.items) |e| {
            allocator.free(e.tag);
            allocator.free(e.annotation);
        }
        gc_list.deinit(allocator);
    }

    var gs_list = std.ArrayList(GsEntry){};
    errdefer {
        for (gs_list.items) |e| {
            allocator.free(e.seq_name);
            allocator.free(e.tag);
            allocator.free(e.value);
        }
        gs_list.deinit(allocator);
    }

    var gr_list = std.ArrayList(GrEntry){};
    errdefer {
        for (gr_list.items) |e| {
            allocator.free(e.seq_name);
            allocator.free(e.tag);
            allocator.free(e.annotation);
        }
        gr_list.deinit(allocator);
    }

    // For interleaved GC/GR, accumulate annotation fragments by key.
    // Key for GC: owned copy of tag.
    // Key for GR: owned string "<seqname>\x00<tag>".
    var gc_parts = std.StringHashMap(std.ArrayList(u8)).init(allocator);
    defer {
        var it = gc_parts.iterator();
        while (it.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            entry.value_ptr.deinit(allocator);
        }
        gc_parts.deinit();
    }

    var gr_parts = std.StringHashMap(std.ArrayList(u8)).init(allocator);
    defer {
        var it = gr_parts.iterator();
        while (it.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            entry.value_ptr.deinit(allocator);
        }
        gr_parts.deinit();
    }

    // Insertion-order tracking for GC/GR keys (slices into map keys — no extra alloc).
    var gc_order = std.ArrayList([]const u8){};
    defer gc_order.deinit(allocator);

    var gr_order = std.ArrayList([]const u8){};
    defer gr_order.deinit(allocator);

    var ended = false;

    while (lines.next()) |raw_line| {
        const line = std.mem.trimRight(u8, raw_line, "\r");

        // End of alignment record.
        if (std.mem.eql(u8, line, "//")) {
            ended = true;
            break;
        }

        // Blank line — block separator.
        if (line.len == 0) continue;

        // Markup lines.
        if (line[0] == '#') {
            if (std.mem.startsWith(u8, line, "#=GF ") or std.mem.startsWith(u8, line, "#=GF\t")) {
                const rest = trimLeft(line[5..]);
                const tag_end = indexOfWhitespace(rest) orelse rest.len;
                const tag = rest[0..tag_end];
                const value = if (tag_end < rest.len) trimLeft(rest[tag_end..]) else "";

                // Store in dedicated fields for ID/AC/DE.
                if (std.mem.eql(u8, tag, "ID")) {
                    if (msa_name) |old| allocator.free(old);
                    msa_name = try allocator.dupe(u8, value);
                } else if (std.mem.eql(u8, tag, "AC")) {
                    if (msa_accession) |old| allocator.free(old);
                    msa_accession = try allocator.dupe(u8, value);
                } else if (std.mem.eql(u8, tag, "DE")) {
                    if (msa_description) |old| allocator.free(old);
                    msa_description = try allocator.dupe(u8, value);
                }

                // Always store in gf_list for round-trip.
                const owned_tag = try allocator.dupe(u8, tag);
                errdefer allocator.free(owned_tag);
                const owned_value = try allocator.dupe(u8, value);
                errdefer allocator.free(owned_value);
                try gf_list.append(allocator, .{ .tag = owned_tag, .value = owned_value });
            } else if (std.mem.startsWith(u8, line, "#=GC ") or std.mem.startsWith(u8, line, "#=GC\t")) {
                const rest = trimLeft(line[5..]);
                const tag_end = indexOfWhitespace(rest) orelse rest.len;
                const tag = rest[0..tag_end];
                const annot = if (tag_end < rest.len) trimLeft(rest[tag_end..]) else "";

                // Accumulate annotation fragments for this tag.
                if (gc_parts.getPtr(tag)) |buf| {
                    try buf.appendSlice(allocator, annot);
                } else {
                    const owned_key = try allocator.dupe(u8, tag);
                    errdefer allocator.free(owned_key);
                    var buf = std.ArrayList(u8){};
                    errdefer buf.deinit(allocator);
                    try buf.appendSlice(allocator, annot);
                    try gc_parts.put(owned_key, buf);
                    try gc_order.append(allocator, owned_key);
                }
            } else if (std.mem.startsWith(u8, line, "#=GS ") or std.mem.startsWith(u8, line, "#=GS\t")) {
                const rest = trimLeft(line[5..]);
                // rest: "<seqname> <tag> <value>"
                const sn_end = indexOfWhitespace(rest) orelse rest.len;
                const seq_name = rest[0..sn_end];
                const after_sn = if (sn_end < rest.len) trimLeft(rest[sn_end..]) else "";
                const tag_end = indexOfWhitespace(after_sn) orelse after_sn.len;
                const tag = after_sn[0..tag_end];
                const value = if (tag_end < after_sn.len) trimLeft(after_sn[tag_end..]) else "";

                const owned_sn = try allocator.dupe(u8, seq_name);
                errdefer allocator.free(owned_sn);
                const owned_tag = try allocator.dupe(u8, tag);
                errdefer allocator.free(owned_tag);
                const owned_value = try allocator.dupe(u8, value);
                errdefer allocator.free(owned_value);
                try gs_list.append(allocator, .{ .seq_name = owned_sn, .tag = owned_tag, .value = owned_value });
            } else if (std.mem.startsWith(u8, line, "#=GR ") or std.mem.startsWith(u8, line, "#=GR\t")) {
                const rest = trimLeft(line[5..]);
                // rest: "<seqname> <tag> <annotation>"
                const sn_end = indexOfWhitespace(rest) orelse rest.len;
                const seq_name = rest[0..sn_end];
                const after_sn = if (sn_end < rest.len) trimLeft(rest[sn_end..]) else "";
                const tag_end = indexOfWhitespace(after_sn) orelse after_sn.len;
                const tag = after_sn[0..tag_end];
                const annot = if (tag_end < after_sn.len) trimLeft(after_sn[tag_end..]) else "";

                // Build composite key "seqname\x00tag".
                const composite_key = try std.fmt.allocPrint(allocator, "{s}\x00{s}", .{ seq_name, tag });
                errdefer allocator.free(composite_key);

                if (gr_parts.getPtr(composite_key)) |buf| {
                    allocator.free(composite_key); // key already stored in map
                    try buf.appendSlice(allocator, annot);
                } else {
                    var buf = std.ArrayList(u8){};
                    errdefer buf.deinit(allocator);
                    try buf.appendSlice(allocator, annot);
                    try gr_parts.put(composite_key, buf);
                    try gr_order.append(allocator, composite_key);
                }
            }
            continue;
        }

        // Sequence line: "<name> <sequence>"
        const name_end = indexOfWhitespace(line) orelse return error.InvalidFormat;
        const seq_name = line[0..name_end];
        const seq_text = trimLeft(line[name_end..]);

        if (seq_text.len == 0) return error.InvalidFormat;

        if (name_index.get(seq_name)) |idx| {
            // Append to existing sequence.
            try seq_parts.items[idx].appendSlice(allocator, seq_text);
        } else {
            // New sequence.
            const idx = name_list.items.len;
            const owned_name = try allocator.dupe(u8, seq_name);
            errdefer allocator.free(owned_name);

            try name_list.append(allocator, owned_name);
            errdefer _ = name_list.pop();

            try name_index.put(owned_name, idx);

            var parts = std.ArrayList(u8){};
            errdefer parts.deinit(allocator);

            try parts.appendSlice(allocator, seq_text);
            try seq_parts.append(allocator, parts);
        }
    }

    if (!ended) return error.InvalidFormat;
    if (name_list.items.len == 0) return error.InvalidInput;

    // Build slices for Msa.fromText.
    const n = name_list.items.len;

    const names_slice = try allocator.alloc([]const u8, n);
    defer allocator.free(names_slice);
    for (name_list.items, 0..) |name, i| {
        names_slice[i] = name;
    }

    const text_seqs = try allocator.alloc([]const u8, n);
    defer allocator.free(text_seqs);
    for (seq_parts.items, 0..) |*parts, i| {
        text_seqs[i] = parts.items;
    }

    var msa = try Msa.fromText(allocator, abc, names_slice, text_seqs);
    errdefer msa.deinit();

    // Transfer ownership of metadata strings.
    msa.name = msa_name;
    msa_name = null;
    msa.accession = msa_accession;
    msa_accession = null;
    msa.description = msa_description;
    msa_description = null;

    // Transfer gf_markup.
    if (gf_list.items.len > 0) {
        msa.gf_markup = try gf_list.toOwnedSlice(allocator);
    }
    gf_list.deinit(allocator);

    // Build gc_markup from gc_parts in insertion order.
    for (gc_order.items) |key| {
        const annot_buf = gc_parts.get(key).?;
        const owned_tag = try allocator.dupe(u8, key);
        errdefer allocator.free(owned_tag);
        const owned_annot = try allocator.dupe(u8, annot_buf.items);
        errdefer allocator.free(owned_annot);
        try gc_list.append(allocator, .{ .tag = owned_tag, .annotation = owned_annot });

        // Populate convenience fields.
        if (std.mem.eql(u8, key, "SS_cons")) {
            if (msa.consensus_ss) |old| allocator.free(old);
            msa.consensus_ss = try allocator.dupe(u8, annot_buf.items);
        } else if (std.mem.eql(u8, key, "RF")) {
            if (msa.reference) |old| allocator.free(old);
            msa.reference = try allocator.dupe(u8, annot_buf.items);
        }
    }
    if (gc_list.items.len > 0) {
        msa.gc_markup = try gc_list.toOwnedSlice(allocator);
    }
    gc_list.deinit(allocator);

    // Transfer gs_markup.
    if (gs_list.items.len > 0) {
        msa.gs_markup = try gs_list.toOwnedSlice(allocator);
    }
    gs_list.deinit(allocator);

    // Build gr_markup from gr_parts in insertion order.
    for (gr_order.items) |key| {
        const annot_buf = gr_parts.get(key).?;
        // Decode composite key: "seqname\x00tag"
        const sep = std.mem.indexOfScalar(u8, key, 0).?;
        const seq_nm = key[0..sep];
        const tag = key[sep + 1 ..];
        const owned_sn = try allocator.dupe(u8, seq_nm);
        errdefer allocator.free(owned_sn);
        const owned_tag = try allocator.dupe(u8, tag);
        errdefer allocator.free(owned_tag);
        const owned_annot = try allocator.dupe(u8, annot_buf.items);
        errdefer allocator.free(owned_annot);
        try gr_list.append(allocator, .{ .seq_name = owned_sn, .tag = owned_tag, .annotation = owned_annot });
    }
    if (gr_list.items.len > 0) {
        msa.gr_markup = try gr_list.toOwnedSlice(allocator);
    }
    gr_list.deinit(allocator);

    return msa;
}

/// Parse ALL MSA records from a multi-MSA Stockholm file.
/// Each record is delimited by "//". Returns a slice of Msa objects.
/// Caller owns the returned slice and each Msa within it.
pub fn parseAll(allocator: Allocator, abc: *const Alphabet, data: []const u8) ![]Msa {
    var result = std.ArrayList(Msa){};
    errdefer {
        for (result.items) |*m| m.deinit();
        result.deinit(allocator);
    }

    var pos: usize = 0;
    while (pos < data.len) {
        // Skip whitespace between records.
        while (pos < data.len and (data[pos] == '\n' or data[pos] == '\r' or data[pos] == ' ' or data[pos] == '\t')) {
            pos += 1;
        }
        if (pos >= data.len) break;

        // Find the end of this record ("//") at start of a line.
        const record_start = pos;
        var record_end: ?usize = null;
        var scan = pos;
        while (scan < data.len) {
            if (data[scan] == '/' and scan + 1 < data.len and data[scan + 1] == '/') {
                if (scan == 0 or data[scan - 1] == '\n') {
                    record_end = scan + 2;
                    if (record_end.? < data.len and data[record_end.?] == '\r') record_end = record_end.? + 1;
                    if (record_end.? < data.len and data[record_end.?] == '\n') record_end = record_end.? + 1;
                    break;
                }
            }
            scan += 1;
        }

        if (record_end == null) break;

        const record = data[record_start..record_end.?];
        var msa = try parse(allocator, abc, record);
        errdefer msa.deinit();
        try result.append(allocator, msa);
        pos = record_end.?;
    }

    return result.toOwnedSlice(allocator);
}

/// Write an Msa in Stockholm format.
/// Writes all markup stored in gf_markup/gc_markup/gs_markup/gr_markup.
/// If the dedicated name/accession/description fields are set but not covered by
/// gf_markup, they are written as #=GF ID/AC/DE lines.
pub fn write(dest: std.io.AnyWriter, m: Msa) !void {
    try dest.writeAll("# STOCKHOLM 1.0\n");

    const has_gf = m.gf_markup != null and m.gf_markup.?.len > 0;

    if (has_gf) {
        for (m.gf_markup.?) |e| {
            try dest.print("#=GF {s}   {s}\n", .{ e.tag, e.value });
        }
        // Emit dedicated fields only if absent from gf_markup.
        if (m.name) |nm| {
            if (!gfMarkupHasTag(m.gf_markup.?, "ID")) {
                try dest.print("#=GF ID   {s}\n", .{nm});
            }
        }
        if (m.accession) |acc| {
            if (!gfMarkupHasTag(m.gf_markup.?, "AC")) {
                try dest.print("#=GF AC   {s}\n", .{acc});
            }
        }
        if (m.description) |desc| {
            if (!gfMarkupHasTag(m.gf_markup.?, "DE")) {
                try dest.print("#=GF DE   {s}\n", .{desc});
            }
        }
    } else {
        // No gf_markup — fall back to dedicated fields.
        if (m.name) |nm| try dest.print("#=GF ID   {s}\n", .{nm});
        if (m.accession) |acc| try dest.print("#=GF AC   {s}\n", .{acc});
        if (m.description) |desc| try dest.print("#=GF DE   {s}\n", .{desc});
    }

    // Write GS lines before the blank separator and sequences.
    if (m.gs_markup) |entries| {
        for (entries) |e| {
            try dest.print("#=GS {s} {s} {s}\n", .{ e.seq_name, e.tag, e.value });
        }
    }

    try dest.writeByte('\n');

    // Compute column width for name padding.
    var max_name_len: usize = 0;
    for (m.names) |name| {
        if (name.len > max_name_len) max_name_len = name.len;
    }
    const col_width = max_name_len + 2;

    for (0..m.nseq()) |i| {
        const name = m.names[i];
        const text = try m.abc.textize(m.allocator, m.seqs[i]);
        defer m.allocator.free(text);

        try dest.writeAll(name);
        try writePadding(dest, col_width - name.len);
        try dest.writeAll(text);
        try dest.writeByte('\n');

        // Write #=GR lines for this sequence immediately after it.
        if (m.gr_markup) |entries| {
            for (entries) |e| {
                if (std.mem.eql(u8, e.seq_name, name)) {
                    try dest.print("#=GR {s} {s} {s}\n", .{ e.seq_name, e.tag, e.annotation });
                }
            }
        }
    }

    // Write #=GC lines after all sequences.
    if (m.gc_markup) |entries| {
        for (entries) |e| {
            try dest.print("#=GC {s} {s}\n", .{ e.tag, e.annotation });
        }
    }

    try dest.writeAll("//\n");
}

// --- Helpers ---

fn gfMarkupHasTag(entries: []const GfEntry, tag: []const u8) bool {
    for (entries) |e| {
        if (std.mem.eql(u8, e.tag, tag)) return true;
    }
    return false;
}

fn writePadding(dest: std.io.AnyWriter, count: usize) !void {
    var i: usize = 0;
    while (i < count) : (i += 1) {
        try dest.writeByte(' ');
    }
}

fn trimLeft(s: []const u8) []const u8 {
    var i: usize = 0;
    while (i < s.len and (s[i] == ' ' or s[i] == '\t')) i += 1;
    return s[i..];
}

fn indexOfWhitespace(s: []const u8) ?usize {
    for (s, 0..) |c, i| {
        if (c == ' ' or c == '\t') return i;
    }
    return null;
}

// --- Tests ---

const alphabet_mod = @import("../alphabet.zig");
const testing = std.testing;

test "parse: simple alignment" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  ACGA
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expectEqual(@as(usize, 2), msa.nseq());
    try testing.expectEqual(@as(usize, 4), msa.alen);
    try testing.expectEqualStrings("seq1", msa.names[0]);
    try testing.expectEqualStrings("seq2", msa.names[1]);

    const text0 = try abc.textize(allocator, msa.seqs[0]);
    defer allocator.free(text0);
    try testing.expectEqualStrings("ACGT", text0);

    const text1 = try abc.textize(allocator, msa.seqs[1]);
    defer allocator.free(text1);
    try testing.expectEqualStrings("ACGA", text1);
}

test "parse: with GF metadata" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\#=GF ID   myalignment
        \\#=GF AC   PF00001
        \\#=GF DE   Test alignment
        \\
        \\seq1  ACGT
        \\seq2  ACGA
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expectEqual(@as(usize, 2), msa.nseq());
    try testing.expectEqualStrings("myalignment", msa.name.?);
    try testing.expectEqualStrings("PF00001", msa.accession.?);
    try testing.expectEqualStrings("Test alignment", msa.description.?);
}

test "parse: interleaved multi-block alignment" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  ACGA
        \\
        \\seq1  TTTT
        \\seq2  CCCC
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expectEqual(@as(usize, 2), msa.nseq());
    try testing.expectEqual(@as(usize, 8), msa.alen);

    const text0 = try abc.textize(allocator, msa.seqs[0]);
    defer allocator.free(text0);
    try testing.expectEqualStrings("ACGTTTTT", text0);

    const text1 = try abc.textize(allocator, msa.seqs[1]);
    defer allocator.free(text1);
    try testing.expectEqualStrings("ACGACCCC", text1);
}

test "parse: GS markup stored" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\#=GS seq1 DE Human hemoglobin
        \\#=GS seq2 DE Mouse hemoglobin
        \\seq1  ACGT
        \\seq2  TTTT
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expectEqual(@as(usize, 2), msa.nseq());
    try testing.expect(msa.gs_markup != null);
    const gs = msa.gs_markup.?;
    try testing.expectEqual(@as(usize, 2), gs.len);
    try testing.expectEqualStrings("seq1", gs[0].seq_name);
    try testing.expectEqualStrings("DE", gs[0].tag);
    try testing.expectEqualStrings("Human hemoglobin", gs[0].value);
    try testing.expectEqualStrings("seq2", gs[1].seq_name);
    try testing.expectEqualStrings("Mouse hemoglobin", gs[1].value);
}

test "parse: GC markup stored and consensus_ss set" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  TTTT
        \\#=GC SS_cons ....
        \\#=GC RF      xxxx
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expect(msa.gc_markup != null);
    try testing.expect(msa.consensus_ss != null);
    try testing.expect(msa.reference != null);
    try testing.expectEqualStrings("....", msa.consensus_ss.?);
    try testing.expectEqualStrings("xxxx", msa.reference.?);

    const gc = msa.gc_markup.?;
    try testing.expectEqual(@as(usize, 2), gc.len);
    try testing.expectEqualStrings("SS_cons", gc[0].tag);
    try testing.expectEqualStrings("....", gc[0].annotation);
    try testing.expectEqualStrings("RF", gc[1].tag);
    try testing.expectEqualStrings("xxxx", gc[1].annotation);
}

test "parse: GR markup stored" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\#=GR seq1 SS ....
        \\seq2  TTTT
        \\#=GR seq2 SS ....
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expect(msa.gr_markup != null);
    const gr = msa.gr_markup.?;
    try testing.expectEqual(@as(usize, 2), gr.len);
    try testing.expectEqualStrings("seq1", gr[0].seq_name);
    try testing.expectEqualStrings("SS", gr[0].tag);
    try testing.expectEqualStrings("....", gr[0].annotation);
    try testing.expectEqualStrings("seq2", gr[1].seq_name);
    try testing.expectEqualStrings("....", gr[1].annotation);
}

test "parse: interleaved GC and GR annotations concatenated" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\#=GR seq1 SS ..((
        \\#=GC SS_cons ..((
        \\
        \\seq1  TTTT
        \\#=GR seq1 SS ))..
        \\#=GC SS_cons ))..
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    // Block 1: "..((" concatenated with block 2: ")).." = "..(()). ."  (8 chars: ..(()).. )
    try testing.expectEqualStrings("..(())", msa.consensus_ss.?[0..6]);
    try testing.expectEqualStrings("..(())", msa.gr_markup.?[0].annotation[0..6]);
}

test "parse: no markup — all markup fields null" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  TTTT
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expect(msa.gf_markup == null);
    try testing.expect(msa.gc_markup == null);
    try testing.expect(msa.gs_markup == null);
    try testing.expect(msa.gr_markup == null);
    try testing.expect(msa.consensus_ss == null);
    try testing.expect(msa.reference == null);
}

test "parse: missing header returns error" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data = "seq1  ACGT\n//\n";
    try testing.expectError(error.InvalidFormat, parse(allocator, abc, data));
}

test "parse: missing terminator returns error" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data = "# STOCKHOLM 1.0\n\nseq1  ACGT\n";
    try testing.expectError(error.InvalidFormat, parse(allocator, abc, data));
}

test "write: basic output" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGT", "TTTT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa);

    const expected =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  TTTT
        \\//
        \\
    ;
    try testing.expectEqualStrings(expected, buf.items);
}

test "write: with metadata" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"seq1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();
    msa.name = try allocator.dupe(u8, "myaln");
    msa.accession = try allocator.dupe(u8, "PF00001");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa);

    try testing.expect(std.mem.indexOf(u8, buf.items, "#=GF ID   myaln") != null);
    try testing.expect(std.mem.indexOf(u8, buf.items, "#=GF AC   PF00001") != null);
}

test "write: GS, GR, and GC markup written" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGT", "TTTT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    const gs = try allocator.alloc(GsEntry, 1);
    gs[0] = .{
        .seq_name = try allocator.dupe(u8, "seq1"),
        .tag = try allocator.dupe(u8, "DE"),
        .value = try allocator.dupe(u8, "Human"),
    };
    msa.gs_markup = gs;

    const gr = try allocator.alloc(GrEntry, 1);
    gr[0] = .{
        .seq_name = try allocator.dupe(u8, "seq1"),
        .tag = try allocator.dupe(u8, "SS"),
        .annotation = try allocator.dupe(u8, "...."),
    };
    msa.gr_markup = gr;

    const gc = try allocator.alloc(GcEntry, 1);
    gc[0] = .{
        .tag = try allocator.dupe(u8, "SS_cons"),
        .annotation = try allocator.dupe(u8, "...."),
    };
    msa.gc_markup = gc;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa);

    try testing.expect(std.mem.indexOf(u8, buf.items, "#=GS seq1 DE Human") != null);
    try testing.expect(std.mem.indexOf(u8, buf.items, "#=GR seq1 SS ....") != null);
    try testing.expect(std.mem.indexOf(u8, buf.items, "#=GC SS_cons ....") != null);
}

test "round-trip: write then parse" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "alpha", "beta" };
    const seqs = [_][]const u8{ "ACGT-A", "ACGA-T" };

    var original = try Msa.fromText(allocator, abc, &names, &seqs);
    defer original.deinit();
    original.name = try allocator.dupe(u8, "testaln");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), original);

    var restored = try parse(allocator, abc, buf.items);
    defer restored.deinit();

    try testing.expectEqual(original.nseq(), restored.nseq());
    try testing.expectEqual(original.alen, restored.alen);
    try testing.expectEqualStrings(original.names[0], restored.names[0]);
    try testing.expectEqualStrings(original.names[1], restored.names[1]);
    try testing.expectEqualSlices(u8, original.seqs[0], restored.seqs[0]);
    try testing.expectEqualSlices(u8, original.seqs[1], restored.seqs[1]);
    try testing.expectEqualStrings(original.name.?, restored.name.?);
}

test "round-trip: full markup preserved" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGT", "TTTT" };

    var original = try Msa.fromText(allocator, abc, &names, &seqs);
    defer original.deinit();

    original.name = try allocator.dupe(u8, "myaln");
    original.accession = try allocator.dupe(u8, "PF00042");
    original.description = try allocator.dupe(u8, "Test family");

    // gf_markup includes ID/AC/DE plus an extra tag.
    const gf = try allocator.alloc(GfEntry, 4);
    gf[0] = .{ .tag = try allocator.dupe(u8, "ID"), .value = try allocator.dupe(u8, "myaln") };
    gf[1] = .{ .tag = try allocator.dupe(u8, "AC"), .value = try allocator.dupe(u8, "PF00042") };
    gf[2] = .{ .tag = try allocator.dupe(u8, "DE"), .value = try allocator.dupe(u8, "Test family") };
    gf[3] = .{ .tag = try allocator.dupe(u8, "TP"), .value = try allocator.dupe(u8, "Family") };
    original.gf_markup = gf;

    const gs = try allocator.alloc(GsEntry, 1);
    gs[0] = .{
        .seq_name = try allocator.dupe(u8, "seq1"),
        .tag = try allocator.dupe(u8, "DE"),
        .value = try allocator.dupe(u8, "Human seq"),
    };
    original.gs_markup = gs;

    const gr = try allocator.alloc(GrEntry, 1);
    gr[0] = .{
        .seq_name = try allocator.dupe(u8, "seq1"),
        .tag = try allocator.dupe(u8, "SS"),
        .annotation = try allocator.dupe(u8, "...."),
    };
    original.gr_markup = gr;

    const gc = try allocator.alloc(GcEntry, 2);
    gc[0] = .{ .tag = try allocator.dupe(u8, "SS_cons"), .annotation = try allocator.dupe(u8, "....") };
    gc[1] = .{ .tag = try allocator.dupe(u8, "RF"), .annotation = try allocator.dupe(u8, "xxxx") };
    original.gc_markup = gc;
    original.consensus_ss = try allocator.dupe(u8, "....");
    original.reference = try allocator.dupe(u8, "xxxx");

    // Write.
    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);
    try write(buf.writer(allocator).any(), original);

    // Parse back.
    var restored = try parse(allocator, abc, buf.items);
    defer restored.deinit();

    // Sequence data.
    try testing.expectEqual(original.nseq(), restored.nseq());
    try testing.expectEqual(original.alen, restored.alen);
    try testing.expectEqualStrings("seq1", restored.names[0]);
    try testing.expectEqualStrings("seq2", restored.names[1]);

    // Dedicated fields.
    try testing.expectEqualStrings("myaln", restored.name.?);
    try testing.expectEqualStrings("PF00042", restored.accession.?);
    try testing.expectEqualStrings("Test family", restored.description.?);

    // gf_markup.
    try testing.expect(restored.gf_markup != null);
    try testing.expectEqual(@as(usize, 4), restored.gf_markup.?.len);
    try testing.expectEqualStrings("TP", restored.gf_markup.?[3].tag);
    try testing.expectEqualStrings("Family", restored.gf_markup.?[3].value);

    // gs_markup.
    try testing.expect(restored.gs_markup != null);
    try testing.expectEqual(@as(usize, 1), restored.gs_markup.?.len);
    try testing.expectEqualStrings("seq1", restored.gs_markup.?[0].seq_name);
    try testing.expectEqualStrings("Human seq", restored.gs_markup.?[0].value);

    // gr_markup.
    try testing.expect(restored.gr_markup != null);
    try testing.expectEqual(@as(usize, 1), restored.gr_markup.?.len);
    try testing.expectEqualStrings("seq1", restored.gr_markup.?[0].seq_name);
    try testing.expectEqualStrings("SS", restored.gr_markup.?[0].tag);
    try testing.expectEqualStrings("....", restored.gr_markup.?[0].annotation);

    // gc_markup.
    try testing.expect(restored.gc_markup != null);
    try testing.expectEqual(@as(usize, 2), restored.gc_markup.?.len);
    try testing.expectEqualStrings("SS_cons", restored.gc_markup.?[0].tag);
    try testing.expectEqualStrings("....", restored.gc_markup.?[0].annotation);
    try testing.expectEqualStrings("RF", restored.gc_markup.?[1].tag);
    try testing.expectEqualStrings("xxxx", restored.gc_markup.?[1].annotation);

    // Convenience fields.
    try testing.expectEqualStrings("....", restored.consensus_ss.?);
    try testing.expectEqualStrings("xxxx", restored.reference.?);
}

test "parseAll: multiple MSA records" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  TTTT
        \\//
        \\# STOCKHOLM 1.0
        \\
        \\alpha  GGGG
        \\beta   CCCC
        \\//
    ;

    const msas = try parseAll(allocator, abc, data);
    defer {
        for (msas) |*m| m.deinit();
        allocator.free(msas);
    }

    try testing.expectEqual(@as(usize, 2), msas.len);
    try testing.expectEqual(@as(usize, 2), msas[0].nseq());
    try testing.expectEqualStrings("seq1", msas[0].names[0]);
    try testing.expectEqual(@as(usize, 2), msas[1].nseq());
    try testing.expectEqualStrings("alpha", msas[1].names[0]);
    try testing.expectEqualStrings("beta", msas[1].names[1]);
}

test "parseAll: single MSA record" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\//
    ;

    const msas = try parseAll(allocator, abc, data);
    defer {
        for (msas) |*m| m.deinit();
        allocator.free(msas);
    }

    try testing.expectEqual(@as(usize, 1), msas.len);
    try testing.expectEqual(@as(usize, 1), msas[0].nseq());
}

test "parseAll: empty input returns empty slice" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const msas = try parseAll(allocator, abc, "");
    defer allocator.free(msas);

    try testing.expectEqual(@as(usize, 0), msas.len);
}

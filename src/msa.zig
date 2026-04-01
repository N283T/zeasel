// Multiple Sequence Alignment (MSA) in digital encoding.

const std = @import("std");
const Allocator = std.mem.Allocator;
const alphabet_mod = @import("alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;
const Sequence = @import("sequence.zig").Sequence;
const wuss = @import("wuss.zig");

/// A single #=GF line: per-file annotation tag and value.
pub const GfEntry = struct { tag: []const u8, value: []const u8 };

/// A single #=GC line: per-column annotation tag and annotation string.
pub const GcEntry = struct { tag: []const u8, annotation: []const u8 };

/// A single #=GS line: per-sequence annotation (seq name, tag, value).
pub const GsEntry = struct { seq_name: []const u8, tag: []const u8, value: []const u8 };

/// A single #=GR line: per-residue annotation (seq name, tag, annotation string).
pub const GrEntry = struct { seq_name: []const u8, tag: []const u8, annotation: []const u8 };

pub const Msa = struct {
    names: [][]const u8, // Sequence names [nseq], each is allocator-owned
    seqs: [][]u8, // Digital sequences [nseq][alen], each row is allocator-owned
    alen: usize, // Alignment length (columns)
    abc: *const Alphabet,
    allocator: Allocator,

    // Optional metadata (all allocator-owned if non-null)
    name: ?[]const u8 = null,
    accession: ?[]const u8 = null,
    description: ?[]const u8 = null,
    author: ?[]const u8 = null,
    weights: ?[]f64 = null,

    // Per-sequence optional metadata [nseq], each entry is optional.
    seq_accessions: ?[]?[]const u8 = null,
    seq_descriptions: ?[]?[]const u8 = null,

    // Per-sequence annotations indexed by sequence index (allocator-owned if non-null)
    ss: ?[]?[]const u8 = null, // per-sequence secondary structure (#=GR SS)
    sa: ?[]?[]const u8 = null, // per-sequence surface accessibility (#=GR SA)
    pp: ?[]?[]const u8 = null, // per-sequence posterior probability (#=GR PP)

    // Per-column consensus annotations (allocator-owned if non-null)
    consensus_ss: ?[]const u8 = null, // #=GC SS_cons
    reference: ?[]const u8 = null, // #=GC RF
    sa_cons: ?[]const u8 = null, // #=GC SA_cons
    pp_cons: ?[]const u8 = null, // #=GC PP_cons
    mm: ?[]const u8 = null, // #=GC MM

    // Generic markup storage for round-trip fidelity (allocator-owned if non-null)
    gf_markup: ?[]GfEntry = null, // all #=GF lines
    gc_markup: ?[]GcEntry = null, // all #=GC lines (including SS_cons/RF)
    gs_markup: ?[]GsEntry = null, // all #=GS lines
    gr_markup: ?[]GrEntry = null, // all #=GR lines

    pub fn nseq(self: Msa) usize {
        return self.names.len;
    }

    /// Create MSA from text sequences (all must be same length).
    /// This is the primary constructor — the Stockholm parser will use this.
    pub fn fromText(
        allocator: Allocator,
        abc: *const Alphabet,
        names: []const []const u8,
        text_seqs: []const []const u8,
    ) !Msa {
        if (names.len != text_seqs.len) return error.InvalidInput;
        if (names.len == 0) return error.InvalidInput;

        const alen = text_seqs[0].len;
        for (text_seqs) |s| {
            if (s.len != alen) return error.InvalidInput;
        }

        const n = names.len;

        var owned_names = try allocator.alloc([]const u8, n);
        var names_done: usize = 0;
        errdefer {
            for (0..names_done) |i| allocator.free(owned_names[i]);
            allocator.free(owned_names);
        }

        var owned_seqs = try allocator.alloc([]u8, n);
        var seqs_done: usize = 0;
        errdefer {
            for (0..seqs_done) |i| allocator.free(owned_seqs[i]);
            allocator.free(owned_seqs);
        }

        for (0..n) |i| {
            owned_names[i] = try allocator.dupe(u8, names[i]);
            names_done += 1;
            owned_seqs[i] = try abc.digitize(allocator, text_seqs[i]);
            seqs_done += 1;
        }

        return Msa{
            .names = owned_names,
            .seqs = owned_seqs,
            .alen = alen,
            .abc = abc,
            .allocator = allocator,
        };
    }

    // --- Metadata setters ---
    // These replace the current value, freeing the old one.

    pub fn setName(self: *Msa, s: []const u8) !void {
        if (self.name) |old| self.allocator.free(old);
        self.name = try self.allocator.dupe(u8, s);
    }

    pub fn setDescription(self: *Msa, s: []const u8) !void {
        if (self.description) |old| self.allocator.free(old);
        self.description = try self.allocator.dupe(u8, s);
    }

    pub fn setAccession(self: *Msa, s: []const u8) !void {
        if (self.accession) |old| self.allocator.free(old);
        self.accession = try self.allocator.dupe(u8, s);
    }

    pub fn setAuthor(self: *Msa, s: []const u8) !void {
        if (self.author) |old| self.allocator.free(old);
        self.author = try self.allocator.dupe(u8, s);
    }

    /// Set per-sequence accession. Allocates the accessions array on first use.
    pub fn setSeqAccession(self: *Msa, idx: usize, s: []const u8) !void {
        if (idx >= self.names.len) return error.OutOfBounds;
        if (self.seq_accessions == null) {
            self.seq_accessions = try self.allocator.alloc(?[]const u8, self.names.len);
            for (self.seq_accessions.?) |*slot| slot.* = null;
        }
        if (self.seq_accessions.?[idx]) |old| self.allocator.free(old);
        self.seq_accessions.?[idx] = try self.allocator.dupe(u8, s);
    }

    /// Set per-sequence description. Allocates the descriptions array on first use.
    pub fn setSeqDescription(self: *Msa, idx: usize, s: []const u8) !void {
        if (idx >= self.names.len) return error.OutOfBounds;
        if (self.seq_descriptions == null) {
            self.seq_descriptions = try self.allocator.alloc(?[]const u8, self.names.len);
            for (self.seq_descriptions.?) |*slot| slot.* = null;
        }
        if (self.seq_descriptions.?[idx]) |old| self.allocator.free(old);
        self.seq_descriptions.?[idx] = try self.allocator.dupe(u8, s);
    }

    /// Set per-sequence secondary structure annotation. Allocates the ss array on first use.
    pub fn setSeqSS(self: *Msa, idx: usize, s: []const u8) !void {
        if (idx >= self.names.len) return error.OutOfBounds;
        if (self.ss == null) {
            self.ss = try self.allocator.alloc(?[]const u8, self.names.len);
            for (self.ss.?) |*slot| slot.* = null;
        }
        if (self.ss.?[idx]) |old| self.allocator.free(old);
        self.ss.?[idx] = try self.allocator.dupe(u8, s);
    }

    /// Set per-sequence surface accessibility annotation. Allocates the sa array on first use.
    pub fn setSeqSA(self: *Msa, idx: usize, s: []const u8) !void {
        if (idx >= self.names.len) return error.OutOfBounds;
        if (self.sa == null) {
            self.sa = try self.allocator.alloc(?[]const u8, self.names.len);
            for (self.sa.?) |*slot| slot.* = null;
        }
        if (self.sa.?[idx]) |old| self.allocator.free(old);
        self.sa.?[idx] = try self.allocator.dupe(u8, s);
    }

    /// Set per-sequence posterior probability annotation. Allocates the pp array on first use.
    pub fn setSeqPP(self: *Msa, idx: usize, s: []const u8) !void {
        if (idx >= self.names.len) return error.OutOfBounds;
        if (self.pp == null) {
            self.pp = try self.allocator.alloc(?[]const u8, self.names.len);
            for (self.pp.?) |*slot| slot.* = null;
        }
        if (self.pp.?[idx]) |old| self.allocator.free(old);
        self.pp.?[idx] = try self.allocator.dupe(u8, s);
    }

    /// Set all sequence weights to 1.0 (uniform).
    pub fn setDefaultWeights(self: *Msa) !void {
        if (self.weights) |old| self.allocator.free(old);
        const w = try self.allocator.alloc(f64, self.names.len);
        @memset(w, 1.0);
        self.weights = w;
    }

    // --- Annotation management ---

    /// Append a #=GF tag/value pair to the GF markup list.
    pub fn addGF(self: *Msa, tag: []const u8, value: []const u8) !void {
        const new_tag = try self.allocator.dupe(u8, tag);
        errdefer self.allocator.free(new_tag);
        const new_value = try self.allocator.dupe(u8, value);
        errdefer self.allocator.free(new_value);

        const entry = GfEntry{ .tag = new_tag, .value = new_value };
        if (self.gf_markup) |old| {
            const new_list = try self.allocator.alloc(GfEntry, old.len + 1);
            @memcpy(new_list[0..old.len], old);
            new_list[old.len] = entry;
            self.allocator.free(old);
            self.gf_markup = new_list;
        } else {
            const new_list = try self.allocator.alloc(GfEntry, 1);
            new_list[0] = entry;
            self.gf_markup = new_list;
        }
    }

    /// Append a #=GS tag/value pair for a specific sequence.
    pub fn addGS(self: *Msa, seq_name: []const u8, tag: []const u8, value: []const u8) !void {
        const new_seq_name = try self.allocator.dupe(u8, seq_name);
        errdefer self.allocator.free(new_seq_name);
        const new_tag = try self.allocator.dupe(u8, tag);
        errdefer self.allocator.free(new_tag);
        const new_value = try self.allocator.dupe(u8, value);
        errdefer self.allocator.free(new_value);

        const entry = GsEntry{ .seq_name = new_seq_name, .tag = new_tag, .value = new_value };
        if (self.gs_markup) |old| {
            const new_list = try self.allocator.alloc(GsEntry, old.len + 1);
            @memcpy(new_list[0..old.len], old);
            new_list[old.len] = entry;
            self.allocator.free(old);
            self.gs_markup = new_list;
        } else {
            const new_list = try self.allocator.alloc(GsEntry, 1);
            new_list[0] = entry;
            self.gs_markup = new_list;
        }
    }

    /// Append a #=GC per-column annotation.
    pub fn appendGC(self: *Msa, tag: []const u8, annotation: []const u8) !void {
        const new_tag = try self.allocator.dupe(u8, tag);
        errdefer self.allocator.free(new_tag);
        const new_annotation = try self.allocator.dupe(u8, annotation);
        errdefer self.allocator.free(new_annotation);

        const entry = GcEntry{ .tag = new_tag, .annotation = new_annotation };
        if (self.gc_markup) |old| {
            const new_list = try self.allocator.alloc(GcEntry, old.len + 1);
            @memcpy(new_list[0..old.len], old);
            new_list[old.len] = entry;
            self.allocator.free(old);
            self.gc_markup = new_list;
        } else {
            const new_list = try self.allocator.alloc(GcEntry, 1);
            new_list[0] = entry;
            self.gc_markup = new_list;
        }
    }

    /// Append a #=GR per-residue annotation for a specific sequence.
    pub fn appendGR(self: *Msa, seq_name: []const u8, tag: []const u8, annotation: []const u8) !void {
        const new_seq_name = try self.allocator.dupe(u8, seq_name);
        errdefer self.allocator.free(new_seq_name);
        const new_tag = try self.allocator.dupe(u8, tag);
        errdefer self.allocator.free(new_tag);
        const new_annotation = try self.allocator.dupe(u8, annotation);
        errdefer self.allocator.free(new_annotation);

        const entry = GrEntry{ .seq_name = new_seq_name, .tag = new_tag, .annotation = new_annotation };
        if (self.gr_markup) |old| {
            const new_list = try self.allocator.alloc(GrEntry, old.len + 1);
            @memcpy(new_list[0..old.len], old);
            new_list[old.len] = entry;
            self.allocator.free(old);
            self.gr_markup = new_list;
        } else {
            const new_list = try self.allocator.alloc(GrEntry, 1);
            new_list[0] = entry;
            self.gr_markup = new_list;
        }
    }

    // --- Validation and utilities ---

    /// Check that all sequence names in the MSA are unique.
    pub fn checkUniqueNames(self: Msa) !bool {
        var seen = std.StringHashMap(void).init(self.allocator);
        defer seen.deinit();
        for (self.names) |name| {
            const result = try seen.getOrPut(name);
            if (result.found_existing) return false;
        }
        return true;
    }

    /// Build a name-to-index hash map for O(1) sequence lookup.
    /// Caller owns the returned map and must call deinit().
    pub fn hash(self: Msa) !std.StringHashMap(usize) {
        var map = std.StringHashMap(usize).init(self.allocator);
        errdefer map.deinit();
        for (self.names, 0..) |name, i| {
            try map.put(name, i);
        }
        return map;
    }

    /// Compute a 32-bit checksum over alignment data using Jenkins one-at-a-time hash.
    /// Matches Easel's esl_msa_Checksum().
    pub fn checksum(self: Msa) u32 {
        var val: u32 = 0;
        for (self.seqs) |seq| {
            for (seq) |code| {
                val +%= @as(u32, code);
                val +%= val << 10;
                val ^= val >> 6;
            }
        }
        // Finalization
        val +%= val << 3;
        val ^= val >> 11;
        val +%= val << 15;
        return val;
    }

    /// Validate MSA internal consistency.
    /// Returns true if valid, false if problems found.
    pub fn validate(self: Msa) bool {
        if (self.names.len != self.seqs.len) return false;
        if (self.names.len == 0) return false;
        for (self.seqs) |seq| {
            if (seq.len != self.alen) return false;
        }
        if (self.weights) |w| {
            if (w.len != self.names.len) return false;
        }
        if (self.consensus_ss) |ss| {
            if (ss.len != self.alen) return false;
        }
        if (self.reference) |rf| {
            if (rf.len != self.alen) return false;
        }
        if (self.seq_accessions) |accs| {
            if (accs.len != self.names.len) return false;
        }
        if (self.seq_descriptions) |descs| {
            if (descs.len != self.names.len) return false;
        }
        if (self.ss) |arr| {
            if (arr.len != self.names.len) return false;
        }
        if (self.sa) |arr| {
            if (arr.len != self.names.len) return false;
        }
        if (self.pp) |arr| {
            if (arr.len != self.names.len) return false;
        }
        if (self.sa_cons) |v| {
            if (v.len != self.alen) return false;
        }
        if (self.pp_cons) |v| {
            if (v.len != self.alen) return false;
        }
        if (self.mm) |v| {
            if (v.len != self.alen) return false;
        }
        return true;
    }

    /// Compare two MSAs for equality (alignment data and names only).
    pub fn compare(self: Msa, other: Msa) bool {
        if (self.nseq() != other.nseq()) return false;
        if (self.alen != other.alen) return false;
        for (self.names, other.names) |a, b| {
            if (!std.mem.eql(u8, a, b)) return false;
        }
        for (self.seqs, other.seqs) |a, b| {
            if (!std.mem.eql(u8, a, b)) return false;
        }
        return true;
    }

    // --- Manipulation (return new Msa) ---

    /// Extract a single sequence from the alignment, removing gap columns.
    pub fn extractSeq(self: Msa, idx: usize) !Sequence {
        if (idx >= self.names.len) return error.OutOfBounds;

        // Count non-gap residues
        var count: usize = 0;
        for (self.seqs[idx]) |code| {
            if (!self.abc.isGap(code)) count += 1;
        }

        // Build ungapped sequence
        const dsq = try self.allocator.alloc(u8, count);
        errdefer self.allocator.free(dsq);
        var j: usize = 0;
        for (self.seqs[idx]) |code| {
            if (!self.abc.isGap(code)) {
                dsq[j] = code;
                j += 1;
            }
        }

        const name_copy = try self.allocator.dupe(u8, self.names[idx]);
        errdefer self.allocator.free(name_copy);

        return Sequence{
            .name = name_copy,
            .accession = null,
            .description = null,
            .taxonomy_id = null,
            .dsq = dsq,
            .secondary_structure = null,
            .source = null,
            .abc = self.abc,
            .allocator = self.allocator,
        };
    }

    /// Select a subset of columns by boolean mask.
    /// mask.len must equal self.alen. Returns a new Msa with only true columns.
    pub fn selectColumns(self: Msa, mask: []const bool) !Msa {
        if (mask.len != self.alen) return error.InvalidInput;

        // Count selected columns
        var new_alen: usize = 0;
        for (mask) |m| {
            if (m) new_alen += 1;
        }

        const n = self.names.len;

        var new_names = try self.allocator.alloc([]const u8, n);
        var names_done: usize = 0;
        errdefer {
            for (0..names_done) |i| self.allocator.free(new_names[i]);
            self.allocator.free(new_names);
        }

        var new_seqs = try self.allocator.alloc([]u8, n);
        var seqs_done: usize = 0;
        errdefer {
            for (0..seqs_done) |i| self.allocator.free(new_seqs[i]);
            self.allocator.free(new_seqs);
        }

        for (0..n) |i| {
            new_names[i] = try self.allocator.dupe(u8, self.names[i]);
            names_done += 1;

            new_seqs[i] = try self.allocator.alloc(u8, new_alen);
            seqs_done += 1;

            var col: usize = 0;
            for (0..self.alen) |c| {
                if (mask[c]) {
                    new_seqs[i][col] = self.seqs[i][c];
                    col += 1;
                }
            }
        }

        // Filter per-column annotations, fixing broken base pairs in SS first.
        var new_ss: ?[]const u8 = null;
        errdefer if (new_ss) |s| self.allocator.free(s);
        if (self.consensus_ss) |ss| {
            // Fix broken base pairs at original length, then select columns.
            const fixed = try fixBrokenBasepairs(self.allocator, ss, mask);
            defer self.allocator.free(fixed);
            const buf = try self.allocator.alloc(u8, new_alen);
            var col: usize = 0;
            for (0..self.alen) |c| {
                if (mask[c]) {
                    buf[col] = fixed[c];
                    col += 1;
                }
            }
            new_ss = buf;
        }

        var new_rf: ?[]const u8 = null;
        errdefer if (new_rf) |s| self.allocator.free(s);
        if (self.reference) |rf| {
            const buf = try self.allocator.alloc(u8, new_alen);
            var col: usize = 0;
            for (0..self.alen) |c| {
                if (mask[c]) {
                    buf[col] = rf[c];
                    col += 1;
                }
            }
            new_rf = buf;
        }

        // Filter additional per-column consensus annotations
        var new_sa_cons: ?[]const u8 = null;
        errdefer if (new_sa_cons) |s| self.allocator.free(s);
        if (self.sa_cons) |v| {
            const buf = try self.allocator.alloc(u8, new_alen);
            var col: usize = 0;
            for (0..self.alen) |c| {
                if (mask[c]) {
                    buf[col] = v[c];
                    col += 1;
                }
            }
            new_sa_cons = buf;
        }

        var new_pp_cons: ?[]const u8 = null;
        errdefer if (new_pp_cons) |s| self.allocator.free(s);
        if (self.pp_cons) |v| {
            const buf = try self.allocator.alloc(u8, new_alen);
            var col: usize = 0;
            for (0..self.alen) |c| {
                if (mask[c]) {
                    buf[col] = v[c];
                    col += 1;
                }
            }
            new_pp_cons = buf;
        }

        var new_mm: ?[]const u8 = null;
        errdefer if (new_mm) |s| self.allocator.free(s);
        if (self.mm) |v| {
            const buf = try self.allocator.alloc(u8, new_alen);
            var col: usize = 0;
            for (0..self.alen) |c| {
                if (mask[c]) {
                    buf[col] = v[c];
                    col += 1;
                }
            }
            new_mm = buf;
        }

        // Filter per-sequence annotations (ss/sa/pp)
        var new_per_ss: ?[]?[]const u8 = null;
        var per_ss_done: usize = 0;
        errdefer {
            if (new_per_ss) |arr| {
                for (0..per_ss_done) |i| {
                    if (arr[i]) |s| self.allocator.free(s);
                }
                self.allocator.free(arr);
            }
        }
        if (self.ss) |arr| {
            new_per_ss = try self.allocator.alloc(?[]const u8, n);
            for (0..n) |i| {
                if (arr[i]) |ann| {
                    // Fix broken base pairs at original length, then select columns.
                    const fixed = try fixBrokenBasepairs(self.allocator, ann, mask);
                    defer self.allocator.free(fixed);
                    const buf = try self.allocator.alloc(u8, new_alen);
                    var col: usize = 0;
                    for (0..self.alen) |c| {
                        if (mask[c]) {
                            buf[col] = fixed[c];
                            col += 1;
                        }
                    }
                    new_per_ss.?[i] = buf;
                } else {
                    new_per_ss.?[i] = null;
                }
                per_ss_done += 1;
            }
        }

        var new_per_sa: ?[]?[]const u8 = null;
        var per_sa_done: usize = 0;
        errdefer {
            if (new_per_sa) |arr| {
                for (0..per_sa_done) |i| {
                    if (arr[i]) |s| self.allocator.free(s);
                }
                self.allocator.free(arr);
            }
        }
        if (self.sa) |arr| {
            new_per_sa = try self.allocator.alloc(?[]const u8, n);
            for (0..n) |i| {
                if (arr[i]) |ann| {
                    const buf = try self.allocator.alloc(u8, new_alen);
                    var col: usize = 0;
                    for (0..self.alen) |c| {
                        if (mask[c]) {
                            buf[col] = ann[c];
                            col += 1;
                        }
                    }
                    new_per_sa.?[i] = buf;
                } else {
                    new_per_sa.?[i] = null;
                }
                per_sa_done += 1;
            }
        }

        var new_per_pp: ?[]?[]const u8 = null;
        var per_pp_done: usize = 0;
        errdefer {
            if (new_per_pp) |arr| {
                for (0..per_pp_done) |i| {
                    if (arr[i]) |s| self.allocator.free(s);
                }
                self.allocator.free(arr);
            }
        }
        if (self.pp) |arr| {
            new_per_pp = try self.allocator.alloc(?[]const u8, n);
            for (0..n) |i| {
                if (arr[i]) |ann| {
                    const buf = try self.allocator.alloc(u8, new_alen);
                    var col: usize = 0;
                    for (0..self.alen) |c| {
                        if (mask[c]) {
                            buf[col] = ann[c];
                            col += 1;
                        }
                    }
                    new_per_pp.?[i] = buf;
                } else {
                    new_per_pp.?[i] = null;
                }
                per_pp_done += 1;
            }
        }

        // Copy weights
        var new_weights: ?[]f64 = null;
        errdefer if (new_weights) |w| self.allocator.free(w);
        if (self.weights) |w| {
            new_weights = try self.allocator.dupe(f64, w);
        }

        return Msa{
            .names = new_names,
            .seqs = new_seqs,
            .alen = new_alen,
            .abc = self.abc,
            .allocator = self.allocator,
            .name = if (self.name) |v| try self.allocator.dupe(u8, v) else null,
            .accession = if (self.accession) |v| try self.allocator.dupe(u8, v) else null,
            .description = if (self.description) |v| try self.allocator.dupe(u8, v) else null,
            .author = if (self.author) |v| try self.allocator.dupe(u8, v) else null,
            .weights = new_weights,
            .consensus_ss = new_ss,
            .reference = new_rf,
            .sa_cons = new_sa_cons,
            .pp_cons = new_pp_cons,
            .mm = new_mm,
            .ss = new_per_ss,
            .sa = new_per_sa,
            .pp = new_per_pp,
        };
    }

    /// Select a subset of sequences by boolean mask.
    /// useme.len must equal nseq. Returns a new Msa with only selected sequences.
    pub fn sequenceSubset(self: Msa, useme: []const bool) !Msa {
        if (useme.len != self.names.len) return error.InvalidInput;

        var nnew: usize = 0;
        for (useme) |u| {
            if (u) nnew += 1;
        }
        if (nnew == 0) return error.InvalidInput;

        var new_names = try self.allocator.alloc([]const u8, nnew);
        var names_done: usize = 0;
        errdefer {
            for (0..names_done) |i| self.allocator.free(new_names[i]);
            self.allocator.free(new_names);
        }

        var new_seqs = try self.allocator.alloc([]u8, nnew);
        var seqs_done: usize = 0;
        errdefer {
            for (0..seqs_done) |i| self.allocator.free(new_seqs[i]);
            self.allocator.free(new_seqs);
        }

        var new_weights: ?[]f64 = null;
        errdefer if (new_weights) |w| self.allocator.free(w);
        if (self.weights != null) {
            new_weights = try self.allocator.alloc(f64, nnew);
        }

        var new_seq_acc: ?[]?[]const u8 = null;
        var seq_acc_done: usize = 0;
        errdefer {
            if (new_seq_acc) |accs| {
                for (0..seq_acc_done) |i| {
                    if (accs[i]) |a| self.allocator.free(a);
                }
                self.allocator.free(accs);
            }
        }
        if (self.seq_accessions != null) {
            new_seq_acc = try self.allocator.alloc(?[]const u8, nnew);
        }

        var new_seq_desc: ?[]?[]const u8 = null;
        var seq_desc_done: usize = 0;
        errdefer {
            if (new_seq_desc) |descs| {
                for (0..seq_desc_done) |i| {
                    if (descs[i]) |d| self.allocator.free(d);
                }
                self.allocator.free(descs);
            }
        }
        if (self.seq_descriptions != null) {
            new_seq_desc = try self.allocator.alloc(?[]const u8, nnew);
        }

        // Per-sequence SS/SA/PP
        var new_per_ss: ?[]?[]const u8 = null;
        var per_ss_done: usize = 0;
        errdefer {
            if (new_per_ss) |arr| {
                for (0..per_ss_done) |i| {
                    if (arr[i]) |s| self.allocator.free(s);
                }
                self.allocator.free(arr);
            }
        }
        if (self.ss != null) {
            new_per_ss = try self.allocator.alloc(?[]const u8, nnew);
        }

        var new_per_sa: ?[]?[]const u8 = null;
        var per_sa_done: usize = 0;
        errdefer {
            if (new_per_sa) |arr| {
                for (0..per_sa_done) |i| {
                    if (arr[i]) |s| self.allocator.free(s);
                }
                self.allocator.free(arr);
            }
        }
        if (self.sa != null) {
            new_per_sa = try self.allocator.alloc(?[]const u8, nnew);
        }

        var new_per_pp: ?[]?[]const u8 = null;
        var per_pp_done: usize = 0;
        errdefer {
            if (new_per_pp) |arr| {
                for (0..per_pp_done) |i| {
                    if (arr[i]) |s| self.allocator.free(s);
                }
                self.allocator.free(arr);
            }
        }
        if (self.pp != null) {
            new_per_pp = try self.allocator.alloc(?[]const u8, nnew);
        }

        var nidx: usize = 0;
        for (0..self.names.len) |oidx| {
            if (useme[oidx]) {
                new_names[nidx] = try self.allocator.dupe(u8, self.names[oidx]);
                names_done += 1;
                new_seqs[nidx] = try self.allocator.dupe(u8, self.seqs[oidx]);
                seqs_done += 1;
                if (self.weights) |w| {
                    new_weights.?[nidx] = w[oidx];
                }
                if (self.seq_accessions) |accs| {
                    new_seq_acc.?[nidx] = if (accs[oidx]) |a| try self.allocator.dupe(u8, a) else null;
                    seq_acc_done += 1;
                }
                if (self.seq_descriptions) |descs| {
                    new_seq_desc.?[nidx] = if (descs[oidx]) |d| try self.allocator.dupe(u8, d) else null;
                    seq_desc_done += 1;
                }
                if (self.ss) |arr| {
                    new_per_ss.?[nidx] = if (arr[oidx]) |s| try self.allocator.dupe(u8, s) else null;
                    per_ss_done += 1;
                }
                if (self.sa) |arr| {
                    new_per_sa.?[nidx] = if (arr[oidx]) |s| try self.allocator.dupe(u8, s) else null;
                    per_sa_done += 1;
                }
                if (self.pp) |arr| {
                    new_per_pp.?[nidx] = if (arr[oidx]) |s| try self.allocator.dupe(u8, s) else null;
                    per_pp_done += 1;
                }
                nidx += 1;
            }
        }

        return Msa{
            .names = new_names,
            .seqs = new_seqs,
            .alen = self.alen,
            .abc = self.abc,
            .allocator = self.allocator,
            .name = if (self.name) |v| try self.allocator.dupe(u8, v) else null,
            .accession = if (self.accession) |v| try self.allocator.dupe(u8, v) else null,
            .description = if (self.description) |v| try self.allocator.dupe(u8, v) else null,
            .author = if (self.author) |v| try self.allocator.dupe(u8, v) else null,
            .weights = new_weights,
            .consensus_ss = if (self.consensus_ss) |v| try self.allocator.dupe(u8, v) else null,
            .reference = if (self.reference) |v| try self.allocator.dupe(u8, v) else null,
            .sa_cons = if (self.sa_cons) |v| try self.allocator.dupe(u8, v) else null,
            .pp_cons = if (self.pp_cons) |v| try self.allocator.dupe(u8, v) else null,
            .mm = if (self.mm) |v| try self.allocator.dupe(u8, v) else null,
            .seq_accessions = new_seq_acc,
            .seq_descriptions = new_seq_desc,
            .ss = new_per_ss,
            .sa = new_per_sa,
            .pp = new_per_pp,
        };
    }

    /// Remove columns that are all gaps (or all gaps + missing data).
    /// Returns a new Msa with gap-only columns removed.
    /// If consider_rf is true, columns marked in RF annotation are always kept.
    pub fn minimGaps(self: Msa, consider_rf: bool) !Msa {
        const col_mask = try self.allocator.alloc(bool, self.alen);
        defer self.allocator.free(col_mask);

        for (0..self.alen) |col| {
            // If consider_rf and RF is non-gap at this column, keep it
            if (consider_rf) {
                if (self.reference) |rf| {
                    const ch = rf[col];
                    if (ch != '.' and ch != '-' and ch != '_' and ch != '~') {
                        col_mask[col] = true;
                        continue;
                    }
                }
            }
            // Check if any sequence has a non-gap, non-missing residue
            var has_residue = false;
            for (0..self.names.len) |seq| {
                const code = self.seqs[seq][col];
                if (!self.abc.isGap(code) and !self.abc.isMissing(code)) {
                    has_residue = true;
                    break;
                }
            }
            col_mask[col] = has_residue;
        }

        return self.selectColumns(col_mask);
    }

    /// Remove all columns containing any gap character.
    /// Returns a new Msa with no gap-containing columns.
    pub fn noGaps(self: Msa) !Msa {
        const col_mask = try self.allocator.alloc(bool, self.alen);
        defer self.allocator.free(col_mask);

        for (0..self.alen) |col| {
            var has_gap = false;
            for (0..self.names.len) |seq| {
                const code = self.seqs[seq][col];
                if (self.abc.isGap(code) or self.abc.isMissing(code)) {
                    has_gap = true;
                    break;
                }
            }
            col_mask[col] = !has_gap;
        }

        return self.selectColumns(col_mask);
    }

    /// Reverse complement the entire alignment.
    /// Only works for DNA/RNA alphabets. Returns a new Msa.
    pub fn reverseComplement(self: Msa) !Msa {
        const n = self.names.len;

        var new_names = try self.allocator.alloc([]const u8, n);
        var names_done: usize = 0;
        errdefer {
            for (0..names_done) |i| self.allocator.free(new_names[i]);
            self.allocator.free(new_names);
        }

        var new_seqs = try self.allocator.alloc([]u8, n);
        var seqs_done: usize = 0;
        errdefer {
            for (0..seqs_done) |i| self.allocator.free(new_seqs[i]);
            self.allocator.free(new_seqs);
        }

        for (0..n) |i| {
            new_names[i] = try self.allocator.dupe(u8, self.names[i]);
            names_done += 1;
            new_seqs[i] = try self.allocator.dupe(u8, self.seqs[i]);
            seqs_done += 1;
            try self.abc.reverseComplement(new_seqs[i]);
        }

        var result = Msa{
            .names = new_names,
            .seqs = new_seqs,
            .alen = self.alen,
            .abc = self.abc,
            .allocator = self.allocator,
        };

        // Reverse annotation strings (reference, consensus_ss).
        if (self.reference) |ref| {
            const rev = try self.allocator.dupe(u8, ref);
            std.mem.reverse(u8, rev);
            result.reference = rev;
        }
        if (self.consensus_ss) |ss| {
            const rev = try self.allocator.dupe(u8, ss);
            // Reverse bytes and swap WUSS base-pair characters.
            std.mem.reverse(u8, rev);
            for (rev) |*ch| {
                ch.* = switch (ch.*) {
                    '<' => '>',
                    '>' => '<',
                    '(' => ')',
                    ')' => '(',
                    '[' => ']',
                    ']' => '[',
                    '{' => '}',
                    '}' => '{',
                    else => ch.*,
                };
            }
            result.consensus_ss = rev;
        }

        return result;
    }

    /// Convert all degenerate residue codes to the unknown symbol (X/N).
    /// Returns a new Msa.
    pub fn convertDegen2X(self: Msa) !Msa {
        const n = self.names.len;
        const unknown = self.abc.unknownCode();

        var new_names = try self.allocator.alloc([]const u8, n);
        var names_done: usize = 0;
        errdefer {
            for (0..names_done) |i| self.allocator.free(new_names[i]);
            self.allocator.free(new_names);
        }

        var new_seqs = try self.allocator.alloc([]u8, n);
        var seqs_done: usize = 0;
        errdefer {
            for (0..seqs_done) |i| self.allocator.free(new_seqs[i]);
            self.allocator.free(new_seqs);
        }

        for (0..n) |i| {
            new_names[i] = try self.allocator.dupe(u8, self.names[i]);
            names_done += 1;
            new_seqs[i] = try self.allocator.dupe(u8, self.seqs[i]);
            seqs_done += 1;
            for (new_seqs[i]) |*code| {
                if (self.abc.isDegenerate(code.*)) {
                    code.* = unknown;
                }
            }
        }

        return Msa{
            .names = new_names,
            .seqs = new_seqs,
            .alen = self.alen,
            .abc = self.abc,
            .allocator = self.allocator,
        };
    }

    /// Guess alphabet type from the alignment content.
    /// Scans all sequences and uses the alphabet guesser.
    pub fn guessAlphabet(self: Msa) ?alphabet_mod.AlphabetType {
        for (self.seqs) |seq| {
            const text = self.abc.textize(self.allocator, seq) catch continue;
            defer self.allocator.free(text);
            if (alphabet_mod.guessType(text)) |t| {
                if (t == .amino) return .amino;
            }
        }
        // If no amino-only residues found in any sequence, check the first one
        if (self.seqs.len > 0) {
            const text = self.abc.textize(self.allocator, self.seqs[0]) catch return null;
            defer self.allocator.free(text);
            return alphabet_mod.guessType(text);
        }
        return null;
    }

    /// Determine a reasonable RF line marking consensus columns.
    /// Consensus columns have residue occupancy >= symfrac (weighted if weights set).
    /// Returns allocator-owned string of 'x' (consensus) and '.' (insert).
    pub fn reasonableRF(self: Msa, symfrac: f64) ![]u8 {
        const rf = try self.allocator.alloc(u8, self.alen);
        errdefer self.allocator.free(rf);

        for (0..self.alen) |col| {
            var residue_weight: f64 = 0;
            var total_weight: f64 = 0;
            for (0..self.names.len) |seq| {
                const w: f64 = if (self.weights) |wts| wts[seq] else 1.0;
                const code = self.seqs[seq][col];
                if (self.abc.isCanonical(code) or self.abc.isDegenerate(code) or self.abc.isUnknown(code)) {
                    residue_weight += w;
                    total_weight += w;
                } else if (self.abc.isGap(code)) {
                    total_weight += w;
                }
                // missing data is ignored (not counted in total)
            }
            if (total_weight > 0 and residue_weight / total_weight >= symfrac) {
                rf[col] = 'x';
            } else {
                rf[col] = '.';
            }
        }

        return rf;
    }

    /// Free all owned memory.
    pub fn deinit(self: *Msa) void {
        for (self.names) |n| self.allocator.free(n);
        self.allocator.free(self.names);
        for (self.seqs) |s| self.allocator.free(s);
        self.allocator.free(self.seqs);
        if (self.name) |n| self.allocator.free(n);
        if (self.accession) |a| self.allocator.free(a);
        if (self.description) |d| self.allocator.free(d);
        if (self.author) |a| self.allocator.free(a);
        if (self.weights) |w| self.allocator.free(w);
        if (self.consensus_ss) |ss| self.allocator.free(ss);
        if (self.reference) |rf| self.allocator.free(rf);
        if (self.sa_cons) |v| self.allocator.free(v);
        if (self.pp_cons) |v| self.allocator.free(v);
        if (self.mm) |v| self.allocator.free(v);
        if (self.ss) |arr| {
            for (arr) |entry| {
                if (entry) |s| self.allocator.free(s);
            }
            self.allocator.free(arr);
        }
        if (self.sa) |arr| {
            for (arr) |entry| {
                if (entry) |s| self.allocator.free(s);
            }
            self.allocator.free(arr);
        }
        if (self.pp) |arr| {
            for (arr) |entry| {
                if (entry) |s| self.allocator.free(s);
            }
            self.allocator.free(arr);
        }
        if (self.seq_accessions) |accs| {
            for (accs) |acc| {
                if (acc) |a| self.allocator.free(a);
            }
            self.allocator.free(accs);
        }
        if (self.seq_descriptions) |descs| {
            for (descs) |desc| {
                if (desc) |d| self.allocator.free(d);
            }
            self.allocator.free(descs);
        }
        if (self.gf_markup) |entries| {
            for (entries) |e| {
                self.allocator.free(e.tag);
                self.allocator.free(e.value);
            }
            self.allocator.free(entries);
        }
        if (self.gc_markup) |entries| {
            for (entries) |e| {
                self.allocator.free(e.tag);
                self.allocator.free(e.annotation);
            }
            self.allocator.free(entries);
        }
        if (self.gs_markup) |entries| {
            for (entries) |e| {
                self.allocator.free(e.seq_name);
                self.allocator.free(e.tag);
                self.allocator.free(e.value);
            }
            self.allocator.free(entries);
        }
        if (self.gr_markup) |entries| {
            for (entries) |e| {
                self.allocator.free(e.seq_name);
                self.allocator.free(e.tag);
                self.allocator.free(e.annotation);
            }
            self.allocator.free(entries);
        }
    }
};

/// Fix broken base pairs in a WUSS SS string given a column-keep mask.
/// `ss` is the original (full-length) SS annotation. `keep` has the same length.
/// Returns a new string of the same length as `ss` with broken pairs replaced by '.'.
/// A pair is broken if either partner is in a removed column (!keep).
fn fixBrokenBasepairs(allocator: Allocator, ss: []const u8, keep: []const bool) ![]u8 {
    std.debug.assert(ss.len == keep.len);

    const pairs = wuss.parseToPairs(allocator, ss) catch {
        // If SS is already malformed, just return a copy.
        return try allocator.dupe(u8, ss);
    };
    defer allocator.free(pairs);

    const result = try allocator.alloc(u8, ss.len);
    @memcpy(result, ss);

    for (0..ss.len) |i| {
        if (pairs[i] >= 0) {
            const partner: usize = @intCast(pairs[i]);
            if (!keep[i] or !keep[partner]) {
                result[i] = '.';
                result[partner] = '.';
            }
        }
    }

    return result;
}

// --- Tests ---

test "fromText: basic construction" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2", "seq3" };
    const seqs = [_][]const u8{ "ACGT-", "A-GT-", "AC-T-" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 3), msa.nseq());
    try std.testing.expectEqual(@as(usize, 5), msa.alen);
}

test "fromText: mismatched sequence lengths return error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGT", "ACG" }; // different lengths

    try std.testing.expectError(
        error.InvalidInput,
        Msa.fromText(allocator, abc, &names, &seqs),
    );
}

test "fromText: empty input returns error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{};
    const seqs = [_][]const u8{};

    try std.testing.expectError(
        error.InvalidInput,
        Msa.fromText(allocator, abc, &names, &seqs),
    );
}

test "extractSeq: gaps are removed" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"seq1"};
    const seqs = [_][]const u8{"AC-GT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var seq = try msa.extractSeq(0);
    defer seq.deinit();

    const text = try seq.toText();
    defer allocator.free(text);

    try std.testing.expectEqualStrings("ACGT", text);
    try std.testing.expectEqualStrings("seq1", seq.name);
}

test "extractSeq: out of bounds returns error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"seq1"};
    const seqs = [_][]const u8{"ACGT-"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expectError(error.OutOfBounds, msa.extractSeq(1));
}

test "selectColumns: selects correct columns" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGTA", "TGCAT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    // Select columns 0, 2, 4 (true, false, true, false, true)
    const mask = [_]bool{ true, false, true, false, true };
    var sub = try msa.selectColumns(&mask);
    defer sub.deinit();

    try std.testing.expectEqual(@as(usize, 2), sub.nseq());
    try std.testing.expectEqual(@as(usize, 3), sub.alen);

    const text0 = try sub.abc.textize(allocator, sub.seqs[0]);
    defer allocator.free(text0);
    try std.testing.expectEqualStrings("AGA", text0);

    const text1 = try sub.abc.textize(allocator, sub.seqs[1]);
    defer allocator.free(text1);
    try std.testing.expectEqualStrings("TCT", text1);
}

test "selectColumns: mask length mismatch returns error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"seq1"};
    const seqs = [_][]const u8{"ACGTA"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    const mask = [_]bool{ true, false, true }; // wrong length
    try std.testing.expectError(error.InvalidInput, msa.selectColumns(&mask));
}

test "deinit: no memory leaks" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2", "seq3" };
    const seqs = [_][]const u8{ "ACGT-", "A-GT-", "AC-T-" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    msa.deinit();
    // std.testing.allocator detects leaks automatically
}

// --- New tests for added functionality ---

test "setName and setDescription" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setName("test_msa");
    try std.testing.expectEqualStrings("test_msa", msa.name.?);

    try msa.setDescription("a test alignment");
    try std.testing.expectEqualStrings("a test alignment", msa.description.?);

    // Replace existing
    try msa.setName("renamed");
    try std.testing.expectEqualStrings("renamed", msa.name.?);
}

test "setAccession and setAuthor" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setAccession("PF00001");
    try std.testing.expectEqualStrings("PF00001", msa.accession.?);

    try msa.setAuthor("test");
    try std.testing.expectEqualStrings("test", msa.author.?);
}

test "setSeqAccession and setSeqDescription" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setSeqAccession(0, "ACC001");
    try msa.setSeqDescription(1, "second sequence");

    try std.testing.expectEqualStrings("ACC001", msa.seq_accessions.?[0].?);
    try std.testing.expectEqual(@as(?[]const u8, null), msa.seq_accessions.?[1]);
    try std.testing.expectEqual(@as(?[]const u8, null), msa.seq_descriptions.?[0]);
    try std.testing.expectEqualStrings("second sequence", msa.seq_descriptions.?[1].?);
}

test "setSeqAccession: out of bounds" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expectError(error.OutOfBounds, msa.setSeqAccession(5, "x"));
}

test "setDefaultWeights" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AAAA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setDefaultWeights();
    try std.testing.expectEqual(@as(usize, 3), msa.weights.?.len);
    for (msa.weights.?) |w| {
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), w, 1e-9);
    }
}

test "addGF and addGS" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.addGF("AU", "Smith J.");
    try msa.addGF("DE", "Test family");
    try std.testing.expectEqual(@as(usize, 2), msa.gf_markup.?.len);
    try std.testing.expectEqualStrings("AU", msa.gf_markup.?[0].tag);
    try std.testing.expectEqualStrings("Test family", msa.gf_markup.?[1].value);

    try msa.addGS("s1", "AC", "ACC001");
    try std.testing.expectEqual(@as(usize, 1), msa.gs_markup.?.len);
}

test "appendGC and appendGR" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.appendGC("SA_cons", "....");
    try std.testing.expectEqual(@as(usize, 1), msa.gc_markup.?.len);

    try msa.appendGR("s1", "PP", "****");
    try std.testing.expectEqual(@as(usize, 1), msa.gr_markup.?.len);
}

test "checkUniqueNames: unique names" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AAAA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expect(try msa.checkUniqueNames());
}

test "checkUniqueNames: duplicate names" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s1", "s3" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AAAA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expect(!try msa.checkUniqueNames());
}

test "hash: name to index lookup" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "alpha", "beta", "gamma" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AAAA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var map = try msa.hash();
    defer map.deinit();

    try std.testing.expectEqual(@as(usize, 0), map.get("alpha").?);
    try std.testing.expectEqual(@as(usize, 1), map.get("beta").?);
    try std.testing.expectEqual(@as(usize, 2), map.get("gamma").?);
    try std.testing.expectEqual(@as(?usize, null), map.get("delta"));
}

test "checksum: deterministic" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    const c1 = msa.checksum();
    const c2 = msa.checksum();
    try std.testing.expectEqual(c1, c2);
}

test "checksum: different alignments differ" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names1 = [_][]const u8{"s1"};
    const seqs1 = [_][]const u8{"ACGT"};
    var msa1 = try Msa.fromText(allocator, abc, &names1, &seqs1);
    defer msa1.deinit();

    const names2 = [_][]const u8{"s1"};
    const seqs2 = [_][]const u8{"TGCA"};
    var msa2 = try Msa.fromText(allocator, abc, &names2, &seqs2);
    defer msa2.deinit();

    try std.testing.expect(msa1.checksum() != msa2.checksum());
}

test "validate: valid MSA" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expect(msa.validate());
}

test "compare: identical MSAs" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var msa1 = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa1.deinit();

    var msa2 = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa2.deinit();

    try std.testing.expect(msa1.compare(msa2));
}

test "compare: different MSAs" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names1 = [_][]const u8{ "s1", "s2" };
    const seqs1 = [_][]const u8{ "ACGT", "TGCA" };
    var msa1 = try Msa.fromText(allocator, abc, &names1, &seqs1);
    defer msa1.deinit();

    const names2 = [_][]const u8{ "s1", "s3" };
    const seqs2 = [_][]const u8{ "ACGT", "AAAA" };
    var msa2 = try Msa.fromText(allocator, abc, &names2, &seqs2);
    defer msa2.deinit();

    try std.testing.expect(!msa1.compare(msa2));
}

test "sequenceSubset: select 2 of 3" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AAAA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    const useme = [_]bool{ true, false, true };
    var sub = try msa.sequenceSubset(&useme);
    defer sub.deinit();

    try std.testing.expectEqual(@as(usize, 2), sub.nseq());
    try std.testing.expectEqual(@as(usize, 4), sub.alen);
    try std.testing.expectEqualStrings("s1", sub.names[0]);
    try std.testing.expectEqualStrings("s3", sub.names[1]);
}

test "sequenceSubset: no sequences selected returns error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    const useme = [_]bool{ false, false };
    try std.testing.expectError(error.InvalidInput, msa.sequenceSubset(&useme));
}

test "sequenceSubset: preserves weights" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AAAA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setDefaultWeights();
    msa.weights.?[0] = 2.0;
    msa.weights.?[2] = 0.5;

    const useme = [_]bool{ true, false, true };
    var sub = try msa.sequenceSubset(&useme);
    defer sub.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, 2.0), sub.weights.?[0], 1e-9);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), sub.weights.?[1], 1e-9);
}

test "minimGaps: removes all-gap columns" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "A-C-T", "A--GT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var trimmed = try msa.minimGaps(false);
    defer trimmed.deinit();

    // Column 1 has '-' in both -> removed. Column 3: s1='-' s2='G' -> kept
    try std.testing.expectEqual(@as(usize, 4), trimmed.alen);
}

test "minimGaps: keeps RF-marked columns when consider_rf=true" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "A--T", "A--T" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();
    // Mark column 1 as RF even though it's all-gap
    msa.reference = try allocator.dupe(u8, "x.xX");

    var trimmed = try msa.minimGaps(true);
    defer trimmed.deinit();

    // Columns 1,2 are all gap. Column 1 has RF='.' -> removable. Column 2 has RF='x' -> kept
    try std.testing.expectEqual(@as(usize, 3), trimmed.alen);
}

test "noGaps: removes all columns with any gap" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "A-GT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var trimmed = try msa.noGaps();
    defer trimmed.deinit();

    try std.testing.expectEqual(@as(usize, 3), trimmed.alen);
    const text = try trimmed.abc.textize(allocator, trimmed.seqs[0]);
    defer allocator.free(text);
    try std.testing.expectEqualStrings("AGT", text);
}

test "reverseComplement: DNA alignment" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var rc = try msa.reverseComplement();
    defer rc.deinit();

    const text = try rc.abc.textize(allocator, rc.seqs[0]);
    defer allocator.free(text);
    // ACGT is a palindrome
    try std.testing.expectEqualStrings("ACGT", text);
}

test "reverseComplement: non-palindrome" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"AACG"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var rc = try msa.reverseComplement();
    defer rc.deinit();

    const text = try rc.abc.textize(allocator, rc.seqs[0]);
    defer allocator.free(text);
    try std.testing.expectEqualStrings("CGTT", text);
}

test "reverseComplement: amino returns error" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.amino;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACDE"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try std.testing.expectError(error.NoComplement, msa.reverseComplement());
}

test "convertDegen2X: replaces degeneracies" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACRG"}; // R is degenerate

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var converted = try msa.convertDegen2X();
    defer converted.deinit();

    const text = try converted.abc.textize(allocator, converted.seqs[0]);
    defer allocator.free(text);
    try std.testing.expectEqualStrings("ACNG", text); // R -> N (unknown for DNA)
}

test "reasonableRF: basic consensus" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    // 4 seqs: column 0 all A, column 1 half gap, column 2 all gap, column 3 all T
    const names = [_][]const u8{ "s1", "s2", "s3", "s4" };
    const seqs = [_][]const u8{ "A-GT", "A--T", "ACGT", "AC-T" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    const rf = try msa.reasonableRF(0.5);
    defer allocator.free(rf);

    try std.testing.expectEqual(@as(u8, 'x'), rf[0]); // col 0: 4/4 = 1.0 >= 0.5
    try std.testing.expectEqual(@as(u8, 'x'), rf[1]); // col 1: 2/4 = 0.5 >= 0.5
    try std.testing.expectEqual(@as(u8, 'x'), rf[2]); // col 2: 2/4 = 0.5 >= 0.5 (s1=G, s3=G)
    try std.testing.expectEqual(@as(u8, 'x'), rf[3]); // col 3: 4/4 = 1.0 >= 0.5
}

test "selectColumns: preserves consensus_ss and reference" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();
    msa.consensus_ss = try allocator.dupe(u8, "<..>");
    msa.reference = try allocator.dupe(u8, "xx.x");

    const mask = [_]bool{ true, false, false, true };
    var sub = try msa.selectColumns(&mask);
    defer sub.deinit();

    try std.testing.expectEqualStrings("<>", sub.consensus_ss.?);
    try std.testing.expectEqualStrings("xx", sub.reference.?);
}

test "setSeqSS, setSeqSA, setSeqPP" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setSeqSS(0, "<..>");
    try msa.setSeqSA(0, "1234");
    try msa.setSeqPP(1, "****");

    try std.testing.expectEqualStrings("<..>", msa.ss.?[0].?);
    try std.testing.expectEqual(@as(?[]const u8, null), msa.ss.?[1]);
    try std.testing.expectEqualStrings("1234", msa.sa.?[0].?);
    try std.testing.expectEqual(@as(?[]const u8, null), msa.sa.?[1]);
    try std.testing.expectEqual(@as(?[]const u8, null), msa.pp.?[0]);
    try std.testing.expectEqualStrings("****", msa.pp.?[1].?);
}

test "selectColumns: propagates per-sequence SS/SA/PP and consensus sa_cons/pp_cons/mm" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGTA", "TGCAT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setSeqSS(0, "<<.>>");
    try msa.setSeqSA(0, "12345");
    try msa.setSeqPP(1, "*****");
    msa.sa_cons = try allocator.dupe(u8, "abcde");
    msa.pp_cons = try allocator.dupe(u8, "56789");
    msa.mm = try allocator.dupe(u8, "mMmMm");

    // Select columns 0, 2, 4
    const mask = [_]bool{ true, false, true, false, true };
    var sub = try msa.selectColumns(&mask);
    defer sub.deinit();

    try std.testing.expectEqual(@as(usize, 3), sub.alen);
    try std.testing.expectEqualStrings("12345"[0..1] ++ "12345"[2..3] ++ "12345"[4..5], sub.sa.?[0].?);
    try std.testing.expectEqualStrings("*****"[0..1] ++ "*****"[2..3] ++ "*****"[4..5], sub.pp.?[1].?);
    try std.testing.expectEqualStrings("ace", sub.sa_cons.?);
    try std.testing.expectEqualStrings("579", sub.pp_cons.?);
    try std.testing.expectEqualStrings("mmm", sub.mm.?);
}

test "selectColumns: fixes broken base pairs in consensus SS" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGTAA"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();
    // SS: <<..>> means 0-5 paired, 1-4 paired
    msa.consensus_ss = try allocator.dupe(u8, "<<..>>");

    // Remove column 0: pair (0,5) is broken
    const mask = [_]bool{ false, true, true, true, true, true };
    var sub = try msa.selectColumns(&mask);
    defer sub.deinit();

    // Position 0 was removed; its partner (5) should now be '.'
    // Result after fixing and selecting: position 1 paired with 4 still intact
    // Original after fix: ".<..>." -> select columns 1..5 -> "<..>."
    try std.testing.expectEqual(@as(usize, 5), sub.alen);
    try std.testing.expectEqualStrings("<..>.", sub.consensus_ss.?);
}

test "selectColumns: fixes broken base pairs in per-sequence SS" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"s1"};
    const seqs = [_][]const u8{"ACGTAA"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setSeqSS(0, "<<..>>");

    // Remove column 5: pair (0,5) is broken
    const mask = [_]bool{ true, true, true, true, true, false };
    var sub = try msa.selectColumns(&mask);
    defer sub.deinit();

    // After fix: ".<..>." -> select 0..4 -> ".<..>"
    try std.testing.expectEqual(@as(usize, 5), sub.alen);
    try std.testing.expectEqualStrings(".<..>", sub.ss.?[0].?);
}

test "sequenceSubset: propagates per-sequence SS/SA/PP" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2", "s3" };
    const seqs = [_][]const u8{ "ACGT", "TGCA", "AAAA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    try msa.setSeqSS(0, "<..>");
    try msa.setSeqSS(2, "(())");
    try msa.setSeqPP(1, "****");

    const useme = [_]bool{ true, false, true };
    var sub = try msa.sequenceSubset(&useme);
    defer sub.deinit();

    try std.testing.expectEqual(@as(usize, 2), sub.nseq());
    try std.testing.expect(sub.ss != null);
    try std.testing.expectEqualStrings("<..>", sub.ss.?[0].?);
    try std.testing.expectEqualStrings("(())", sub.ss.?[1].?);
    // PP was only on s2 which was removed; s1 and s3 had no PP
    try std.testing.expect(sub.pp != null);
    try std.testing.expectEqual(@as(?[]const u8, null), sub.pp.?[0]);
    try std.testing.expectEqual(@as(?[]const u8, null), sub.pp.?[1]);
}

test "sequenceSubset: propagates consensus sa_cons/pp_cons/mm" {
    const allocator = std.testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();
    msa.sa_cons = try allocator.dupe(u8, "abcd");
    msa.pp_cons = try allocator.dupe(u8, "1234");
    msa.mm = try allocator.dupe(u8, "mmmm");

    const useme = [_]bool{ true, false };
    var sub = try msa.sequenceSubset(&useme);
    defer sub.deinit();

    try std.testing.expectEqualStrings("abcd", sub.sa_cons.?);
    try std.testing.expectEqualStrings("1234", sub.pp_cons.?);
    try std.testing.expectEqualStrings("mmmm", sub.mm.?);
}

// Random biological sequence generation and in-place shuffling.

const std = @import("std");
const Allocator = std.mem.Allocator;
const alphabet_mod = @import("../alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;
const Sequence = @import("../sequence.zig").Sequence;
const Random = @import("random.zig").Random;

/// Generate a random i.i.d. sequence from a residue frequency distribution.
///
/// freq[0..K-1] gives the probability of each canonical residue.
/// The caller owns the returned Sequence and must call deinit() on it.
pub fn randomIid(
    allocator: Allocator,
    rng: *Random,
    abc: *const Alphabet,
    name: []const u8,
    length: usize,
    freq: []const f64,
) !Sequence {
    const dsq = try allocator.alloc(u8, length);
    errdefer allocator.free(dsq);

    for (0..length) |i| {
        dsq[i] = @intCast(rng.choose(freq));
    }

    const name_copy = try allocator.dupe(u8, name);
    errdefer allocator.free(name_copy);

    return Sequence{
        .name = name_copy,
        .accession = null,
        .description = null,
        .taxonomy_id = null,
        .dsq = dsq,
        .secondary_structure = null,
        .source = null,
        .abc = abc,
        .allocator = allocator,
    };
}

/// Generate a random sequence with uniform residue frequencies over K canonical residues.
///
/// The caller owns the returned Sequence and must call deinit() on it.
pub fn randomUniform(
    allocator: Allocator,
    rng: *Random,
    abc: *const Alphabet,
    name: []const u8,
    length: usize,
) !Sequence {
    const k = abc.k;
    const freq_value = 1.0 / @as(f64, @floatFromInt(k));

    // Build uniform frequency array on the stack for small K, heap for large K.
    // Amino acids have K=20 — always fits on the stack.
    var freq_buf: [32]f64 = undefined;
    if (k > freq_buf.len) return error.AlphabetTooLarge;
    for (0..k) |i| freq_buf[i] = freq_value;

    return randomIid(allocator, rng, abc, name, length, freq_buf[0..k]);
}

/// Shuffle a digital sequence in place using the Fisher-Yates algorithm.
/// Preserves mono-residue composition.
pub fn shuffle(rng: *Random, dsq: []u8) void {
    var i = dsq.len;
    while (i > 1) {
        i -= 1;
        const j = rng.uniformInt(@intCast(i + 1));
        const tmp = dsq[i];
        dsq[i] = dsq[j];
        dsq[j] = tmp;
    }
}

/// Shuffle preserving mono-residue composition.
/// Alias for shuffle — included for API clarity matching Easel naming.
pub const shuffleMono = shuffle;

/// Maximum alphabet size supported by the dinucleotide shuffle.
/// Amino acid alphabets have K=20; this leaves headroom.
const MAX_K: usize = 32;

/// Shuffle preserving dinucleotide composition (Altschul-Erickson algorithm).
///
/// Maintains the count of each dinucleotide pair while randomizing the
/// sequence.  `k` is the alphabet size (number of canonical residues).
///
/// For sequences of length 0 or 1, the sequence is unchanged.
/// For length 2, the sequence is unchanged (only one dinucleotide).
///
/// Algorithm (Altschul & Erickson 1985; Easel esl_rsq_XShuffleDP):
///   1. Count directed-edge frequencies E[i][j] from consecutive pairs.
///   2. For each node with outgoing edges, randomly mark one edge as
///      "last" (reserve it so an Eulerian path is guaranteed).
///   3. Traverse from dsq[0]: at each node, pick a random non-last edge;
///      when only the last edge remains, use it.  Append the destination
///      residue and continue until all edges are consumed.
pub fn shuffleDi(rng: *Random, dsq: []u8, k: u8) void {
    const K: usize = @intCast(k);
    if (K > MAX_K) return; // safety guard
    if (dsq.len <= 2) return; // nothing to shuffle

    // --- Step 1: count edge frequencies E[i][j] ---
    var E: [MAX_K][MAX_K]u16 = .{.{0} ** MAX_K} ** MAX_K;
    for (0..dsq.len - 1) |pos| {
        E[dsq[pos]][dsq[pos + 1]] += 1;
    }

    // --- Step 2: for each node, randomly choose one outgoing edge to be "last" ---
    // iE[i] = the destination j of the reserved last edge from node i.
    // We decrement E[i][j] for that edge so it won't be picked during normal traversal.
    var iE: [MAX_K]u8 = .{0} ** MAX_K;
    for (0..K) |i| {
        // Count total outgoing edges from node i.
        var total: u32 = 0;
        for (0..K) |j| total += E[i][j];
        if (total == 0) continue;

        // Pick a random edge uniformly among all outgoing edges.
        var pick: u32 = rng.uniformInt(total);
        for (0..K) |j| {
            if (E[i][j] == 0) continue;
            if (pick < E[i][j]) {
                iE[i] = @intCast(j);
                E[i][j] -= 1;
                break;
            }
            pick -= E[i][j];
        }
    }

    // --- Step 3: traverse Eulerian path ---
    var pos: usize = 0;
    dsq[pos] = dsq[0]; // first residue stays the same
    var x: usize = dsq[0];

    while (pos < dsq.len - 2) {
        // Count remaining (non-last) outgoing edges from x.
        var total: u32 = 0;
        for (0..K) |j| total += E[x][j];

        if (total == 0) {
            // Only the reserved last edge remains; use it.
            pos += 1;
            dsq[pos] = iE[x];
            x = iE[x];
        } else {
            // Pick a random edge from the remaining ones.
            var pick: u32 = rng.uniformInt(total);
            for (0..K) |j| {
                if (E[x][j] == 0) continue;
                if (pick < E[x][j]) {
                    pos += 1;
                    dsq[pos] = @intCast(j);
                    E[x][j] -= 1;
                    x = j;
                    break;
                }
                pick -= E[x][j];
            }
        }
    }

    // The last residue must match the original last residue.
    // With a correct Eulerian path it is guaranteed by the reserved last edges,
    // but we set it explicitly for clarity.  The final "last edge" from the
    // penultimate node leads to the original last residue.
    dsq[dsq.len - 1] = iE[x];
}

// --- Tests ---

test "randomIid: length and all codes < K for DNA" {
    const allocator = std.testing.allocator;
    var rng = Random.init(1);
    const freq = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    var seq = try randomIid(allocator, &rng, &alphabet_mod.dna, "test", 100, &freq);
    defer seq.deinit();

    try std.testing.expectEqual(@as(usize, 100), seq.len());
    for (seq.dsq) |code| {
        try std.testing.expect(code < alphabet_mod.dna.k);
    }
}

test "randomIid: zero-length sequence" {
    const allocator = std.testing.allocator;
    var rng = Random.init(2);
    const freq = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    var seq = try randomIid(allocator, &rng, &alphabet_mod.dna, "empty", 0, &freq);
    defer seq.deinit();

    try std.testing.expectEqual(@as(usize, 0), seq.len());
}

test "randomUniform: amino acid — correct length and canonical codes" {
    const allocator = std.testing.allocator;
    var rng = Random.init(3);
    var seq = try randomUniform(allocator, &rng, &alphabet_mod.amino, "prot", 200);
    defer seq.deinit();

    try std.testing.expectEqual(@as(usize, 200), seq.len());
    for (seq.dsq) |code| {
        try std.testing.expect(code < alphabet_mod.amino.k);
    }
}

test "randomUniform: DNA has correct length and canonical codes" {
    const allocator = std.testing.allocator;
    var rng = Random.init(4);
    var seq = try randomUniform(allocator, &rng, &alphabet_mod.dna, "dna_seq", 50);
    defer seq.deinit();

    try std.testing.expectEqual(@as(usize, 50), seq.len());
    for (seq.dsq) |code| {
        try std.testing.expect(code < alphabet_mod.dna.k);
    }
}

test "shuffle: preserves mono-residue composition" {
    var rng = Random.init(5);
    var dsq = [_]u8{ 0, 0, 1, 1, 2, 2, 3, 3 };

    var counts_before = [_]usize{0} ** 4;
    for (dsq) |c| counts_before[c] += 1;

    shuffle(&rng, &dsq);

    var counts_after = [_]usize{0} ** 4;
    for (dsq) |c| counts_after[c] += 1;

    try std.testing.expectEqualSlices(usize, &counts_before, &counts_after);
}

test "shuffle: length preserved" {
    var rng = Random.init(6);
    var dsq = [_]u8{ 0, 1, 2, 3, 0, 1 };
    const original_len = dsq.len;
    shuffle(&rng, &dsq);
    try std.testing.expectEqual(original_len, dsq.len);
}

test "shuffle: empty slice does not panic" {
    var rng = Random.init(7);
    var dsq = [_]u8{};
    shuffle(&rng, &dsq);
}

test "shuffle: single element unchanged" {
    var rng = Random.init(8);
    var dsq = [_]u8{3};
    shuffle(&rng, &dsq);
    try std.testing.expectEqual(@as(u8, 3), dsq[0]);
}

test "reproducibility: same seed produces same random sequence" {
    const allocator = std.testing.allocator;

    var rng1 = Random.init(42);
    var seq1 = try randomUniform(allocator, &rng1, &alphabet_mod.dna, "s1", 50);
    defer seq1.deinit();

    var rng2 = Random.init(42);
    var seq2 = try randomUniform(allocator, &rng2, &alphabet_mod.dna, "s1", 50);
    defer seq2.deinit();

    try std.testing.expectEqualSlices(u8, seq1.dsq, seq2.dsq);
}

test "reproducibility: different seeds produce different sequences" {
    const allocator = std.testing.allocator;

    var rng1 = Random.init(100);
    var seq1 = try randomUniform(allocator, &rng1, &alphabet_mod.dna, "s1", 50);
    defer seq1.deinit();

    var rng2 = Random.init(200);
    var seq2 = try randomUniform(allocator, &rng2, &alphabet_mod.dna, "s2", 50);
    defer seq2.deinit();

    // Extremely unlikely to match by chance with length 50.
    try std.testing.expect(!std.mem.eql(u8, seq1.dsq, seq2.dsq));
}

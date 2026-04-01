// WUSS (Washington University Secondary Structure) notation parser.
// Supports the four standard bracket types: <>, (), [], {}
// and pseudoknot alphabetic pairs Aa..Zz.
// Unpaired positions: . - _ ~

const std = @import("std");
const Allocator = std.mem.Allocator;

fn isOpenBracket(ch: u8) bool {
    return switch (ch) {
        '<', '(', '[', '{' => true,
        'a'...'z' => true,
        else => false,
    };
}

fn isCloseBracket(ch: u8) bool {
    return switch (ch) {
        '>', ')', ']', '}' => true,
        'A'...'Z' => true,
        else => false,
    };
}

/// Map a bracket character to an index identifying its type.
/// Opening and closing brackets of the same type share the same index.
fn bracketKind(ch: u8) u8 {
    return switch (ch) {
        '<', '>' => 0,
        '(', ')' => 1,
        '[', ']' => 2,
        '{', '}' => 3,
        'a'...'z' => 4 + (ch - 'a'),
        'A'...'Z' => 4 + (ch - 'A'),
        else => 0xff,
    };
}

/// Parse WUSS notation into a base-pair table.
/// pairs[i] = j means position i is paired with j; pairs[i] = -1 means unpaired.
/// Returns error.UnbalancedBrackets if brackets do not match.
pub fn parseToPairs(allocator: Allocator, wuss: []const u8) ![]i32 {
    const pairs = try allocator.alloc(i32, wuss.len);
    errdefer allocator.free(pairs);
    @memset(pairs, -1);

    // One stack per bracket kind (0..3 standard + 26 alphabetic pseudoknot kinds).
    const MAX_KINDS = 30;
    var stacks: [MAX_KINDS]std.ArrayList(usize) = @splat(.empty);
    defer for (&stacks) |*s| s.deinit(allocator);

    for (wuss, 0..) |ch, i| {
        // Reject non-printable characters, matching Easel's validation.
        if (ch < 0x20 or ch > 0x7E) return error.InvalidCharacter;

        if (isOpenBracket(ch)) {
            const kind = bracketKind(ch);
            if (kind >= MAX_KINDS) continue;
            try stacks[kind].append(allocator, i);
        } else if (isCloseBracket(ch)) {
            const kind = bracketKind(ch);
            if (kind >= MAX_KINDS) continue;
            const open_pos = stacks[kind].pop() orelse return error.UnbalancedBrackets;
            pairs[open_pos] = @intCast(i);
            pairs[i] = @intCast(open_pos);
        }
        // Unpaired printable characters (. - _ ~ etc.) leave pairs[i] = -1.
    }

    // All stacks must be empty for balanced notation.
    for (stacks) |s| {
        if (s.items.len != 0) return error.UnbalancedBrackets;
    }

    return pairs;
}

/// Convert a base-pair table back to simple dot-bracket notation using '(' and ')'.
pub fn pairsToSimple(allocator: Allocator, pairs: []const i32) ![]u8 {
    const buf = try allocator.alloc(u8, pairs.len);
    for (pairs, 0..) |partner, i| {
        if (partner < 0) {
            buf[i] = '.';
        } else if (@as(usize, @intCast(partner)) > i) {
            buf[i] = '(';
        } else {
            buf[i] = ')';
        }
    }
    return buf;
}

/// Count the number of base pairs (each pair counted once, via open brackets).
/// Returns error.InvalidCharacter for non-printable characters.
pub fn countPairs(wuss: []const u8) !usize {
    var count: usize = 0;
    for (wuss) |ch| {
        if (ch < 0x20 or ch > 0x7E) return error.InvalidCharacter;
        if (isOpenBracket(ch)) count += 1;
    }
    return count;
}

/// Return true if all brackets in wuss are balanced.
/// Returns error.InvalidCharacter for non-printable characters.
pub fn validate(wuss: []const u8) !bool {
    const MAX_KINDS = 30;
    var depth: [MAX_KINDS]usize = @splat(0);
    for (wuss) |ch| {
        if (ch < 0x20 or ch > 0x7E) return error.InvalidCharacter;
        if (isOpenBracket(ch)) {
            const kind = bracketKind(ch);
            if (kind < MAX_KINDS) depth[kind] += 1;
        } else if (isCloseBracket(ch)) {
            const kind = bracketKind(ch);
            if (kind >= MAX_KINDS) continue;
            if (depth[kind] == 0) return false;
            depth[kind] -= 1;
        }
    }
    for (depth) |d| if (d != 0) return false;
    return true;
}

/// Remove pseudoknot annotations from WUSS notation.
/// Pseudoknot pairs (alphabetic Aa..Zz) are replaced with '.'.
/// Standard bracket pairs (<>, (), [], {}) are preserved.
pub fn noPseudo(allocator: Allocator, wuss: []const u8) ![]u8 {
    const buf = try allocator.alloc(u8, wuss.len);
    for (wuss, 0..) |ch, i| {
        if ((ch >= 'a' and ch <= 'z') or (ch >= 'A' and ch <= 'Z')) {
            buf[i] = '.';
        } else {
            buf[i] = ch;
        }
    }
    return buf;
}

/// Reverse a WUSS string: reverse order and swap open/close brackets.
/// Useful for reverse complement of RNA secondary structure.
pub fn reverse(allocator: Allocator, wuss: []const u8) ![]u8 {
    const buf = try allocator.alloc(u8, wuss.len);
    for (wuss, 0..) |ch, i| {
        const swapped: u8 = switch (ch) {
            '<' => '>',
            '>' => '<',
            '(' => ')',
            ')' => '(',
            '[' => ']',
            ']' => '[',
            '{' => '}',
            '}' => '{',
            'a'...'z' => ch - 32, // lowercase -> uppercase
            'A'...'Z' => ch + 32, // uppercase -> lowercase
            else => ch,
        };
        buf[wuss.len - 1 - i] = swapped;
    }
    return buf;
}

/// Convert a base-pair table to rich WUSS notation with pseudoknot detection.
/// The main (non-crossing) structure uses '(' and ')'. Pairs that cross the
/// main structure get '<' and '>', further crossings get '[' ']', then '{' '}'.
/// Up to 4 bracket levels are supported; additional levels return error.TooManyPseudoknots.
pub fn pairsToWuss(allocator: Allocator, pairs: []const i32) ![]u8 {
    const len = pairs.len;
    const buf = try allocator.alloc(u8, len);
    errdefer allocator.free(buf);
    @memset(buf, '.');

    if (len == 0) return buf;

    const open_chars = [_]u8{ '(', '<', '[', '{' };
    const close_chars = [_]u8{ ')', '>', ']', '}' };
    const max_levels = open_chars.len;

    // Collect all base pairs (i < j) sorted by opening position.
    var pair_count: usize = 0;
    for (pairs) |p| {
        if (p >= 0) pair_count += 1;
    }
    pair_count /= 2; // Each pair counted twice.

    const pair_list = try allocator.alloc([2]usize, pair_count);
    defer allocator.free(pair_list);
    var idx: usize = 0;
    for (pairs, 0..) |p, i| {
        if (p >= 0 and @as(usize, @intCast(p)) > i) {
            pair_list[idx] = .{ i, @intCast(p) };
            idx += 1;
        }
    }

    // Assign bracket levels. level[k] is the bracket type for pair k.
    const levels = try allocator.alloc(u8, pair_count);
    defer allocator.free(levels);
    @memset(levels, 0xff); // Unassigned sentinel.

    // Greedy assignment: iterate through bracket levels and assign each pair
    // to the lowest level where it does not cross any previously-assigned pair
    // at the same level.
    for (0..max_levels) |level| {
        for (0..pair_count) |k| {
            if (levels[k] != 0xff) continue; // Already assigned.

            // Check if pair k crosses any pair already assigned to this level.
            var crosses = false;
            for (0..pair_count) |m| {
                if (levels[m] != level) continue;
                // Two pairs (a,b) and (c,d) cross iff a < c < b < d or c < a < d < b.
                const a = pair_list[k][0];
                const b = pair_list[k][1];
                const c = pair_list[m][0];
                const d = pair_list[m][1];
                if ((a < c and c < b and b < d) or (c < a and a < d and d < b)) {
                    crosses = true;
                    break;
                }
            }
            if (!crosses) {
                levels[k] = @intCast(level);
            }
        }
    }

    // Check that all pairs were assigned a level.
    for (levels) |l| {
        if (l == 0xff) return error.TooManyPseudoknots;
    }

    // Write the WUSS string.
    for (0..pair_count) |k| {
        const level = levels[k];
        buf[pair_list[k][0]] = open_chars[level];
        buf[pair_list[k][1]] = close_chars[level];
    }

    return buf;
}

/// Assign pairs to non-crossing layers using greedy planarity.
/// Layer 0 is the main (non-crossing) structure; higher layers are pseudoknots.
fn assignPairLayers(allocator: Allocator, all_pairs: [][2]usize, max_layers: u8) ![]u8 {
    const pair_count = all_pairs.len;
    const pair_layer = try allocator.alloc(u8, pair_count);
    errdefer allocator.free(pair_layer);
    @memset(pair_layer, 0xff);

    for (0..max_layers) |layer| {
        for (0..pair_count) |k| {
            if (pair_layer[k] != 0xff) continue;
            var crosses = false;
            for (0..pair_count) |m| {
                if (pair_layer[m] != layer) continue;
                const a = all_pairs[k][0];
                const b_val = all_pairs[k][1];
                const c = all_pairs[m][0];
                const d = all_pairs[m][1];
                if ((a < c and c < b_val and b_val < d) or (c < a and a < d and d < b_val)) {
                    crosses = true;
                    break;
                }
            }
            if (!crosses) {
                pair_layer[k] = @intCast(layer);
            }
        }
    }

    for (pair_layer) |l| {
        if (l == 0xff) {
            allocator.free(pair_layer);
            return error.TooManyPseudoknots;
        }
    }

    return pair_layer;
}

/// Collect validated base pairs (i < j) from a connectivity table.
fn collectPairs(allocator: Allocator, ct: []const ?usize, len: usize) !std.ArrayList([2]usize) {
    var pairs: std.ArrayList([2]usize) = .empty;
    errdefer pairs.deinit(allocator);

    for (0..len) |i| {
        if (ct[i]) |j| {
            if (j < len and j > i) {
                if (ct[j]) |back| {
                    if (back == i) {
                        try pairs.append(allocator, .{ i, j });
                    }
                }
            }
        }
    }
    return pairs;
}

/// Convert a connectivity table (CT) to full WUSS notation with structural context.
///
/// CT format: `ct[i] = j` means position i is paired with position j (0-indexed),
/// `ct[i] = null` means unpaired.
///
/// Paired positions use bracket types based on structural context:
/// - `<>` for stems closing hairpin loops (innermost pairs with no nested pairs)
/// - `()` for stems enclosing other stems (internal loops, multifurcations)
/// - `[]` and `{}` for pseudoknot pairs (crossing pairs at increasing depth)
///
/// Unpaired positions use context characters:
/// - `_` hairpin loop (between closing pair halves with no nested pairs inside)
/// - `-` bulge or interior loop (between stems, one or both sides)
/// - `,` multifurcation loop (between two or more stems inside enclosing pair)
/// - `:` external (outside all base pairs)
pub fn ct2wuss(allocator: Allocator, ct: []const ?usize, n: usize) ![]u8 {
    if (n == 0) return try allocator.alloc(u8, 0);

    const len = @min(n, ct.len);
    const buf = try allocator.alloc(u8, len);
    errdefer allocator.free(buf);
    @memset(buf, '.');

    var all_pairs_list = try collectPairs(allocator, ct, len);
    defer all_pairs_list.deinit(allocator);

    const pair_count = all_pairs_list.items.len;
    if (pair_count == 0) {
        @memset(buf, ':');
        return buf;
    }

    const all_pairs = all_pairs_list.items;

    const pair_layer = try assignPairLayers(allocator, all_pairs, 4);
    defer allocator.free(pair_layer);

    // Build a lookup: for each position, which pair index (or null).
    const pos_to_pair = try allocator.alloc(?usize, len);
    defer allocator.free(pos_to_pair);
    @memset(pos_to_pair, null);

    for (all_pairs, 0..) |pair, k| {
        pos_to_pair[pair[0]] = k;
        pos_to_pair[pair[1]] = k;
    }

    // For layer-0 pairs, classify structural context.
    for (0..pair_count) |k| {
        if (pair_layer[k] != 0) continue;
        const pi = all_pairs[k][0];
        const pj = all_pairs[k][1];

        var nested_stems: usize = 0;
        for (0..pair_count) |m| {
            if (m == k) continue;
            if (pair_layer[m] != 0) continue;
            if (all_pairs[m][0] > pi and all_pairs[m][1] < pj) {
                nested_stems += 1;
            }
        }

        if (nested_stems == 0) {
            buf[pi] = '<';
            buf[pj] = '>';
        } else {
            buf[pi] = '(';
            buf[pj] = ')';
        }
    }

    // Pseudoknot layers get [] and {}
    const pk_open = [_]u8{ '[', '{' };
    const pk_close = [_]u8{ ']', '}' };
    for (0..pair_count) |k| {
        const layer = pair_layer[k];
        if (layer == 0) continue;
        const bracket_idx = layer - 1;
        if (bracket_idx < pk_open.len) {
            buf[all_pairs[k][0]] = pk_open[bracket_idx];
            buf[all_pairs[k][1]] = pk_close[bracket_idx];
        } else {
            buf[all_pairs[k][0]] = '[';
            buf[all_pairs[k][1]] = ']';
        }
    }

    // Classify unpaired positions.
    for (0..len) |i| {
        if (pos_to_pair[i] != null) continue;

        // Find enclosing layer-0 base pair (innermost).
        var best_k: ?usize = null;
        var best_span: usize = std.math.maxInt(usize);
        for (all_pairs, 0..) |pair, k| {
            if (pair_layer[k] != 0) continue;
            if (pair[0] < i and i < pair[1]) {
                const span = pair[1] - pair[0];
                if (span < best_span) {
                    best_span = span;
                    best_k = k;
                }
            }
        }

        if (best_k == null) {
            buf[i] = ':';
        } else {
            const ek = best_k.?;
            const ei = all_pairs[ek][0];
            const ej = all_pairs[ek][1];

            // Count directly nested layer-0 stems inside (ei, ej).
            var direct_stems: usize = 0;
            for (all_pairs, 0..) |pair, m| {
                if (pair_layer[m] != 0) continue;
                if (m == ek) continue;
                const mi = pair[0];
                const mj = pair[1];
                if (mi <= ei or mj >= ej) continue;

                var is_direct = true;
                for (all_pairs, 0..) |other, o| {
                    if (pair_layer[o] != 0) continue;
                    if (o == ek or o == m) continue;
                    if (other[0] > ei and other[1] < ej and other[0] < mi and mj < other[1]) {
                        is_direct = false;
                        break;
                    }
                }
                if (is_direct) direct_stems += 1;
            }

            if (direct_stems == 0) {
                buf[i] = '_';
            } else if (direct_stems == 1) {
                buf[i] = '-';
            } else {
                buf[i] = ',';
            }
        }
    }

    return buf;
}

/// Convert a connectivity table (CT) to simplified WUSS notation.
///
/// All base pairs use `<>`, pseudoknot pairs use `Aa`..`Zz` alphabetic notation,
/// and unpaired positions are `.`.
///
/// CT format: `ct[i] = j` means position i is paired with position j (0-indexed),
/// `ct[i] = null` means unpaired.
pub fn ct2simplewuss(allocator: Allocator, ct: []const ?usize, n: usize) ![]u8 {
    if (n == 0) return try allocator.alloc(u8, 0);

    const len = @min(n, ct.len);
    const buf = try allocator.alloc(u8, len);
    errdefer allocator.free(buf);
    @memset(buf, '.');

    var all_pairs_list = try collectPairs(allocator, ct, len);
    defer all_pairs_list.deinit(allocator);

    const pair_count = all_pairs_list.items.len;
    if (pair_count == 0) return buf;

    const all_pairs = all_pairs_list.items;

    const pair_layer = try assignPairLayers(allocator, all_pairs, 27);
    defer allocator.free(pair_layer);

    for (0..pair_count) |k| {
        const layer = pair_layer[k];
        const pi = all_pairs[k][0];
        const pj = all_pairs[k][1];
        if (layer == 0) {
            buf[pi] = '<';
            buf[pj] = '>';
        } else {
            const letter: u8 = @intCast(layer - 1);
            if (letter < 26) {
                buf[pi] = 'a' + letter;
                buf[pj] = 'A' + letter;
            } else {
                buf[pi] = 'z';
                buf[pj] = 'Z';
            }
        }
    }

    return buf;
}

// --- Tests ---

test "parseToPairs: ((...)) gives expected pair table" {
    // pairs[0]=6, pairs[1]=5, pairs[2..4]=-1, pairs[5]=1, pairs[6]=0
    const allocator = std.testing.allocator;
    const pairs = try parseToPairs(allocator, "((...))");
    defer allocator.free(pairs);
    try std.testing.expectEqual(@as(i32, 6), pairs[0]);
    try std.testing.expectEqual(@as(i32, 5), pairs[1]);
    try std.testing.expectEqual(@as(i32, -1), pairs[2]);
    try std.testing.expectEqual(@as(i32, -1), pairs[3]);
    try std.testing.expectEqual(@as(i32, -1), pairs[4]);
    try std.testing.expectEqual(@as(i32, 1), pairs[5]);
    try std.testing.expectEqual(@as(i32, 0), pairs[6]);
}

test "parseToPairs: ((....)) nested pairs" {
    const allocator = std.testing.allocator;
    const pairs = try parseToPairs(allocator, "((....))");
    defer allocator.free(pairs);
    try std.testing.expectEqual(@as(i32, 7), pairs[0]);
    try std.testing.expectEqual(@as(i32, 6), pairs[1]);
    for (2..6) |i| try std.testing.expectEqual(@as(i32, -1), pairs[i]);
    try std.testing.expectEqual(@as(i32, 1), pairs[6]);
    try std.testing.expectEqual(@as(i32, 0), pairs[7]);
}

test "parseToPairs: all unpaired" {
    const allocator = std.testing.allocator;
    const pairs = try parseToPairs(allocator, "....");
    defer allocator.free(pairs);
    for (pairs) |p| try std.testing.expectEqual(@as(i32, -1), p);
}

test "parseToPairs: mixed bracket types" {
    const allocator = std.testing.allocator;
    // "<()>" -> 0 pairs with 3, 1 pairs with 2
    const pairs = try parseToPairs(allocator, "<()>");
    defer allocator.free(pairs);
    try std.testing.expectEqual(@as(i32, 3), pairs[0]);
    try std.testing.expectEqual(@as(i32, 2), pairs[1]);
    try std.testing.expectEqual(@as(i32, 1), pairs[2]);
    try std.testing.expectEqual(@as(i32, 0), pairs[3]);
}

test "parseToPairs: unbalanced returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(error.UnbalancedBrackets, parseToPairs(allocator, "((.)"));
    try std.testing.expectError(error.UnbalancedBrackets, parseToPairs(allocator, ".)"));
}

test "parseToPairs: non-printable characters rejected" {
    const allocator = std.testing.allocator;
    // Control character (0x01) should be rejected
    try std.testing.expectError(error.InvalidCharacter, parseToPairs(allocator, "((\x01..))"));
    // DEL (0x7F) should be rejected
    try std.testing.expectError(error.InvalidCharacter, parseToPairs(allocator, "((\x7F..))"));
    // High byte (0x80) should be rejected
    try std.testing.expectError(error.InvalidCharacter, parseToPairs(allocator, "((\x80..))"));
}

test "countPairs: two nested pairs" {
    try std.testing.expectEqual(@as(usize, 2), try countPairs("((....))"));
}

test "countPairs: zero pairs" {
    try std.testing.expectEqual(@as(usize, 0), try countPairs("...."));
}

test "countPairs: mixed brackets" {
    try std.testing.expectEqual(@as(usize, 2), try countPairs("<()>"));
}

test "countPairs: non-printable character rejected" {
    try std.testing.expectError(error.InvalidCharacter, countPairs("((\x01..))"));
}

test "validate: balanced strings" {
    try std.testing.expect(try validate("((....))"));
    try std.testing.expect(try validate("...."));
    try std.testing.expect(try validate("<()>"));
    try std.testing.expect(try validate(""));
}

test "validate: unbalanced strings" {
    try std.testing.expect(!try validate("((.)"));
    try std.testing.expect(!try validate(".)"));
}

test "validate: non-printable character rejected" {
    try std.testing.expectError(error.InvalidCharacter, validate("((\x7F..))"));
}

test "pairsToSimple: round-trip from ((....))'" {
    const allocator = std.testing.allocator;
    const pairs = try parseToPairs(allocator, "((....))");
    defer allocator.free(pairs);
    const simple = try pairsToSimple(allocator, pairs);
    defer allocator.free(simple);
    try std.testing.expectEqualStrings("((....))", simple);
}

test "pairsToSimple: all unpaired" {
    const allocator = std.testing.allocator;
    const raw = [_]i32{ -1, -1, -1 };
    const simple = try pairsToSimple(allocator, &raw);
    defer allocator.free(simple);
    try std.testing.expectEqualStrings("...", simple);
}

test "pairsToWuss: simple nested structure uses parens" {
    const allocator = std.testing.allocator;
    const pairs = try parseToPairs(allocator, "((....))");
    defer allocator.free(pairs);
    const wuss = try pairsToWuss(allocator, pairs);
    defer allocator.free(wuss);
    try std.testing.expectEqualStrings("((....))", wuss);
}

test "pairsToWuss: all unpaired" {
    const allocator = std.testing.allocator;
    const wuss = try pairsToWuss(allocator, &[_]i32{ -1, -1, -1, -1 });
    defer allocator.free(wuss);
    try std.testing.expectEqualStrings("....", wuss);
}

test "pairsToWuss: empty input" {
    const allocator = std.testing.allocator;
    const wuss = try pairsToWuss(allocator, &[_]i32{});
    defer allocator.free(wuss);
    try std.testing.expectEqualStrings("", wuss);
}

test "pairsToWuss: pseudoknot gets angle brackets" {
    // (0,5) and (2,7) cross because 0 < 2 < 5 < 7.
    //   pos:   0  1  2  3  4  5  6  7
    //   pairs: 5 -1  7 -1 -1  0 -1  2
    const allocator = std.testing.allocator;
    const pairs = [_]i32{ 5, -1, 7, -1, -1, 0, -1, 2 };
    const wuss = try pairsToWuss(allocator, &pairs);
    defer allocator.free(wuss);
    // (0,5) -> '(' ')' at level 0, (2,7) crosses it -> '<' '>' at level 1.
    try std.testing.expectEqualStrings("(.<..).>", wuss);
}

test "pairsToWuss: round-trip non-crossing pairs all get parens" {
    // Parse "((..<<..>>..))". The pairs (0,13) (1,12) (4,9) (5,8) are all nested,
    // none cross, so pairsToWuss should assign all to level 0 (parens).
    const allocator = std.testing.allocator;
    const input = "((..<<..>>..))";
    const pairs = try parseToPairs(allocator, input);
    defer allocator.free(pairs);
    const wuss = try pairsToWuss(allocator, pairs);
    defer allocator.free(wuss);
    try std.testing.expectEqualStrings("((..((..))..))", wuss);
}

test "ct2wuss: simple hairpin" {
    // CT: 0-8, 1-7, 2-6 paired, 3,4,5 unpaired (hairpin loop)
    const allocator = std.testing.allocator;
    const ct = [_]?usize{ 8, 7, 6, null, null, null, 2, 1, 0 };
    const result = try ct2wuss(allocator, &ct, 9);
    defer allocator.free(result);
    // Inner pair (2,6) encloses no other pairs -> <>, outer pairs enclose stems -> ()
    // Positions 3,4,5 are in hairpin -> '_'
    try std.testing.expectEqualStrings("((<___>))", result);
}

test "ct2wuss: external unpaired" {
    // Leading/trailing unpaired should be ':'
    const allocator = std.testing.allocator;
    const ct = [_]?usize{ null, null, null, 11, 10, 9, null, null, null, 5, 4, 3, null, null, null };
    const result = try ct2wuss(allocator, &ct, 15);
    defer allocator.free(result);
    try std.testing.expectEqualStrings(":::((<___>)):::", result);
}

test "ct2wuss: bulge/interior loop" {
    // Pair (0,11) (1,10) enclose pair (3,8)(4,7), positions 2 and 9 are bulge/interior
    const allocator = std.testing.allocator;
    const ct = [_]?usize{ 11, 10, null, 8, 7, null, null, 4, 3, null, 1, 0 };
    const result = try ct2wuss(allocator, &ct, 12);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("((-(<__>)-))", result);
}

test "ct2wuss: multifurcation loop" {
    // Two stems inside one enclosing pair
    const allocator = std.testing.allocator;
    const ct = [_]?usize{
        19,   null, 8,    7,    null, null, null, 3,    2,    null,
        null, 17,   16,   null, null, null, 12,   11,   null, 0,
    };
    const result = try ct2wuss(allocator, &ct, 20);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("(,(<___>),,(<___>),)", result);
}

test "ct2wuss: all unpaired" {
    const allocator = std.testing.allocator;
    const ct = [_]?usize{ null, null, null, null };
    const result = try ct2wuss(allocator, &ct, 4);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("::::", result);
}

test "ct2wuss: empty input" {
    const allocator = std.testing.allocator;
    const ct = [_]?usize{};
    const result = try ct2wuss(allocator, &ct, 0);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("", result);
}

test "ct2wuss: single hairpin no nesting" {
    const allocator = std.testing.allocator;
    const ct = [_]?usize{ 5, null, null, null, null, 0 };
    const result = try ct2wuss(allocator, &ct, 6);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("<____>", result);
}

test "ct2simplewuss: basic pairs" {
    const allocator = std.testing.allocator;
    const ct = [_]?usize{ 8, 7, 6, null, null, null, 2, 1, 0 };
    const result = try ct2simplewuss(allocator, &ct, 9);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("<<<...>>>", result);
}

test "ct2simplewuss: all unpaired" {
    const allocator = std.testing.allocator;
    const ct = [_]?usize{ null, null, null };
    const result = try ct2simplewuss(allocator, &ct, 3);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("...", result);
}

test "ct2simplewuss: pseudoknot uses alphabetic pairs" {
    // Crossing pairs: (0,5) and (2,7) cross.
    const allocator = std.testing.allocator;
    const ct = [_]?usize{ 5, null, 7, null, null, 0, null, 2 };
    const result = try ct2simplewuss(allocator, &ct, 8);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("<.a..>.A", result);
}

test "ct2simplewuss: empty input" {
    const allocator = std.testing.allocator;
    const ct = [_]?usize{};
    const result = try ct2simplewuss(allocator, &ct, 0);
    defer allocator.free(result);
    try std.testing.expectEqualStrings("", result);
}

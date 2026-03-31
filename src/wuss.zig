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
        // Unpaired characters (. - _ ~) leave pairs[i] = -1.
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
pub fn countPairs(wuss: []const u8) usize {
    var count: usize = 0;
    for (wuss) |ch| {
        if (isOpenBracket(ch)) count += 1;
    }
    return count;
}

/// Return true if all brackets in wuss are balanced.
pub fn validate(wuss: []const u8) bool {
    const MAX_KINDS = 30;
    var depth: [MAX_KINDS]usize = @splat(0);
    for (wuss) |ch| {
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

// --- Tests ---

test "parseToPairs: ((...)) gives expected pair table" {
    // pairs[0]=6, pairs[1]=5, pairs[2..4]=-1, pairs[5]=1, pairs[6]=0
    const allocator = std.testing.allocator;
    const pairs = try parseToPairs(allocator, "((...))" );
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

test "countPairs: two nested pairs" {
    try std.testing.expectEqual(@as(usize, 2), countPairs("((....))"));
}

test "countPairs: zero pairs" {
    try std.testing.expectEqual(@as(usize, 0), countPairs("...."));
}

test "countPairs: mixed brackets" {
    try std.testing.expectEqual(@as(usize, 2), countPairs("<()>"));
}

test "validate: balanced strings" {
    try std.testing.expect(validate("((....))"));
    try std.testing.expect(validate("...."));
    try std.testing.expect(validate("<()>"));
    try std.testing.expect(validate(""));
}

test "validate: unbalanced strings" {
    try std.testing.expect(!validate("((.)"));
    try std.testing.expect(!validate(".)"));
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

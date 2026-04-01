// Bitwise pack/unpack for the Easel dsqdata packet format.
//
// Packet format (uint32):
//   Bit 31:    EOD flag (end-of-data)
//   Bit 30:    type (0 = 2-bit encoding, 1 = 5-bit encoding)
//   Bits 29-0: packed residue data
//
// 5-bit mode: 6 residues per full packet at bit positions 25, 20, 15, 10, 5, 0.
//             EOD packet contains 0-6 residues padded with sentinel 31.
//
// 2-bit mode: 15 residues per full packet at bit positions 28, 26, ..., 2, 0.
//             Switches to 5-bit mid-sequence when degenerate residues (code >= 4)
//             are encountered. 2-bit EOD packets are always 15 residues (full).

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const EOD: u32 = 1 << 31;
pub const FIVEBIT: u32 = 1 << 30;

// Sentinel value used for padding in 5-bit EOD packets.
const SENTINEL_5BIT: u8 = 31;

// Number of residues packed per full packet in each mode.
const RESIDUES_PER_5BIT_PACKET: usize = 6;
const RESIDUES_PER_2BIT_PACKET: usize = 15;

// --- 5-bit pack/unpack ---

/// Pack a digital sequence into 5-bit packets.
/// Empty sequence produces a single all-sentinel EOD packet.
/// Caller owns the returned slice.
pub fn pack5(allocator: Allocator, dsq: []const u8) ![]u32 {
    // Number of packets: ceil(len / 6), minimum 1 (for empty or any input).
    const n_packets = if (dsq.len == 0) 1 else (dsq.len + RESIDUES_PER_5BIT_PACKET - 1) / RESIDUES_PER_5BIT_PACKET;

    const packets = try allocator.alloc(u32, n_packets);
    errdefer allocator.free(packets);

    var i: usize = 0; // index into dsq
    var p: usize = 0; // packet index

    while (p < n_packets) : (p += 1) {
        const is_last = (p == n_packets - 1);

        // Collect up to 6 residues for this packet.
        var residues: [RESIDUES_PER_5BIT_PACKET]u8 = [_]u8{SENTINEL_5BIT} ** RESIDUES_PER_5BIT_PACKET;
        var count: usize = 0;
        while (count < RESIDUES_PER_5BIT_PACKET and i < dsq.len) {
            residues[count] = dsq[i];
            count += 1;
            i += 1;
        }

        // Pack residues into the 30-bit data field.
        // Residue 0 at bits 25-21, residue 1 at 20-16, ..., residue 5 at 4-0.
        var packet: u32 = 0;
        for (residues, 0..) |r, slot| {
            const shift: u5 = @intCast((RESIDUES_PER_5BIT_PACKET - 1 - slot) * 5);
            packet |= @as(u32, r) << shift;
        }

        packet |= FIVEBIT;
        if (is_last) packet |= EOD;

        packets[p] = packet;
    }

    return packets;
}

/// Unpack 5-bit packets into a digital sequence.
/// seq_len must match the actual number of residues encoded.
/// Caller owns the returned slice.
pub fn unpack5(allocator: Allocator, packets: []const u32, seq_len: u64) ![]u8 {
    const dsq = try allocator.alloc(u8, @intCast(seq_len));
    errdefer allocator.free(dsq);

    var out: usize = 0;

    for (packets) |packet| {
        const is_eod = (packet & EOD) != 0;

        for (0..RESIDUES_PER_5BIT_PACKET) |slot| {
            const shift: u5 = @intCast((RESIDUES_PER_5BIT_PACKET - 1 - slot) * 5);
            const r: u8 = @intCast((packet >> shift) & 0x1F);

            if (is_eod and r == SENTINEL_5BIT) break; // padding sentinel

            if (out >= seq_len) break;
            dsq[out] = r;
            out += 1;
        }

        if (is_eod) break;
    }

    return dsq;
}

// --- 2-bit pack/unpack ---

/// Pack a digital sequence into 2-bit packets (with 5-bit fallback for
/// degenerate residues). Empty sequence produces a single all-sentinel
/// EOD 5-bit packet (same as pack5 empty).
/// Caller owns the returned slice.
pub fn pack2(allocator: Allocator, dsq: []const u8) ![]u32 {
    if (dsq.len == 0) {
        return pack5(allocator, dsq);
    }

    var packets: std.ArrayList(u32) = .empty;
    errdefer packets.deinit(allocator);

    var i: usize = 0;

    while (i < dsq.len) {
        const remaining = dsq.len - i;

        // Decide 2-bit vs 5-bit for this chunk:
        // Use 2-bit only when >=15 residues remain AND none of the next 15 are degenerate.
        const use_2bit = blk: {
            if (remaining < RESIDUES_PER_2BIT_PACKET) break :blk false;
            for (0..RESIDUES_PER_2BIT_PACKET) |k| {
                if (dsq[i + k] > 3) break :blk false;
            }
            break :blk true;
        };

        if (use_2bit) {
            // Pack 15 residues into bits 29-0 (pairs at positions 28-29, 26-27, ..., 0-1).
            var packet: u32 = 0;
            for (0..RESIDUES_PER_2BIT_PACKET) |k| {
                const shift: u5 = @intCast((RESIDUES_PER_2BIT_PACKET - 1 - k) * 2);
                packet |= @as(u32, dsq[i + k]) << shift;
            }

            const is_last = (i + RESIDUES_PER_2BIT_PACKET >= dsq.len);
            if (is_last) packet |= EOD;
            // FIVEBIT bit is 0 for 2-bit packets.

            try packets.append(allocator, packet);
            i += RESIDUES_PER_2BIT_PACKET;
        } else {
            // 5-bit packet: pack up to 6 residues.
            var residues: [RESIDUES_PER_5BIT_PACKET]u8 = [_]u8{SENTINEL_5BIT} ** RESIDUES_PER_5BIT_PACKET;
            var count: usize = 0;
            while (count < RESIDUES_PER_5BIT_PACKET and i + count < dsq.len) {
                residues[count] = dsq[i + count];
                count += 1;
            }

            var packet: u32 = 0;
            for (residues, 0..) |r, slot| {
                const shift: u5 = @intCast((RESIDUES_PER_5BIT_PACKET - 1 - slot) * 5);
                packet |= @as(u32, r) << shift;
            }

            packet |= FIVEBIT;

            i += count;
            const is_last = (i >= dsq.len);
            if (is_last) packet |= EOD;

            try packets.append(allocator, packet);
        }
    }

    return packets.toOwnedSlice(allocator);
}

/// Unpack mixed 2-bit/5-bit packets into a digital sequence.
/// seq_len must match the actual number of residues encoded.
/// Caller owns the returned slice.
pub fn unpack2(allocator: Allocator, packets: []const u32, seq_len: u64) ![]u8 {
    const dsq = try allocator.alloc(u8, @intCast(seq_len));
    errdefer allocator.free(dsq);

    var out: usize = 0;

    for (packets) |packet| {
        const is_eod = (packet & EOD) != 0;
        const is_5bit = (packet & FIVEBIT) != 0;

        if (is_5bit) {
            // 5-bit packet: 6 residues, stop at sentinel 31 if EOD.
            for (0..RESIDUES_PER_5BIT_PACKET) |slot| {
                const shift: u5 = @intCast((RESIDUES_PER_5BIT_PACKET - 1 - slot) * 5);
                const r: u8 = @intCast((packet >> shift) & 0x1F);

                if (is_eod and r == SENTINEL_5BIT) break;
                if (out >= seq_len) break;
                dsq[out] = r;
                out += 1;
            }
        } else {
            // 2-bit packet: always 15 residues (including EOD packets).
            for (0..RESIDUES_PER_2BIT_PACKET) |k| {
                const shift: u5 = @intCast((RESIDUES_PER_2BIT_PACKET - 1 - k) * 2);
                const r: u8 = @intCast((packet >> shift) & 0x3);

                if (out >= seq_len) break;
                dsq[out] = r;
                out += 1;
            }
        }

        if (is_eod) break;
    }

    return dsq;
}

// --- Tests ---

test "unpack5: 3 amino acids (single EOD packet)" {
    const allocator = std.testing.allocator;

    // 3 residues [1, 2, 3] + 3 sentinels (31), with EOD | FIVEBIT.
    // Bit layout: r0=1 at 25, r1=2 at 20, r2=3 at 15, r3=31 at 10, r4=31 at 5, r5=31 at 0.
    const packet: u32 = EOD | FIVEBIT |
        (@as(u32, 1) << 25) |
        (@as(u32, 2) << 20) |
        (@as(u32, 3) << 15) |
        (@as(u32, 31) << 10) |
        (@as(u32, 31) << 5) |
        (@as(u32, 31) << 0);

    const result = try unpack5(allocator, &[_]u32{packet}, 3);
    defer allocator.free(result);

    try std.testing.expectEqualSlices(u8, &[_]u8{ 1, 2, 3 }, result);
}

test "unpack2: 15 canonical DNA residues (single 2-bit EOD packet)" {
    const allocator = std.testing.allocator;

    // 15 residues: 0,1,2,3 repeated, then 0,1,2.
    const dsq = [_]u8{ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2 };
    var packet: u32 = EOD; // 2-bit (no FIVEBIT)
    for (dsq, 0..) |r, k| {
        const shift: u5 = @intCast((RESIDUES_PER_2BIT_PACKET - 1 - k) * 2);
        packet |= @as(u32, r) << shift;
    }

    const result = try unpack2(allocator, &[_]u32{packet}, 15);
    defer allocator.free(result);

    try std.testing.expectEqualSlices(u8, &dsq, result);
}

test "unpack2: mixed 2-bit and 5-bit packets (degenerate in middle)" {
    const allocator = std.testing.allocator;

    // Build a sequence: 15 canonical residues, then degenerate code 5, then 5 more canonical.
    // This would produce: one 2-bit packet, then one 5-bit packet with [5,0,1,2,3,0], EOD.
    // We hand-craft the packets to verify the decoder.

    // Packet 0: 2-bit, 15 x residue 1 (no FIVEBIT, no EOD).
    var p0: u32 = 0;
    for (0..15) |k| {
        const shift: u5 = @intCast((14 - k) * 2);
        p0 |= @as(u32, 1) << shift;
    }

    // Packet 1: 5-bit EOD with residues [5, 0, 1, 2, 3] and one sentinel.
    const p1: u32 = EOD | FIVEBIT |
        (@as(u32, 5) << 25) |
        (@as(u32, 0) << 20) |
        (@as(u32, 1) << 15) |
        (@as(u32, 2) << 10) |
        (@as(u32, 3) << 5) |
        (@as(u32, SENTINEL_5BIT) << 0);

    const result = try unpack2(allocator, &[_]u32{ p0, p1 }, 20);
    defer allocator.free(result);

    const expected = [_]u8{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0, 1, 2, 3 };
    try std.testing.expectEqualSlices(u8, &expected, result);
}

test "pack5/unpack5 round-trip: amino acid sequence (11 residues)" {
    const allocator = std.testing.allocator;

    const dsq = [_]u8{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

    const packets = try pack5(allocator, &dsq);
    defer allocator.free(packets);

    // 11 residues -> ceil(11/6) = 2 packets.
    try std.testing.expectEqual(@as(usize, 2), packets.len);

    // First packet: not EOD, is FIVEBIT.
    try std.testing.expect((packets[0] & EOD) == 0);
    try std.testing.expect((packets[0] & FIVEBIT) != 0);

    // Last packet: EOD and FIVEBIT.
    try std.testing.expect((packets[1] & EOD) != 0);
    try std.testing.expect((packets[1] & FIVEBIT) != 0);

    const result = try unpack5(allocator, packets, dsq.len);
    defer allocator.free(result);

    try std.testing.expectEqualSlices(u8, &dsq, result);
}

test "pack2/unpack2 round-trip: canonical DNA (30 residues)" {
    const allocator = std.testing.allocator;

    // 30 canonical residues (codes 0-3).
    const dsq = [_]u8{ 0, 1, 2, 3 } ** 7 ++ [_]u8{ 0, 1 };

    const packets = try pack2(allocator, &dsq);
    defer allocator.free(packets);

    // 30 residues -> exactly 2 full 2-bit packets.
    try std.testing.expectEqual(@as(usize, 2), packets.len);
    try std.testing.expect((packets[0] & FIVEBIT) == 0);
    try std.testing.expect((packets[1] & FIVEBIT) == 0);
    try std.testing.expect((packets[1] & EOD) != 0);

    const result = try unpack2(allocator, packets, dsq.len);
    defer allocator.free(result);

    try std.testing.expectEqualSlices(u8, &dsq, result);
}

test "pack2/unpack2 round-trip: DNA with degenerate at position 14" {
    const allocator = std.testing.allocator;

    // 17 residues: positions 0-13 canonical (0-3), position 14 degenerate (code 5),
    // positions 15-16 canonical.
    var dsq: [17]u8 = undefined;
    for (0..14) |k| dsq[k] = @intCast(k % 4);
    dsq[14] = 5; // degenerate
    dsq[15] = 2;
    dsq[16] = 3;

    const packets = try pack2(allocator, &dsq);
    defer allocator.free(packets);

    // Positions 0-13: only 14 canonical before degenerate at 14, < 15 -> falls to 5-bit.
    // All packets should be 5-bit.
    for (packets) |pkt| {
        try std.testing.expect((pkt & FIVEBIT) != 0);
    }

    const result = try unpack2(allocator, packets, dsq.len);
    defer allocator.free(result);

    try std.testing.expectEqualSlices(u8, &dsq, result);
}

test "pack5/unpack5 round-trip: empty sequence" {
    const allocator = std.testing.allocator;

    const packets = try pack5(allocator, &[_]u8{});
    defer allocator.free(packets);

    // Empty -> 1 packet, all sentinels, EOD | FIVEBIT.
    try std.testing.expectEqual(@as(usize, 1), packets.len);
    try std.testing.expect((packets[0] & EOD) != 0);
    try std.testing.expect((packets[0] & FIVEBIT) != 0);

    const result = try unpack5(allocator, packets, 0);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 0), result.len);
}

test "pack5/unpack5 round-trip: exact 6-residue boundary" {
    const allocator = std.testing.allocator;

    const dsq = [_]u8{ 0, 1, 2, 3, 4, 5 };

    const packets = try pack5(allocator, &dsq);
    defer allocator.free(packets);

    // 6 residues -> exactly 1 full packet (EOD | FIVEBIT, no sentinels needed).
    try std.testing.expectEqual(@as(usize, 1), packets.len);
    try std.testing.expect((packets[0] & EOD) != 0);
    try std.testing.expect((packets[0] & FIVEBIT) != 0);

    const result = try unpack5(allocator, packets, dsq.len);
    defer allocator.free(result);

    try std.testing.expectEqualSlices(u8, &dsq, result);
}

test "pack2/unpack2 round-trip: exact 15-residue boundary" {
    const allocator = std.testing.allocator;

    // 15 canonical residues -> exactly 1 full 2-bit packet.
    const dsq = [_]u8{ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2 };

    const packets = try pack2(allocator, &dsq);
    defer allocator.free(packets);

    try std.testing.expectEqual(@as(usize, 1), packets.len);
    try std.testing.expect((packets[0] & FIVEBIT) == 0);
    try std.testing.expect((packets[0] & EOD) != 0);

    const result = try unpack2(allocator, packets, dsq.len);
    defer allocator.free(result);

    try std.testing.expectEqualSlices(u8, &dsq, result);
}

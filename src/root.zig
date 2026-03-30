pub const alphabet = @import("alphabet.zig");
pub const sequence = @import("sequence.zig");
pub const io = @import("io.zig");
pub const msa = @import("msa.zig");
pub const simd = @import("simd.zig");
pub const stats = @import("stats.zig");
pub const score_matrix = @import("score_matrix.zig");
pub const util = @import("util.zig");

// Include tests from all sub-modules in `zig build test`.
test {
    _ = alphabet;
    _ = sequence;
    _ = io;
    _ = msa;
    _ = simd;
    _ = stats;
    _ = score_matrix;
    _ = util;
}

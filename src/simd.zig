/// SIMD-accelerated vector operations.
/// Provides VectorOps(T) — a comptime generic namespace covering element-wise
/// operations, reductions, and (for float types) statistical operations such as
/// normalize, entropy, and log-sum-exp.
pub const vector_ops = @import("simd/vector_ops.zig");

pub const VectorOps = vector_ops.VectorOps;
pub const VecF64 = vector_ops.VecF64;
pub const VecF32 = vector_ops.VecF32;
pub const VecI32 = vector_ops.VecI32;
pub const VecU8 = vector_ops.VecU8;

/// SIMD primitive operations for HMMER-style striped DP algorithms:
/// lane shifts, horizontal reductions, comparisons, select/blend,
/// and saturating arithmetic.
pub const primitives = @import("simd/primitives.zig");

pub const rightShift = primitives.rightShift;
pub const leftShift = primitives.leftShift;
pub const hmax = primitives.hmax;
pub const hmin = primitives.hmin;
pub const hsum = primitives.hsum;
pub const anyGt = primitives.anyGt;
pub const anyGte = primitives.anyGte;
pub const select = primitives.select;
pub const addSat = primitives.addSat;
pub const subSat = primitives.subSat;

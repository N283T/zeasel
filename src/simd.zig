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

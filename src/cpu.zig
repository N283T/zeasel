// CPU feature detection for runtime SIMD dispatch.
//
// At comptime, Zig knows the target features. This module provides
// both comptime and runtime detection. Runtime detection uses CPUID
// on x86_64 (following Easel's esl_cpu.c approach) and assumes
// NEON on aarch64 (always available).

const std = @import("std");
const builtin = @import("builtin");

/// SIMD instruction set support flags.
pub const SimdSupport = struct {
    sse: bool = false,
    sse2: bool = false,
    sse4_1: bool = false,
    avx: bool = false,
    avx2: bool = false,
    avx512f: bool = false,
    neon: bool = false,

    /// Return the name of the best available SIMD implementation.
    pub fn bestImplementation(self: SimdSupport) []const u8 {
        if (self.avx512f) return "AVX512";
        if (self.avx2 or self.avx) return "AVX";
        if (self.sse4_1) return "SSE4";
        if (self.sse2 or self.sse) return "SSE";
        if (self.neon) return "NEON";
        return "none";
    }

    /// Return the best SIMD vector width in bytes.
    pub fn vectorWidth(self: SimdSupport) usize {
        if (self.avx512f) return 64;
        if (self.avx2 or self.avx) return 32;
        if (self.sse or self.sse2 or self.sse4_1 or self.neon) return 16;
        return 8; // scalar fallback
    }
};

// ---------------------------------------------------------------------------
// x86/x86_64 CPUID helpers
// ---------------------------------------------------------------------------

const CpuidResult = struct {
    eax: u32,
    ebx: u32,
    ecx: u32,
    edx: u32,
};

/// Execute CPUID instruction with given leaf and sub-leaf.
/// Only available on x86/x86_64.
inline fn cpuid(leaf: u32, sub_leaf: u32) CpuidResult {
    var eax: u32 = undefined;
    var ebx: u32 = undefined;
    var ecx: u32 = undefined;
    var edx: u32 = undefined;

    asm volatile ("cpuid"
        : [_] "={eax}" (eax),
          [_] "={ebx}" (ebx),
          [_] "={ecx}" (ecx),
          [_] "={edx}" (edx),
        : [_] "{eax}" (leaf),
          [_] "{ecx}" (sub_leaf),
    );

    return .{ .eax = eax, .ebx = ebx, .ecx = ecx, .edx = edx };
}

/// Read XCR0 via XGETBV (extended control register 0).
inline fn getXCR0() u32 {
    return asm volatile (
        \\ xor %%ecx, %%ecx
        \\ xgetbv
        : [_] "={eax}" (-> u32),
        :
        : .{ .edx = true, .ecx = true });
}

/// Check OS support for YMM registers (AVX) via XGETBV.
inline fn xgetbvChecksYmm() bool {
    const xcr0 = getXCR0();
    const ymm_xmm: u32 = (1 << 2) | (1 << 1);
    return (xcr0 & ymm_xmm) == ymm_xmm;
}

/// Check OS support for ZMM registers (AVX-512) via XGETBV.
inline fn xgetbvChecksZmm() bool {
    const xcr0 = getXCR0();
    const zmm_ymm_xmm: u32 = (7 << 5) | (1 << 2) | (1 << 1);
    return (xcr0 & zmm_ymm_xmm) == zmm_ymm_xmm;
}

// ---------------------------------------------------------------------------
// Runtime detection (x86_64)
// ---------------------------------------------------------------------------

/// Detect SSE + SSE2 support at runtime via CPUID leaf 1.
pub fn runtimeHasSSE() bool {
    if (comptime builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .x86)
        return false;
    const r = cpuid(1, 0);
    const sse_mask: u32 = (1 << 25); // EDX bit 25 = SSE
    const sse2_mask: u32 = (1 << 26); // EDX bit 26 = SSE2
    const mask = sse_mask | sse2_mask;
    return (r.edx & mask) == mask;
}

/// Detect SSE4.1 support at runtime (implies SSE + SSE2).
pub fn runtimeHasSSE4() bool {
    if (comptime builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .x86)
        return false;
    if (!runtimeHasSSE()) return false;
    const r = cpuid(1, 0);
    const sse41_mask: u32 = (1 << 19); // ECX bit 19 = SSE4.1
    return (r.ecx & sse41_mask) == sse41_mask;
}

/// Detect AVX + AVX2 support at runtime via CPUID + XGETBV.
/// Follows Easel's approach: checks OSXSAVE, YMM OS support,
/// and AVX2 in leaf 7.
pub fn runtimeHasAVX() bool {
    if (comptime builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .x86)
        return false;
    const r1 = cpuid(1, 0);
    // Check OSXSAVE (bit 27) — required for XGETBV
    const osxsave_mask: u32 = (1 << 27);
    if ((r1.ecx & osxsave_mask) != osxsave_mask) return false;
    // Check OS supports YMM state save
    if (!xgetbvChecksYmm()) return false;
    // Check AVX (bit 28 in ECX of leaf 1)
    const avx_mask: u32 = (1 << 28);
    return (r1.ecx & avx_mask) == avx_mask;
}

/// Detect AVX2 support at runtime.
/// Also checks FMA3, BMI1, and BMI2 which are required by AVX2 code paths
/// (matching Easel's esl_cpu_has_avx2 checks).
pub fn runtimeHasAVX2() bool {
    if (comptime builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .x86)
        return false;
    if (!runtimeHasAVX()) return false;
    // FMA3 = leaf 1, ECX bit 12
    const r1 = cpuid(1, 0);
    const fma_mask: u32 = (1 << 12);
    if ((r1.ecx & fma_mask) != fma_mask) return false;
    const r7 = cpuid(7, 0);
    // AVX2 = EBX bit 5, BMI1 = EBX bit 3, BMI2 = EBX bit 8
    const avx2_bmi_mask: u32 = (1 << 5) | (1 << 3) | (1 << 8);
    return (r7.ebx & avx2_bmi_mask) == avx2_bmi_mask;
}

/// Detect AVX-512 Foundation support at runtime.
/// Checks OSXSAVE, ZMM OS support, and AVX512F/DQ/BW in leaf 7.
pub fn runtimeHasAVX512() bool {
    if (comptime builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .x86)
        return false;
    const r1 = cpuid(1, 0);
    const osxsave_mask: u32 = (1 << 27);
    if ((r1.ecx & osxsave_mask) != osxsave_mask) return false;
    if (!xgetbvChecksZmm()) return false;
    const r7 = cpuid(7, 0);
    // AVX-512F = EBX bit 16, AVX-512DQ = bit 17, AVX-512BW = bit 30
    const avx512_mask: u32 = (1 << 16) | (1 << 17) | (1 << 30);
    return (r7.ebx & avx512_mask) == avx512_mask;
}

/// Detect NEON support at runtime.
/// On aarch64, NEON (Advanced SIMD) is always available.
pub fn runtimeHasNeon() bool {
    return comptime builtin.cpu.arch == .aarch64;
}

// ---------------------------------------------------------------------------
// Aggregate detection
// ---------------------------------------------------------------------------

/// Detect all SIMD features at runtime using CPUID (x86) or
/// architecture guarantees (aarch64).
pub fn runtimeDetect() SimdSupport {
    var support = SimdSupport{};

    switch (comptime builtin.cpu.arch) {
        .x86_64, .x86 => {
            support.sse = runtimeHasSSE();
            support.sse2 = runtimeHasSSE();
            support.sse4_1 = runtimeHasSSE4();
            support.avx = runtimeHasAVX();
            support.avx2 = runtimeHasAVX2();
            support.avx512f = runtimeHasAVX512();
        },
        .aarch64 => {
            support.neon = true;
        },
        else => {},
    }

    return support;
}

/// Detect available SIMD features using comptime target info.
/// (Zig cross-compilation model — the compiler already knows the target.)
pub fn comptimeDetect() SimdSupport {
    var support = SimdSupport{};

    switch (builtin.cpu.arch) {
        .x86_64, .x86 => {
            const features = builtin.cpu.features;
            support.sse = features.isEnabled(@intFromEnum(std.Target.x86.Feature.sse));
            support.sse2 = features.isEnabled(@intFromEnum(std.Target.x86.Feature.sse2));
            support.sse4_1 = features.isEnabled(@intFromEnum(std.Target.x86.Feature.sse4_1));
            support.avx = features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx));
            support.avx2 = features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx2));
            support.avx512f = features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx512f));
        },
        .aarch64 => {
            support.neon = true;
        },
        else => {},
    }

    return support;
}

/// Backwards-compatible alias: detect using comptime info.
pub const detect = comptimeDetect;

/// Return the number of available CPU cores.
pub fn cpuCount() usize {
    return std.Thread.getCpuCount() catch 1;
}

/// Return the best SIMD vector width in bytes for this platform.
/// Uses comptime detection (for backward compatibility).
pub fn bestVectorWidth() usize {
    return comptimeDetect().vectorWidth();
}

/// Return the best SIMD vector width in bytes detected at runtime.
pub fn runtimeBestVectorWidth() usize {
    return runtimeDetect().vectorWidth();
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "comptimeDetect: returns valid struct" {
    const simd = comptimeDetect();
    if (builtin.cpu.arch == .x86_64) {
        try std.testing.expect(simd.sse2);
    }
    if (builtin.cpu.arch == .aarch64) {
        try std.testing.expect(simd.neon);
    }
}

test "runtimeDetect: returns valid struct" {
    const simd = runtimeDetect();
    // On x86_64, at minimum SSE2 should be available (baseline since AMD64)
    if (builtin.cpu.arch == .x86_64) {
        try std.testing.expect(simd.sse);
        try std.testing.expect(simd.sse2);
    }
    // On aarch64, NEON should be available
    if (builtin.cpu.arch == .aarch64) {
        try std.testing.expect(simd.neon);
    }
}

test "runtimeDetect: consistency — comptime features are subset of runtime" {
    // If the compiler was told the target has a feature, the actual CPU
    // must also have it (assuming native target).
    if (builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .aarch64) return;

    const ct = comptimeDetect();
    const rt = runtimeDetect();

    // Every comptime-enabled feature should also be detected at runtime
    // (this holds when building for the native CPU).
    if (ct.sse) try std.testing.expect(rt.sse);
    if (ct.sse2) try std.testing.expect(rt.sse2);
    if (ct.sse4_1) try std.testing.expect(rt.sse4_1);
    if (ct.avx) try std.testing.expect(rt.avx);
    if (ct.avx2) try std.testing.expect(rt.avx2);
    if (ct.avx512f) try std.testing.expect(rt.avx512f);
    if (ct.neon) try std.testing.expect(rt.neon);
}

test "runtimeDetect: hierarchical consistency" {
    const rt = runtimeDetect();
    // AVX-512 implies AVX2 implies AVX implies SSE
    if (rt.avx512f) try std.testing.expect(rt.avx2);
    if (rt.avx2) try std.testing.expect(rt.avx);
    if (rt.avx) try std.testing.expect(rt.sse);
    if (rt.sse4_1) try std.testing.expect(rt.sse);
    if (rt.sse2) try std.testing.expect(rt.sse);
}

test "bestImplementation: returns known string" {
    const rt = runtimeDetect();
    const name = rt.bestImplementation();
    const valid = std.mem.eql(u8, name, "AVX512") or
        std.mem.eql(u8, name, "AVX") or
        std.mem.eql(u8, name, "SSE4") or
        std.mem.eql(u8, name, "SSE") or
        std.mem.eql(u8, name, "NEON") or
        std.mem.eql(u8, name, "none");
    try std.testing.expect(valid);
}

test "vectorWidth: reasonable values" {
    const rt = runtimeDetect();
    const w = rt.vectorWidth();
    try std.testing.expect(w == 8 or w == 16 or w == 32 or w == 64);
}

test "cpuCount: at least 1" {
    try std.testing.expect(cpuCount() >= 1);
}

test "bestVectorWidth: at least 8" {
    try std.testing.expect(bestVectorWidth() >= 8);
}

test "runtimeBestVectorWidth: at least 8" {
    try std.testing.expect(runtimeBestVectorWidth() >= 8);
}

test "individual runtime functions" {
    // Just exercise each function to ensure no crashes
    if (builtin.cpu.arch == .x86_64 or builtin.cpu.arch == .x86) {
        _ = runtimeHasSSE();
        _ = runtimeHasSSE4();
        _ = runtimeHasAVX();
        _ = runtimeHasAVX2();
        _ = runtimeHasAVX512();
    }
    _ = runtimeHasNeon();
}

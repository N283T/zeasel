// Seeded pseudo-random number generator wrapper.
//
// Wraps std.Random.Pcg to provide a reproducible, seedable RNG with
// utility methods matching Easel-style usage patterns.

const std = @import("std");
const math = std.math;

pub const Random = struct {
    pcg: std.Random.Pcg,

    /// Create a new RNG with a specific seed.
    pub fn init(seed: u64) Random {
        return Random{ .pcg = std.Random.Pcg.init(seed) };
    }

    /// Return the underlying std.Random interface.
    fn rng(self: *Random) std.Random {
        return self.pcg.random();
    }

    /// Uniform random float in [0, 1).
    pub fn uniform(self: *Random) f64 {
        return self.pcg.random().float(f64);
    }

    /// Uniform random integer in [0, n).
    pub fn uniformInt(self: *Random, n: u32) u32 {
        return self.pcg.random().intRangeLessThan(u32, 0, n);
    }

    /// Gaussian (normal) random variate with mean 0, standard deviation 1.
    /// Uses the Box-Muller transform. Consumes two uniform samples per call.
    pub fn gaussian(self: *Random) f64 {
        // Box-Muller: given u1, u2 uniform in (0, 1),
        // z = sqrt(-2 * ln(u1)) * cos(2*pi*u2) is standard normal.
        // Guard against u1 == 0 to avoid ln(0).
        var s: f64 = 0.0;
        while (s == 0.0) {
            s = self.uniform();
        }
        const t = self.uniform();
        const mag = math.sqrt(-2.0 * math.log(f64, math.e, s));
        return mag * math.cos(2.0 * math.pi * t);
    }

    /// Random choice from a probability distribution p[0..n-1].
    /// p must sum to approximately 1.0 and all values must be non-negative.
    /// Returns the index of the chosen residue.
    pub fn choose(self: *Random, p: []const f64) usize {
        const u = self.uniform();
        var cumulative: f64 = 0.0;
        for (p, 0..) |prob, i| {
            cumulative += prob;
            if (u < cumulative) return i;
        }
        // Return last index if floating-point rounding pushes u >= sum(p).
        return p.len - 1;
    }
};

// --- Tests ---

test "init: reproducibility — same seed produces same sequence" {
    var r1 = Random.init(42);
    var r2 = Random.init(42);

    for (0..20) |_| {
        try std.testing.expectEqual(r1.uniform(), r2.uniform());
    }
}

test "init: different seeds produce different sequences" {
    var r1 = Random.init(1);
    var r2 = Random.init(2);

    // Extremely unlikely all 10 values match by chance.
    var all_equal = true;
    for (0..10) |_| {
        if (r1.uniform() != r2.uniform()) {
            all_equal = false;
            break;
        }
    }
    try std.testing.expect(!all_equal);
}

test "uniform: all values in [0, 1)" {
    var rng = Random.init(123);
    for (0..1000) |_| {
        const v = rng.uniform();
        try std.testing.expect(v >= 0.0);
        try std.testing.expect(v < 1.0);
    }
}

test "uniformInt: all values in [0, n)" {
    var rng = Random.init(456);
    const n: u32 = 7;
    for (0..1000) |_| {
        const v = rng.uniformInt(n);
        try std.testing.expect(v < n);
    }
}

test "uniformInt: n=1 always returns 0" {
    var rng = Random.init(789);
    for (0..100) |_| {
        try std.testing.expectEqual(@as(u32, 0), rng.uniformInt(1));
    }
}

test "choose: uniform distribution — all 4 choices appear over 1000 draws" {
    var rng = Random.init(999);
    const freq = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    var counts = [_]usize{0} ** 4;

    for (0..1000) |_| {
        const idx = rng.choose(&freq);
        counts[idx] += 1;
    }

    for (counts) |c| {
        try std.testing.expect(c > 0);
    }
}

test "choose: deterministic with biased distribution" {
    var rng = Random.init(0);
    // Always pick index 2 — probability 1.0
    const freq = [_]f64{ 0.0, 0.0, 1.0, 0.0 };
    for (0..100) |_| {
        try std.testing.expectEqual(@as(usize, 2), rng.choose(&freq));
    }
}

test "gaussian: mean near 0 and std near 1 over 1000 samples" {
    var rng = Random.init(314159);
    var sum: f64 = 0.0;
    var sum_sq: f64 = 0.0;
    const n: f64 = 1000.0;

    for (0..@intFromFloat(n)) |_| {
        const v = rng.gaussian();
        sum += v;
        sum_sq += v * v;
    }

    const mean = sum / n;
    const variance = sum_sq / n - mean * mean;
    const std_dev = @sqrt(variance);

    // With 1000 samples the CLT guarantees these tolerances hold with very high
    // probability (mean within 0.15, std_dev within 0.15 of 1.0).
    try std.testing.expect(@abs(mean) < 0.15);
    try std.testing.expect(@abs(std_dev - 1.0) < 0.15);
}

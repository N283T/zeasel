/// Statistics and probability distributions used by HMMER for E-value calculations.
pub const gumbel = @import("stats/gumbel.zig");
pub const exponential = @import("stats/exponential.zig");
pub const gamma = @import("stats/gamma.zig");
pub const dirichlet = @import("stats/dirichlet.zig");
pub const minimizer = @import("stats/minimizer.zig");
pub const rootfinder = @import("stats/rootfinder.zig");

test {
    _ = gumbel;
    _ = exponential;
    _ = gamma;
    _ = dirichlet;
    _ = minimizer;
    _ = rootfinder;
}

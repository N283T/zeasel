// I/O subsystem: parsers and writers for biological sequence file formats.

pub const fasta = @import("io/fasta.zig");
pub const stockholm = @import("io/stockholm.zig");
pub const genbank = @import("io/genbank.zig");
pub const clustal = @import("io/clustal.zig");
pub const afa = @import("io/afa.zig");
pub const phylip = @import("io/phylip.zig");
pub const a2m = @import("io/a2m.zig");
pub const psiblast = @import("io/psiblast.zig");
pub const selex = @import("io/selex.zig");
pub const Reader = @import("io/reader.zig").Reader;
pub const Writer = @import("io/writer.zig").Writer;
pub const Format = @import("io/reader.zig").Format;

const reader_mod = @import("io/reader.zig");
const writer_mod = @import("io/writer.zig");
const stockholm_mod = @import("io/stockholm.zig");
const genbank_mod = @import("io/genbank.zig");
const clustal_mod = @import("io/clustal.zig");
const afa_mod = @import("io/afa.zig");

// Include tests from all I/O sub-modules in `zig build test`.
test {
    _ = fasta;
    _ = stockholm_mod;
    _ = reader_mod;
    _ = writer_mod;
    _ = genbank_mod;
    _ = clustal_mod;
    _ = afa_mod;
    _ = @import("io/phylip.zig");
    _ = @import("io/a2m.zig");
    _ = @import("io/psiblast.zig");
    _ = @import("io/selex.zig");
}

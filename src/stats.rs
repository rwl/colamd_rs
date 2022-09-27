// Knob and statistics definitions. //

/// Size of the knobs[] array. Only knobs[0..1] are currently used.
pub(crate) const KNOBS: usize = 20;

/// Number of output statistics. Only stats[0..6] are currently used.
pub const STATS: usize = 20;

/// Dense row knob and output statistic.
pub(crate) const DENSE_ROW: usize = 0;

/// Dense column knob and output statistic.
pub(crate) const DENSE_COL: usize = 1;

/// Aggressive absorption.
pub(crate) const AGGRESSIVE: usize = 2;

/// Memory defragmentation count output statistic.
pub(crate) const DEFRAG_COUNT: usize = 2;

/// Zero OK, > 0 warning or notice, < 0 error.
pub(crate) const STATUS: usize = 3;

// Error info, or info on jumbled columns.
pub(crate) const INFO1: usize = 4;
pub(crate) const INFO2: usize = 5;
pub(crate) const INFO3: usize = 6;

pub(crate) fn default_knobs() -> [f64; KNOBS] {
    let mut knobs = [0.0; KNOBS];
    for i in 0..KNOBS {
        knobs[i] = 0.0
    }
    knobs[DENSE_ROW] = 10.0;
    knobs[DENSE_COL] = 10.0;
    knobs[AGGRESSIVE] = 1.0; // Default to aggressive absorption.
    knobs
}

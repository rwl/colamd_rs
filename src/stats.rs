// Knob and statistics definitions. //

/// Size of the `knobs[]` array. Only `knobs[0..1]` are currently used.
pub const KNOBS: usize = 20;

/// Number of output statistics. Only `stats[0..6]` are currently used.
pub const STATS: usize = 20;

/// Dense row knob and output statistic.
pub const DENSE_ROW: usize = 0;

/// Dense column knob and output statistic.
pub const DENSE_COL: usize = 1;

/// Aggressive absorption.
pub const AGGRESSIVE: usize = 2;

/// Memory defragmentation count output statistic.
pub const DEFRAG_COUNT: usize = 2;

/// Zero OK, > 0 warning or notice, < 0 error.
pub const STATUS: usize = 3;

/// Error info, or info on jumbled columns.
pub const INFO1: usize = 4;
pub const INFO2: usize = 5;
pub const INFO3: usize = 6;

// Error codes returned in `stats[3]`.
pub const OK: i32 = 0;
pub const OK_BUT_JUMBLED: i32 = 1;
// pub const ERROR_A_NOT_PRESENT: i32 = -1;
// pub const ERROR_P_NOT_PRESENT: i32 = -2;
pub const ERROR_NROW_NEGATIVE: i32 = -3;
pub const ERROR_NCOL_NEGATIVE: i32 = -4;
pub const ERROR_NNZ_NEGATIVE: i32 = -5;
pub const ERROR_P0_NONZERO: i32 = -6;
pub const ERROR_A_TOO_SMALL: i32 = -7;
pub const ERROR_COL_LENGTH_NEGATIVE: i32 = -8;
pub const ERROR_ROW_INDEX_OUT_OF_BOUNDS: i32 = -9;
pub const ERROR_OUT_OF_MEMORY: i32 = -10;

/// Returns the default values of the user-controllable parameters for
/// `colamd` and `symamd`.
///
/// Colamd: rows with more than `max(16, knobs[0] * sqrt(n_col))`
/// entries are removed prior to ordering. Columns with more than
/// `max(16, knobs[1] * sqrt(min(n_row, n_col)))` entries are removed
/// prior to ordering, and placed last in the output column ordering.
///
/// Symamd: Rows and columns with more than `max(16, knobs[0] * sqrt(n))`
/// entries are removed prior to ordering, and placed last in the
/// output ordering.
///
/// ```txt
/// knobs[0]  dense row control
/// knobs[1]  dense column control
/// knobs[2]  if nonzero, do aggresive absorption
/// knobs[3..19]  unused, but future versions might use this
/// knobs [0] and knobs [1] control dense row and col detection:
/// ```
///
/// Colamd: rows with more than `max(16, knobs[DENSE_ROW] * sqrt(n_col))`
/// entries are removed prior to ordering. Columns with more than
/// `max(16, knobs[DENSE_COL] * sqrt(min(n_row,n_col)))` entries are removed
/// prior to ordering, and placed last in the output column ordering.
///
/// Symamd: uses only `knobs[DENSE_ROW]`, which is `knobs[0]`.
/// Rows and columns with more than `max(16, knobs[DENSE_ROW] * sqrt(n))`
/// entries are removed prior to ordering, and placed last in the
/// output ordering.
///
/// `DENSE_ROW` and `DENSE_COL` are defined as 0 and 1, respectively. Default
/// values of these two knobs are both 10. Currently, only `knobs[0]` and
/// `knobs[1]` are used, but future versions may use more knobs. If so, they
/// will be properly set to their defaults by the future version of
/// `default_knobs`, so that the code that calls `colamd` will not need to
/// change, assuming that you either use `default_knobs`, or pass `None` as
/// the knobs array to colamd or symamd.
///
/// ```txt
/// knobs[2]: aggressive absorption
/// knobs[AGGRESSIVE] controls whether or not to do
/// aggressive absorption during the ordering. Default is TRUE.
/// ```
pub fn default_knobs() -> [f64; KNOBS] {
    let mut knobs = [0.0; KNOBS];
    for i in 0..KNOBS {
        knobs[i] = 0.0
    }
    knobs[DENSE_ROW] = 10.0;
    knobs[DENSE_COL] = 10.0;
    knobs[AGGRESSIVE] = 1.0; // Default to aggressive absorption.
    knobs
}

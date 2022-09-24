// Error codes returned in stats[3].
pub(crate) const OK: i32 = 0;
pub(crate) const OK_BUT_JUMBLED: i32 = 1;
// pub(crate) const ERROR_A_NOT_PRESENT: i32 = -1;
// pub(crate) const ERROR_P_NOT_PRESENT: i32 = -2;
pub(crate) const ERROR_NROW_NEGATIVE: i32 = -3;
pub(crate) const ERROR_NCOL_NEGATIVE: i32 = -4;
pub(crate) const ERROR_NNZ_NEGATIVE: i32 = -5;
pub(crate) const ERROR_P0_NONZERO: i32 = -6;
pub(crate) const ERROR_A_TOO_SMALL: i32 = -7;
pub(crate) const ERROR_COL_LENGTH_NEGATIVE: i32 = -8;
pub(crate) const ERROR_ROW_INDEX_OUT_OF_BOUNDS: i32 = -9;
pub(crate) const ERROR_OUT_OF_MEMORY: i32 = -10;

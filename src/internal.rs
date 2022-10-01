use crate::col::Col;
use crate::row::{Row, RowShared2};

#[cfg(not(feature = "i64"))]
pub type Int = i32;

#[cfg(feature = "i64")]
pub type Int = i64;

pub(crate) fn dense_degree(alpha: f64, n: Int) -> Int {
    f64::max(16.0, alpha * (n as f64).sqrt()) as Int
}

pub(crate) fn ones_complement(r: Int) -> Int {
    return -r - 1;
}

pub(crate) const EMPTY: Int = -1;

// Row and column status.
pub(crate) const ALIVE: Int = 0;
pub(crate) const DEAD: Int = -1;

// Column status.
pub(crate) const DEAD_PRINCIPAL: Int = -1;
pub(crate) const DEAD_NON_PRINCIPAL: Int = -2;

// Row and column status update and checking.
pub(crate) fn row_is_dead(row: &[Row], r: usize) -> bool {
    row_is_marked_dead(row[r].shared2.mark())
}

pub(crate) fn row_is_marked_dead(row_mark: Int) -> bool {
    row_mark < ALIVE
}

pub(crate) fn row_is_alive(row: &[Row], r: usize) -> bool {
    row[r].shared2.mark() >= ALIVE
}

pub(crate) fn col_is_dead(col: &[Col], c: usize) -> bool {
    col[c].start < ALIVE
}

pub(crate) fn col_is_alive(col: &[Col], c: usize) -> bool {
    col[c].start >= ALIVE
}

pub(crate) fn col_is_dead_principal(col: &mut [Col], c: usize) -> bool {
    col[c].start == DEAD_PRINCIPAL
}

pub(crate) fn kill_row(row: &mut [Row], r: usize) {
    row[r].shared2 = RowShared2::Mark(DEAD)
}

pub(crate) fn kill_principal_col(col: &mut [Col], c: usize) {
    col[c].start = DEAD_PRINCIPAL
}

pub(crate) fn kill_non_principal_col(col: &mut [Col], c: usize) {
    col[c].start = DEAD_NON_PRINCIPAL
}

// Feature: debug

#[cfg(feature = "debug")]
macro_rules! assert_debug {
    ($cond:expr) => {
        assert!($cond)
    };
}

#[cfg(not(feature = "debug"))]
macro_rules! assert_debug {
    ($cond:expr) => {};
}

#[cfg(feature = "debug")]
macro_rules! debug0 {
    ($( $args:expr ),*) => {
        #[cfg(not(feature = "nprint"))]
        println!( $( $args ),* );
    }
}

#[cfg(not(feature = "debug"))]
macro_rules! debug0 {
    ($( $args:expr ),*) => {};
}

// Feature: debug1

#[cfg(feature = "debug1")]
macro_rules! debug1 {
    ($( $args:expr ),*) => {
        #[cfg(not(feature = "nprint"))]
        println!( $( $args ),* );
    }
}

#[cfg(not(feature = "debug1"))]
macro_rules! debug1 {
    ($( $args:expr ),*) => {};
}

// Feature: debug2

#[cfg(feature = "debug2")]
macro_rules! debug2 {
    ($( $args:expr ),*) => {
        #[cfg(not(feature = "nprint"))]
        println!( $( $args ),* );
    }
}

#[cfg(not(feature = "debug2"))]
macro_rules! debug2 {
    ($( $args:expr ),*) => {};
}

// Feature: debug3

#[cfg(feature = "debug3")]
macro_rules! debug3 {
    ($( $args:expr ),*) => {
        #[cfg(not(feature = "nprint"))]
        println!( $( $args ),* );
    }
}

#[cfg(not(feature = "debug3"))]
macro_rules! debug3 {
    ($( $args:expr ),*) => {};
}

// Feature: debug4

#[cfg(feature = "debug4")]
macro_rules! debug4 {
    ($( $args:expr ),*) => {
        #[cfg(not(feature = "nprint"))]
        println!( $( $args ),* );
    }
}

#[cfg(not(feature = "debug4"))]
macro_rules! debug4 {
    ($( $args:expr ),*) => {};
}

pub(crate) use {assert_debug, debug0, debug1, debug2, debug3, debug4};

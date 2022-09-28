pub(crate) fn dense_degree(alpha: f64, n: i32) -> i32 {
    f64::max(16.0, alpha * (n as f64).sqrt()) as i32
}

pub(crate) fn ones_complement(r: i32) -> i32 {
    return -r - 1;
}

pub(crate) const EMPTY: i32 = -1;

// Row and column status.
pub(crate) const ALIVE: i32 = 0;
pub(crate) const DEAD: i32 = -1;

// Column status.
pub(crate) const DEAD_PRINCIPAL: i32 = -1;
pub(crate) const DEAD_NON_PRINCIPAL: i32 = -2;

// Row and column status update and checking.
pub(crate) fn row_is_dead(row: &[Row], r: usize) -> bool {
    row_is_marked_dead(row[r].mark())
}

pub(crate) fn row_is_marked_dead(row_mark: i32) -> bool {
    row_mark < ALIVE
}

pub(crate) fn row_is_alive(row: &[Row], r: usize) -> bool {
    row[r].mark() >= ALIVE
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
    row[r].set_mark(DEAD)
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

use crate::col::Col;
use crate::row::Row;
pub(crate) use {assert_debug, debug0, debug1, debug2, debug3, debug4};

// Copyright (c) 1998-2007, Timothy A. Davis, All Rights Reserved.

use crate::codes::*;
use crate::col::Col;
use crate::colamd2::{init_rows_cols, init_scoring};
use crate::internal::*;
use crate::row::Row;
use crate::stats::*;
use std::cmp;

pub fn colamd(
    n_row: i32,
    n_col: i32,
    Alen: i32,
    A: &[i32],
    p: &mut [i32],
    knobs: Option<[f64; KNOBS]>,
    stats: &mut [i32; STATS],
) -> bool {
    let mut nnz: i32; // nonzeros in A
    let mut Row_size: i32; // size of Row[], in integers
    let mut Col_size: i32; // size of Col[], in integers
    let mut need: i32; // minimum required length of A
    let mut rows: Vec<Row>; // pointer into A of Row[0..n_row] array
    let mut cols: Vec<Col>; // pointer into A of Col[0..n_col] array
    let mut n_col2: i32; // number of non-dense, non-empty columns
    let mut n_row2: i32; // number of non-dense, non-empty rows
    let mut ngarbage: i32; // number of garbage collections performed
    let mut max_deg: i32; // maximum row degree
                          // let mut default_knobs = vec![0.0; KNOBS]; // default knobs array
    let mut aggressive: i32; // do aggressive absorption
    let mut ok = false;

    // Check the input arguments.
    // #[cfg(feature = "debug")]
    // if stats.is_none() {
    //     debug0!("colamd: stats not present");
    //     return false;
    // }
    for i in 0..STATS {
        stats[i] = 0;
    }
    stats[STATUS] = OK;
    stats[INFO1] = -1;
    stats[INFO2] = -1;

    // if A.is_none() {
    //     // A is not present.
    //     stats[STATUS] = ERROR_A_NOT_PRESENT;
    //     debug0!("colamd: A not present");
    //     return false;
    // }
    //
    // if p.is_none() {
    //     // p is not present.
    //     stats[STATUS] = ERROR_P_NOT_PRESENT;
    //     debug0!("colamd: p not present");
    //     return false;
    // }

    if n_row < 0 {
        // n_row must be >= 0
        stats[STATUS] = ERROR_NROW_NEGATIVE;
        stats[INFO1] = n_row;
        debug0!("colamd: nrow negative {}", n_row);
        return false;
    }

    if n_col < 0 {
        // n_col must be >= 0
        stats[STATUS] = ERROR_NCOL_NEGATIVE;
        stats[INFO1] = n_col;
        debug0!("colamd: ncol negative {}", n_col);
        return false;
    }

    let nnz = p[n_col as usize];
    if nnz < 0 {
        // nnz must be >= 0
        stats[STATUS] = ERROR_NNZ_NEGATIVE;
        stats[INFO1] = nnz;
        debug0!("colamd: number of entries negative {}", nnz);
        return false;
    }

    if p[0] != 0 {
        stats[STATUS] = ERROR_P0_NONZERO;
        stats[INFO1] = p[0];
        debug0!("colamd: p[0] not zero {}", p[0]);
        return false;
    }

    // If no knobs, set default knobs.
    let knobs = if knobs.is_none() {
        default_knobs()
    } else {
        knobs.unwrap()
    };

    let aggressive = knobs[AGGRESSIVE] != 0.0;

    // Allocate the Row and Col arrays from array A.
    ok = true;
    Col_size = tc(n_col, &mut ok); // Size of Col array of structs
    Row_size = tr(n_row, &mut ok); // Size of Row array of structs

    // need = 2*nnz + n_col + Col_size + Row_size
    need = tmult(nnz, 2, &mut ok);
    need = tadd(need, n_col, &mut ok);
    // need = t_add (need, Col_size, ok)
    // need = t_add (need, Row_size, ok)

    if !ok || need > Alen || need > i32::MAX {
        // Not enough space in array A to perform the ordering.
        stats[STATUS] = ERROR_A_TOO_SMALL;
        stats[INFO1] = need;
        stats[INFO2] = Alen;
        debug0!("colamd: Need Alen >= {}, given only Alen = {}", need, Alen);
        return false;
    }

    // Alen -= Col_size + Row_size
    let mut cols = vec![Col::default(); Col_size as usize]; //A[Alen]
    let mut rows = vec![Row::default(); Row_size as usize]; //A[Alen + Col_size]

    // Construct the row and column data structures.

    if !init_rows_cols(n_row, n_col, &mut rows, &mut cols, A, p, stats) {
        // Input matrix is invalid.
        debug0!("colamd: Matrix invalid");
        return false;
    }

    // Initialize scores, kill dense rows/columns.
    init_scoring(
        n_row,
        n_col,
        &rows,
        &cols,
        A,
        p,
        knobs,
        &mut n_row2,
        &mut n_col2,
        &mut max_deg,
    );

    // Order the supercolumns.
    // let ngarbage = find_ordering(
    //     n_row,
    //     n_col,
    //     Alen,
    //     rows,
    //     cols,
    //     A,
    //     p,
    //     n_col2,
    //     max_deg,
    //     2 * nnz,
    //     aggressive,
    // );

    // Order the non-principal columns.
    // order_children(n_col, cols, p);

    // Return statistics in stats.
    stats[DENSE_ROW] = n_row - n_row2;
    stats[DENSE_COL] = n_col - n_col2;
    stats[DEFRAG_COUNT] = ngarbage;

    debug0!("colamd: done.");

    true
}

fn default_knobs() -> [f64; KNOBS] {
    let mut knobs = [0.0; KNOBS];
    for i in 0..KNOBS {
        knobs[i] = 0.0
    }
    knobs[DENSE_ROW] = 10.0;
    knobs[DENSE_COL] = 10.0;
    knobs[AGGRESSIVE] = 1.0; // Default to aggressive absorption.
    knobs
}

// Add two values of type int, and check for integer overflow.
fn tadd(a: i32, b: i32, ok: &mut bool) -> i32 {
    if (*ok != false) && ((a + b) >= cmp::max(a, b)) {
        *ok = true
    } else {
        *ok = false
    }
    if *ok {
        a + b
    } else {
        0
    }
}

// Compute a*k where k is a small integer, and check for integer overflow.
fn tmult(a: i32, k: i32, ok: &mut bool) -> i32 {
    let mut s = 0;
    for i in 0..k {
        s = tadd(s, a, ok);
    }
    s
}

// Size of the Col and Row structures.

fn tc(n_col: i32, ok: &mut bool) -> i32 {
    // return ((t_mult (t_add (n_col, 1, ok), sizeof (Colamd_Col), ok) / sizeof (int))) ;
    tadd(n_col, 1, ok)
}

fn tr(n_row: i32, ok: &mut bool) -> i32 {
    // return ((t_mult (t_add (n_row, 1, ok), sizeof (Colamd_Row), ok) / sizeof (int))) ;
    tadd(n_row, 1, ok)
}

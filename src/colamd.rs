use crate::col::Col;
use crate::colamd2::{find_ordering, init_rows_cols, init_scoring, order_children};
use crate::internal::*;
use crate::row::Row;
use crate::stats::*;
use std::cmp;

/// Computes a column ordering `Q` of a sparse matrix
/// `A` such that the LU factorization `P(AQ) = LU` remains sparse,
/// where `P` is selected via partial pivoting. The routine can also
/// be viewed as providing a permutation `Q` such that the Cholesky
/// factorization `(AQ)'(AQ) = LL'` remains sparse.
///
/// Computes a column ordering (`Q`) of `A` such that `P(AQ)=LU` or
/// `(AQ)'AQ=LL'` have less fill-in and require fewer floating point
/// operations than factorizing the unpermuted matrix `A` or `A'A`,
/// respectively.
///
/// `colamd` returns `false` if `n_row` or `n_col` are negative.
/// `a_len >= 2*nnz + 6*(n_col+1) + 4*(n_row+1) + n_col`.
/// `colamd` returns `false` if these conditions are not met.
///
/// Note: this restriction makes an modest assumption regarding
/// the size of the two structures. We do, however, guarantee that
///
///     a_len >= recommended(nnz, n_row, n_col)
///
/// will be sufficient.
///
/// `a_i` is an integer array of size `a_len`. `a_len` must be at least
/// as large as the bare minimum value given above, but this is very
/// low, and can result in excessive run time. For best
/// performance, we recommend that `a_len` be greater than or equal to
/// `recommended(nnz, n_row, n_col)`, which adds `nnz/5` to the bare
/// minimum value given above.
///
/// On input, the row indices of the entries in column `c` of the
/// matrix are held in `a_i[(p[c]) ... (p[c+1]-1)]`. The row indices
/// in a given column `c` need not be in ascending order, and
/// duplicate row indices may be be present. However, `colamd` will
/// work a little faster if both of these conditions are met
/// (`colamd` puts the matrix into this format, if it finds that the
/// the conditions are not met).
///
/// The matrix is 0-based. That is, rows are in the range 0 to
/// `n_row-1`, and columns are in the range 0 to `n_col-1`. Colamd
/// returns `false` if any row index is out of range.
///
/// The contents of `a_i` are modified during ordering, and are
/// undefined on output.
///
/// `p` is an integer array of size `n_col+1`. On input, it holds the
/// "pointers" for the column form of the matrix `A`. Column `c` of
/// the matrix `A` is held in `a_i[(p [c]) ... (p [c+1]-1)]`. The first
/// entry, `p[0]`, must be zero, and `p[c] <= p[c+1]` must hold
/// for all `c` in the range 0 to `n_col-1`. The value `p[n_col]` is
/// thus the total number of entries in the pattern of the matrix `A`.
/// `colamd` returns `false` if these conditions are not met.
///
/// On output, if `colamd` returns `true`, the array `p` holds the column
/// permutation (`Q`, for `P(AQ)=LU` or `(AQ)'(AQ)=LL'`), where `p[0]` is
/// the first column index in the new ordering, and `p[n_col-1]` is
/// the last. That is, `p[k] = j` means that column `j` of `A` is the
/// `k`th pivot column, in `AQ`, where `k` is in the range 0 to `n_col-1`
/// (`p[0] = j` means that column `j` of `A` is the first column in `AQ`).
///
/// If `colamd` returns `false`, then no permutation is returned, and
/// `p` is undefined on output.
///
/// Statistics on the ordering, and error status:
///
/// ```txt
/// stats[0]: number of dense or empty rows ignored.
/// stats[1]: number of dense or empty columns ignored (and
///     ordered last in the output permutation p)
///     Note that a row can become "empty" if it
///     contains only "dense" and/or "empty" columns,
///     and similarly a column can become "empty" if it
///     only contains "dense" and/or "empty" rows.
/// stats[2]: number of garbage collections performed.
///     This can be excessively high if a_len is close
///     to the minimum required value.
/// stats[3]: status code. < 0 is an error code.
///       > 1 is a warning or notice.
///     0: OK. Each column of the input matrix contained
///       row indices in increasing order, with no
///       duplicates.
///     1: OK, but columns of input matrix were jumbled
///       (unsorted columns or duplicate entries). Colamd
///       had to do some extra work to sort the matrix
///       first and remove duplicate entries, but it
///       still was able to return a valid permutation
///       (return value of colamd was TRUE).
///         stats[4]: highest numbered column that
///           is unsorted or has duplicate
///           entries.
///         stats[5]: last seen duplicate or
///           unsorted row index.
///         stats[6]: number of duplicate or
///           unsorted row indices.
///     -1: A is a null pointer
///     -2: p is a null pointer
///     -3: n_row is negative
///         stats[4]: n_row
///     -4: n_col is negative
///         stats[4]: n_col
///     -5: number of nonzeros in matrix is negative
///         stats[4]: number of nonzeros, p [n_col]
///     -6: p[0] is nonzero
///         stats[4]: p[0]
///     -7: A is too small
///         stats[4]: required size
///         stats[5]: actual size (a_len)
///     -8: a column has a negative number of entries
///         stats[4]: column with < 0 entries
///         stats[5]: number of entries in col
///     -9: a row index is out of bounds
///         stats[4]: column with bad row index
///         stats[5]: bad row index
///         stats[6]: n_row, # of rows of matrx
///     -10: (unused; see symamd.go)
///     -999  (unused; see symamd.go)
/// ```
pub fn colamd(
    n_row: i32,
    n_col: i32,
    a_len: i32,
    a_i: &mut [i32],
    p: &mut [i32],
    knobs: Option<[f64; KNOBS]>,
    stats: &mut [i32; STATS],
) -> bool {
    for i in 0..STATS {
        stats[i] = 0;
    }
    stats[STATUS] = OK;
    stats[INFO1] = -1;
    stats[INFO2] = -1;

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
    let mut ok = true;
    let col_size = tc(n_col, &mut ok); // Size of Col array of structs
    let row_size = tr(n_row, &mut ok); // Size of Row array of structs

    // need = 2*nnz + n_col + Col_size + Row_size
    let need = tmult(nnz, 2, &mut ok);
    let need = tadd(need, n_col, &mut ok);
    // need = t_add (need, Col_size, ok)
    // need = t_add (need, Row_size, ok)

    if !ok || need > a_len || need > i32::MAX {
        // Not enough space in array A to perform the ordering.
        stats[STATUS] = ERROR_A_TOO_SMALL;
        stats[INFO1] = need;
        stats[INFO2] = a_len;
        debug0!(
            "colamd: Need a_len >= {}, given only a_len = {}",
            need,
            a_len
        );
        return false;
    }

    // a_len -= Col_size + Row_size
    let mut cols = vec![Col::default(); col_size as usize]; //A[a_len]
    let mut rows = vec![Row::default(); row_size as usize]; //A[a_len + Col_size]

    // Construct the row and column data structures.

    if !init_rows_cols(n_row, n_col, &mut rows, &mut cols, a_i, p, stats) {
        // Input matrix is invalid.
        debug0!("colamd: Matrix invalid");
        return false;
    }

    // Initialize scores, kill dense rows/columns.
    let mut n_col2: i32 = 0; // number of non-dense, non-empty columns
    let mut n_row2: i32 = 0; // number of non-dense, non-empty rows
    let mut max_deg: i32 = 0; // maximum row degree
    init_scoring(
        n_row,
        n_col,
        &mut rows,
        &mut cols,
        a_i,
        p,
        knobs,
        &mut n_row2,
        &mut n_col2,
        &mut max_deg,
    );

    // Order the supercolumns.
    let ngarbage = find_ordering(
        n_row,
        n_col,
        a_len,
        &mut rows,
        &mut cols,
        a_i,
        p,
        n_col2,
        max_deg,
        2 * nnz,
        aggressive,
    );

    // Order the non-principal columns.
    order_children(n_col, &mut cols, p);

    // Return statistics in stats.
    stats[DENSE_ROW] = n_row - n_row2;
    stats[DENSE_COL] = n_col - n_col2;
    stats[DEFRAG_COUNT] = ngarbage;

    debug0!("colamd: done.");

    true
}

/// Recommended returns the suggested size for `a_len`.
///
/// This value has been determined to provide good balance between the number
/// of garbage collections and the memory requirements for `colamd`. If any
/// argument is negative, or if integer overflow occurs, a 0 is returned as
/// an error condition. `2*nnz` space is required for the row and column
/// indices of the matrix. `COLAMD_C(n_col) + COLAMD_R(n_row)` space is
/// required for the `Col` and `Row` arrays, respectively, which are internal to
/// `colamd` (roughly `6*n_col + 4*n_row`). An additional `n_col` space is the
/// minimal amount of "elbow room", and `nnz/5` more space is recommended for
/// run time efficiency.
///
/// `a_len` is approximately `2.2*nnz + 7*n_col + 4*n_row + 10`.
///
/// This function is not needed when using `symamd`.
///
/// `nnz` must be the same value as `p[n_col]` in the call to `colamd` - otherwise
/// you will get a wrong value of the recommended memory to use.
/// Returns the recommended value of `a_len` or 0 if any input argument is
/// negative. The use of this routine is optional. Not needed for `symamd`,
/// which dynamically allocates its own memory.
pub fn recommended(nnz: i32, n_row: i32, n_col: i32) -> i32 {
    let mut ok = true;
    if nnz < 0 || n_row < 0 || n_col < 0 {
        return 0;
    }
    let mut s = tmult(nnz, 2, &mut ok); // 2*nnz

    //c = COLAMD_C(n_col, ok) ;      // size of column structures
    //r = COLAMD_R(n_row, ok) ;      // size of row structures
    //s = t_add(s, c, ok) ;
    //s = t_add(s, r, ok) ;

    s = tadd(s, n_col, &mut ok); // elbow room
    s = tadd(s, nnz / 5, &mut ok); // elbow room
    ok = true; //(s < Int_MAX) ? 1 : 0;
    if ok {
        s
    } else {
        0
    }
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
    for _i in 0..k {
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

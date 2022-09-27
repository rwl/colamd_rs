use crate::codes::*;
use crate::colamd;
use crate::colamd::recommended;
use crate::internal::*;
use crate::stats::*;

/// Computes an ordering `P` of a symmetric sparse
/// matrix `A` such that the Cholesky factorization `PAP' = LL'` remains
/// sparse. It is based on a column ordering of a matrix `M` constructed
/// so that the nonzero pattern of `M'M` is the same as `A`. The matrix `A`
/// is assumed to be symmetric; only the strictly lower triangular part
/// is accessed.
///
/// The row indices of the entries in column `c` of the matrix are
/// held in `A[(p[c])...(p[c+1]-1)]`. The row indices in a
/// given column `c` need not be in ascending order, and duplicate
/// row indices may be present. However, `symamd` will run faster
/// if the columns are in sorted order with no duplicate entries.
///
/// The matrix is 0-based. That is, rows are in the range 0 to
/// `n-1`, and columns are in the range 0 to `n-1`. `symamd`
/// returns `false` if any row index is out of range.
///
/// The contents of `A` are not modified.
/// `A` is an integer array of size `n+1`. On input, it holds the
/// "pointers" for the column form of the matrix `A`. Column `c` of
/// the matrix `A` is held in `A[(p[c])...(p[c+1]-1)]`. The first
/// entry, `p[0]`, must be zero, and `p[c] <= p[c+1]` must hold
/// for all `c` in the range 0 to `n-1`. The value `p[n]` is
/// thus the total number of entries in the pattern of the matrix `A`.
/// `symamd` returns `false` if these conditions are not met.
///
/// The contents of `p` are not modified.
/// On output, if `symamd` returns `true`, the array perm holds the
/// permutation `P`, where `perm[0]` is the first index in the new
/// ordering, and `perm[n-1]` is the last. That is, `perm[k] = j`
/// means that row and column `j` of `A` is the `k`th column in `PAP'`,
/// where `k` is in the range 0 to `n-1` (`perm[0] = j` means
/// that row and column `j` of `A` are the first row and column in
/// `PAP'`). The array is used as a workspace during the ordering,
/// which is why it must be of length `n+1`, not just `n`.
///
///     `stats[0]`: number of dense or empty row and columns ignored
///        (and ordered last in the output permutation
///        perm). Note that a row/column can become
///        "empty" if it contains only "dense" and/or
///        "empty" columns/rows.
///     `stats[1]`: (same as `stats[0]`)
///     `stats[2]`: number of garbage collections performed.
///     `stats[3]`: status code < 0 is an error code.
///          > 1 is a warning or notice.
///
///     0: OK. Each column of the input matrix contained
///       row indices in increasing order, with no
///       duplicates.
///
///     1: OK, but columns of input matrix were jumbled
///       (unsorted columns or duplicate entries). Symamd
///       had to do some extra work to sort the matrix
///       first and remove duplicate entries, but it
///       still was able to return a valid permutation
///       (return value of symamd was true).
///
///         stats[4]: highest numbered column that
///           is unsorted or has duplicate entries.
///         stats[5]: last seen duplicate or unsorted row index.
///         stats[6]: number of duplicate or unsorted row indices.
///     -1: A is a null pointer
///     -2: p is a null pointer
///     -3: (unused, see colamd.go)
///     -4:  n is negative
///         stats[4]: n
///     -5: Number of nonzeros in matrix is negative
///         stats[4]: # of nonzeros (p [n]).
///     -6: p[0] is nonzero
///         stats[4]: p[0]
///     -7: (unused)
///     -8: A column has a negative number of entries
///         stats[4]: column with < 0 entries
///         stats[5]: number of entries in col
///     -9: A row index is out of bounds
///         stats[4]: column with bad row index
///         stats[5]: bad row index
///         stats[6]: n_row, # of rows of matrix
///     -10: Out of memory (unable to allocate temporary
///       workspace for M or count arrays using the
///       "allocate" routine passed into symamd).
pub fn symamd(
    n: i32,
    A: &[i32],
    p: &[i32],
    perm: &mut [i32],
    knobs: Option<[f64; KNOBS]>,
    stats: &mut [i32; STATS],
) -> bool {
    let count: Vec<i32>; // Length of each column of M, and col pointer.
    let mark: Vec<i32>; // Mark array for finding duplicate entries.
                        // let M: Vec<i32>; // Row indices of matrix M.
    let Mlen: i32; // Length of M.
    let nrow: i32; // Number of rows in M.
    let nnz: i32; // Number of entries in A.
                  // let i: i32; // Row index of A.
                  // let j: i32; // Column index of A.
                  // let k: i32; // Row index of M.
    let mnz: i32; // Number of nonzeros in M.
    let pp: i32; // Index into a column of A.
                 // let lastRow: i32; // Last row seen in the current column.
                 // let length: i32; // Number of nonzeros in a column.

    let mut cknobs = [0.0; KNOBS]; // Knobs for colamd.
                                   // let defaultKnobs = [0.0; KNOBS]; // Default knobs for colamd.

    for i in 0..STATS {
        stats[i] = 0;
    }
    stats[STATUS] = OK;
    stats[INFO1] = -1;
    stats[INFO2] = -1;

    if n < 0 {
        // n must be >= 0
        stats[STATUS] = ERROR_NCOL_NEGATIVE;
        stats[INFO1] = n;
        debug0!("symamd: n negative {}", n);
        return false;
    }

    nnz = p[n as usize];
    if nnz < 0 {
        // nnz must be >= 0
        stats[STATUS] = ERROR_NNZ_NEGATIVE;
        stats[INFO1] = nnz;
        debug0!("symamd: number of entries negative %d", nnz);
        return false;
    }

    if p[0] != 0 {
        stats[STATUS] = ERROR_P0_NONZERO;
        stats[INFO1] = p[0];
        debug0!("symamd: p[0] not zero %d", p[0]);
        return false;
    }

    // If no knobs, set default knobs.

    let knobs = if knobs.is_none() {
        default_knobs()
    } else {
        knobs.unwrap()
    };

    // Allocate count and mark.

    let mut count = vec![0; (n + 1) as usize];
    let mut mark = vec![0; (n + 1) as usize];

    // Compute column counts of M, check if A is valid.

    stats[INFO3] = 0; // Number of duplicate or unsorted row indices.

    for i in 0..n as usize {
        mark[i] = -1;
    }

    for j in 0..n {
        let mut lastRow = -1; // Last row seen in the current column.

        let length = p[(j + 1) as usize] - p[j as usize]; // Number of nonzeros in a column.
        if length < 0 {
            // Column pointers must be non-decreasing.
            stats[STATUS] = ERROR_COL_LENGTH_NEGATIVE;
            stats[INFO1] = j;
            stats[INFO2] = length;
            // count = None;
            // mark = None;
            debug0!("symamd: col {} negative length {}", j, length);
            return false;
        }

        for pp in p[j as usize]..p[(j + 1) as usize] {
            let i = A[pp as usize];
            if i < 0 || i >= n {
                // Row index i, in column j, is out of bounds.
                stats[STATUS] = ERROR_ROW_INDEX_OUT_OF_BOUNDS;
                stats[INFO1] = j;
                stats[INFO2] = i;
                stats[INFO3] = n;
                // count = None;
                // mark = None;
                debug0!("symamd: row {} col {} out of bounds", i, j);
                return false;
            }

            if i <= lastRow || mark[i as usize] == j {
                // Row index is unsorted or repeated (or both), thus col
                // is jumbled. This is a notice, not an error condition.
                stats[STATUS] = OK_BUT_JUMBLED;
                stats[INFO1] = j;
                stats[INFO2] = i;
                stats[INFO3] += 1;
                debug1!("symamd: row {} col {} unsorted/duplicate", i, j);
            }

            if i > j && mark[i as usize] != j {
                // Row k of M will contain column indices i and j.
                count[i as usize] += 1;
                count[j as usize] += 1;
            }

            // Mark the row as having been seen in this column.
            mark[i as usize] = j;

            lastRow = i;
        }
    }

    // Compute column pointers of M.

    // Use output permutation, perm, for column pointers of M.
    perm[0] = 0;
    for j in 1..=n as usize {
        perm[j] = perm[j - 1] + count[j - 1];
    }
    for j in 0..n as usize {
        count[j] = perm[j]
    }

    // Construct M.

    mnz = perm[n as usize];
    nrow = mnz / 2;
    Mlen = recommended(mnz, nrow, n);
    let mut M = vec![0; Mlen as usize];
    debug0!(
        "symamd: M is {}-by-{} with {} entries, Mlen = {}",
        nrow,
        n,
        mnz,
        Mlen
    );

    let mut k = 0;

    if stats[STATUS] == OK {
        // Matrix is OK.
        for j in 0..n {
            #[cfg(feature = "debug")]
            {
                assert!(p[j + 1] - p[j] >= 0);
            }
            for pp in p[j as usize]..p[(j + 1) as usize] {
                let i = A[pp as usize];
                #[cfg(feature = "debug")]
                {
                    assert!(i >= 0 && i < n);
                }
                if i > j {
                    // Row k of M contains column indices i and j.
                    M[count[i as usize] as usize] = k;
                    count[i as usize] += 1;
                    M[count[j as usize] as usize] = k;
                    count[j as usize] += 1;
                    k += 1;
                }
            }
        }
    } else {
        // Matrix is jumbled. Do not add duplicates to M. Unsorted cols OK.
        debug0!("symamd: Duplicates in A.");

        for i in 0..n {
            mark[i as usize] = -1;
        }
        for j in 0..n {
            #[cfg(feature = "debug")]
            {
                assert!(p[j + 1] - p[j] >= 0);
            }
            for pp in p[j as usize]..p[(j + 1) as usize] {
                let i = A[pp as usize];
                #[cfg(feature = "debug")]
                {
                    assert!(i >= 0 && i < n);
                }
                if i > j && mark[i as usize] != j {
                    // Row k of M contains column indices i and j.
                    M[count[i as usize] as usize] = k;
                    count[i as usize] += 1;
                    M[count[j as usize] as usize] = k;
                    count[j as usize] += 1;
                    k += 1;
                    mark[i as usize] = j;
                }
            }
        }
    }

    // Count and mark no longer needed.
    // count = nil
    // mark = nil
    #[cfg(feature = "debug")]
    {
        assert!(k == nrow);
    }

    // Adjust the knobs for M.

    for i in 0..KNOBS {
        cknobs[i] = knobs[i];
    }

    // There are no dense rows in M.
    cknobs[DENSE_ROW] = -1.0;
    cknobs[DENSE_COL] = knobs[DENSE_ROW];

    debug0!("symamd: dense col knob for M: {}\n", cknobs[DENSE_COL]);

    // Order the columns of M.

    colamd(nrow, n, Mlen, &mut M, perm, Some(cknobs), stats);

    // Note that the output permutation is now in perm.

    // Get the statistics for symamd from colamd.

    // A dense column in colamd means a dense row and col in symamd.
    stats[DENSE_ROW] = stats[DENSE_COL];

    debug0!("symamd: done.");

    false
}

use crate::codes::*;
use crate::col::Col;
use crate::internal::*;
use crate::row::Row;
use crate::stats::*;

// Takes the column form of the matrix in A and creates the row form of the
// matrix. Also, row and column attributes are stored in the Col and Row
// structs. If the columns are un-sorted or contain duplicate row indices,
// this routine will also sort and remove duplicate row indices from the
// column form of the matrix. Returns false if the matrix is invalid,
// true otherwise.
pub(crate) fn init_rows_cols(
    nrow: i32,
    ncol: i32,
    rows: &mut [Row],
    cols: &mut [Col],
    A: &[i32],
    p: &mut [i32],
    stats: &mut [i32; STATS],
) -> bool {
    // Initialize columns, and check column pointers.
    for col in 0..ncol as usize {
        cols[col].start = p[col];
        cols[col].length = p[col + 1] - p[col];

        if cols[col].length < 0 {
            // Column pointers must be non-decreasing.
            stats[STATUS] = ERROR_COL_LENGTH_NEGATIVE;
            stats[INFO1] = col as i32;
            stats[INFO2] = cols[col].length;
            debug0!("colamd: col {} length {} < 0", col, cols[col].length);
            return false;
        }

        cols[col].set_thickness(1);
        cols[col].set_score(0);
        cols[col].set_prev(empty);
        cols[col].set_degree_next(empty);
    }
    // p[0..n_col] no longer needed, used as "head" in subsequent routines.

    // Scan columns, compute row degrees, and check row indices.

    // Number of duplicate or unsorted row indices.
    stats[INFO3] = 0;

    for row in 0..nrow as usize {
        rows[row].length = 0;
        rows[row].set_mark(-1);
    }

    for col in 0..ncol as usize {
        let mut lastRow = -1;

        let mut cp = p[col] as usize;
        let cpend = p[col + 1] as usize;

        while cp < cpend {
            let row = A[cp];
            cp += 1;

            // Make sure row indices within range.
            if row < 0 || row >= nrow {
                stats[STATUS] = ERROR_ROW_INDEX_OUT_OF_BOUNDS;
                stats[INFO1] = col as i32;
                stats[INFO2] = row;
                stats[INFO3] = nrow;
                debug0!("colamd: row {} col {} out of bounds", row, col);
                return false;
            }

            if row <= lastRow || rows[row].mark() == col {
                // Row index are unsorted or repeated (or both), thus col
                // is jumbled. This is a notice, not an error condition.
                stats[STATUS] = OK_BUT_JUMBLED;
                stats[INFO1] = col as i32;
                stats[INFO2] = row;
                stats[INFO3] += 1;
                debug1!("colamd: row {} col {} unsorted/duplicate", row, col);
            }

            if rows[row].mark() != col {
                rows[row].length += 1;
            } else {
                // This is a repeated entry in the column, it will be removed.
                cols[col].length -= 1;
            }

            // Mark the row as having been seen in this column.
            rows[row].setMark(col);

            lastRow = row;
        }
    }

    // Compute row pointers.

    // Row form of the matrix starts directly after the column
    // form of matrix in A.
    rows[0].start = p[ncol];
    rows[0].setP(rows[0].start);
    rows[0].setMark(-1);
    for row in 1..nrow {
        rows[row].start = rows[row - 1].start + rows[row - 1].length;
        rows[row].setP(rows[row].start);
        rows[row].setMark(-1);
    }

    // Create row form.

    if stats[STATUS] == OK_BUT_JUMBLED {
        // If cols jumbled, watch for repeated row indices.
        for col in 0..ncol {
            let cp = p[col];
            let cpend = p[col + 1];

            while cp < cpend {
                let row = A[cp];
                cp += 1;

                if rows[row].mark() != col {
                    A[rows[row].p()] = col;
                    rows[row].setP(rows[row].p() + 1);
                    rows[row].setMark(col);
                }
            }
        }
    } else {
        // If cols not jumbled, we don't need the mark (this is faster).
        for col in 0..ncol {
            let cp = p[col];
            let cpend = p[col + 1];
            while cp < cpend {
                A[rows[A[cp]].p()] = col;
                rows[A[cp]].setP(rows[A[cp]].p() + 1);
                cp += 1;
            }
        }
    }

    // Clear the row marks and set row degrees.

    for row in 0..nrow {
        rows[row].setMark(0);
        rows[row].setDegree(rows[row].length);
    }

    // See if we need to re-create columns.

    if stats[STATUS] == OK_BUT_JUMBLED {
        #[cfg(feature = "debug")]
        {
            debug0!("colamd: reconstructing column form, matrix jumbled");

            // Make sure column lengths are correct.
            for col in 0..ncol {
                p[col] = cols[col].length;
            }
            for row in 0..nrow {
                let rp = rows[row].start;
                let rpend = rp + rows[row].length;
                while rp < rpend {
                    p[A[rp]] -= 1;
                    rp += 1;
                }
            }
            for col in 0..ncol {
                assert!(p[col] == 0)
            }
            // Now p is all zero (different than when debugging is turned off).
        }
        // ndebug

        // Compute col pointers.

        // Col form of the matrix starts at A [0].
        // Note, we may have a gap between the col form and the row
        // form if there were duplicate entries, if so, it will be
        // removed upon the first garbage collection.
        cols[0].start = 0;
        p[0] = cols[0].start;
        for col in 1..ncol {
            // Note that the lengths here are for pruned columns, i.e.
            // no duplicate row indices will exist for these columns.
            cols[col].start = cols[col - 1].start + cols[col - 1].length;
            p[col] = cols[col].start;
        }

        // Re-create col form.

        for row in 0..nrow {
            let rp = rows[row].start;
            let rpend = rp + rows[row].length;
            while rp < rpend {
                A[p[A[rp]]] = row;
                p[A[rp]] += 1;
                rp += 1;
            }
        }
    }

    // Done. Matrix is not (or no longer) jumbled.

    return true;
}

// Kills dense or empty columns and rows, calculates an initial score for
// each column, and places all columns in the degree lists.
pub(crate) fn init_scoring(
    nrow: i32,
    ncol: i32,
    rows: &[Row],
    cols: &[Col],
    A: &[i32],
    head: &[i32],
    knobs: [f64; KNOBS],
    pnrow2: &mut i32,
    pncol2: &mut i32,
    pmaxDeg: &mut i32,
) {
    let c: i32; // A column index.
    let r: i32; // A row index.
    let row: i32; // A row index.
    let cp: i32; // A column pointer.
    let deg: i32; // Degree of a row or column.
    let cpend: i32; // A pointer to the end of a column.
    let newcp: i32; // New column pointer.
    let colLength: i32; // Length of pruned column.
    let score: i32; // Current column score.
    let ncol2: i32; // Number of non-dense, non-empty columns.
    let nrow2: i32; // Number of non-dense, non-empty rows.
    let denseRowCount: i32; // Remove rows with more entries than this.
    let denseColCount: i32; // Remove cols with more entries than this.
    let minScore: i32; // Smallest column score.
    let maxDeg: i32; // Maximum row degree.
    let nextCol: i32; // Used to add to degree list.

    #[cfg(feature = "debug")]
    let mut debugCount = 0; // Debug only.

    // Extract knobs.

    // Note: if knobs contains a NaN, this is undefined:
    let denseRowCount = if knobs[DENSE_ROW] < 0.0 {
        // Only remove completely dense rows.
        ncol - 1
    } else {
        dense_degree(knobs[DENSE_ROW], ncol)
    };
    let denseColCount = if knobs[DENSE_COL] < 0.0 {
        // Only remove completely dense columns.
        nrow - 1
    } else {
        dense_degree(knobs[DENSE_COL], i32::min(nrow, ncol))
    };

    debug1!("colamd: densecount: {} {}", denseRowCount, denseColCount);
    maxDeg = 0;
    ncol2 = ncol;
    nrow2 = nrow;

    // Kill empty columns.

    // Put the empty columns at the end in their natural order, so that LU
    // factorization can proceed as far as possible.
    // for c = ncol - 1; c >= 0; c-- {  TODO: check
    for c in ((ncol - 1)..=0).rev() {
        deg = cols[c].length;
        if deg == 0 {
            // This is a empty column, kill and order it last.
            ncol2 -= 1;
            cols[c].setOrder(ncol2);
            kill_principal_col(cols, c);
        }
    }
    debug1!("colamd: null columns killed: {}", ncol - ncol2);

    // Kill dense columns.

    // Put the dense columns at the end, in their natural order.
    // for c = ncol - 1; c >= 0; c-- {  TODO: check
    for c in (ncol - 1..=0).rev() {
        // Skip any dead columns.
        if col_is_dead(cols, c) {
            continue;
        }
        deg = cols[c].length;
        if deg > denseColCount {
            // This is a dense column, kill and order it last.
            ncol2 -= 1;
            cols[c].setOrder(ncol2);
            // Decrement the row degrees.
            cp = cols[c].start;
            cpend = cp + cols[c].length;
            while cp < cpend {
                rows[A[cp]].setDegree(rows[A[cp]].degree() - 1);
                cp += 1;
            }
            kill_principal_col(cols, c);
        }
    }
    debug1!("colamd: Dense and null columns killed: {}", ncol - ncol2);

    // Kill dense and empty rows.
    for r in 0..nrow {
        deg = rows[r].degree();
        debug_assert!(deg >= 0 && deg <= ncol);
        if deg > denseRowCount || deg == 0 {
            // Kill a dense or empty row.
            kill_row(rows, r);
            nrow2 -= 1;
        } else {
            // keep track of max degree of remaining rows
            maxDeg = i32::max(maxDeg, deg);
        }
    }
    debug1!("colamd: Dense and null rows killed: {}", nrow - nrow2);
}

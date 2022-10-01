use crate::col::{Col, ColShared1, ColShared2, ColShared3, ColShared4};
use crate::internal::*;
use crate::row::{Row, RowShared1, RowShared2};
use crate::stats::*;

#[cfg(feature = "debug")]
use crate::debug::*;

// Takes the column form of the matrix in A and creates the row form of the
// matrix. Also, row and column attributes are stored in the Col and Row
// structs. If the columns are un-sorted or contain duplicate row indices,
// this routine will also sort and remove duplicate row indices from the
// column form of the matrix. Returns false if the matrix is invalid,
// true otherwise.
pub(crate) fn init_rows_cols(
    nrow: Int,
    ncol: Int,
    rows: &mut [Row],
    cols: &mut [Col],
    a_i: &mut [Int],
    p: &mut [Int],
    stats: &mut [Int; STATS],
) -> bool {
    // Initialize columns, and check column pointers.
    for col in 0..ncol as usize {
        cols[col].start = p[col];
        cols[col].length = p[col + 1] - p[col];

        if cols[col].length < 0 {
            // Column pointers must be non-decreasing.
            stats[STATUS] = ERROR_COL_LENGTH_NEGATIVE;
            stats[INFO1] = col as Int;
            stats[INFO2] = cols[col].length;
            debug0!("colamd: col {} length {} < 0", col, cols[col].length);
            return false;
        }

        cols[col].shared1 = ColShared1::Thickness(1);
        cols[col].shared2 = ColShared2::Score(0);
        cols[col].shared3 = ColShared3::Prev(EMPTY);
        cols[col].shared4 = ColShared4::DegreeNext(EMPTY);
    }
    // p[0..n_col] no longer needed, used as "head" in subsequent routines.

    // Scan columns, compute row degrees, and check row indices.

    // Number of duplicate or unsorted row indices.
    stats[INFO3] = 0;

    for row in 0..nrow as usize {
        rows[row].length = 0;
        rows[row].shared2 = RowShared2::Mark(-1);
    }

    for col in 0..ncol {
        let mut last_row = -1;

        let mut cp = p[col as usize] as usize;
        let cp_end = p[(col + 1) as usize] as usize;

        while cp < cp_end {
            let row = a_i[cp];
            cp += 1;

            // Make sure row indices within range.
            if row < 0 || row >= nrow {
                stats[STATUS] = ERROR_ROW_INDEX_OUT_OF_BOUNDS;
                stats[INFO1] = col;
                stats[INFO2] = row;
                stats[INFO3] = nrow;
                debug0!("colamd: row {} col {} out of bounds", row, col);
                return false;
            }

            if row <= last_row || rows[row as usize].shared2.mark() == col {
                // Row index are unsorted or repeated (or both), thus col
                // is jumbled. This is a notice, not an error condition.
                stats[STATUS] = OK_BUT_JUMBLED;
                stats[INFO1] = col;
                stats[INFO2] = row;
                stats[INFO3] += 1;
                debug1!("colamd: row {} col {} unsorted/duplicate", row, col);
            }

            if rows[row as usize].shared2.mark() != col {
                rows[row as usize].length += 1;
            } else {
                // This is a repeated entry in the column, it will be removed.
                cols[col as usize].length -= 1;
            }

            // Mark the row as having been seen in this column.
            rows[row as usize].shared2 = RowShared2::Mark(col);

            last_row = row;
        }
    }

    // Compute row pointers.

    // Row form of the matrix starts directly after the column
    // form of matrix in A.
    rows[0].start = p[ncol as usize];
    rows[0].shared1 = RowShared1::P(rows[0].start);
    rows[0].shared2 = RowShared2::Mark(-1);
    for row in 1..nrow as usize {
        rows[row].start = rows[row - 1].start + rows[row - 1].length;
        rows[row].shared1 = RowShared1::P(rows[row].start);
        rows[row].shared2 = RowShared2::Mark(-1);
    }

    // Create row form.

    if stats[STATUS] == OK_BUT_JUMBLED {
        // If cols jumbled, watch for repeated row indices.
        for col in 0..ncol {
            let mut cp = p[col as usize];
            let cp_end = p[(col + 1) as usize];

            while cp < cp_end {
                let row = a_i[cp as usize] as usize;
                cp += 1;

                if rows[row].shared2.mark() != col {
                    a_i[rows[row].shared1.p() as usize] = col;
                    rows[row].shared1 = RowShared1::P(rows[row].shared1.p() + 1);
                    rows[row].shared2 = RowShared2::Mark(col);
                }
            }
        }
    } else {
        // If cols not jumbled, we don't need the mark (this is faster).
        for col in 0..ncol {
            let mut cp = p[col as usize] as usize;
            let cp_end = p[(col + 1) as usize] as usize;
            while cp < cp_end {
                a_i[rows[a_i[cp] as usize].shared1.p() as usize] = col;
                rows[a_i[cp] as usize].shared1 =
                    RowShared1::P(rows[a_i[cp] as usize].shared1.p() + 1);
                cp += 1;
            }
        }
    }

    // Clear the row marks and set row degrees.

    for row in 0..nrow as usize {
        rows[row].shared2 = RowShared2::Mark(0);
        rows[row].shared1 = RowShared1::Degree(rows[row].length);
    }

    // See if we need to re-create columns.

    if stats[STATUS] == OK_BUT_JUMBLED {
        #[cfg(feature = "debug")]
        {
            debug0!("colamd: reconstructing column form, matrix jumbled");

            // Make sure column lengths are correct.
            for col in 0..ncol as usize {
                p[col] = cols[col].length;
            }
            for row in 0..nrow as usize {
                let mut rp = rows[row].start as usize;
                let rp_end = rp + rows[row].length as usize;
                while rp < rp_end {
                    p[a_i[rp] as usize] -= 1;
                    rp += 1;
                }
            }
            for col in 0..ncol as usize {
                assert_debug!(p[col] == 0)
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
        for col in 1..ncol as usize {
            // Note that the lengths here are for pruned columns, i.e.
            // no duplicate row indices will exist for these columns.
            cols[col].start = cols[col - 1].start + cols[col - 1].length;
            p[col] = cols[col].start;
        }

        // Re-create col form.

        for row in 0..nrow {
            let mut rp = rows[row as usize].start as usize;
            let rp_end = rp + rows[row as usize].length as usize;
            while rp < rp_end {
                a_i[p[a_i[rp] as usize] as usize] = row;
                p[a_i[rp] as usize] += 1;
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
    nrow: Int,
    ncol: Int,
    rows: &mut [Row],
    cols: &mut [Col],
    a_i: &mut [Int],
    head: &mut [Int],
    knobs: [f64; KNOBS],
    pnrow2: &mut Int,
    pncol2: &mut Int,
    pmax_deg: &mut Int,
) {
    // Extract knobs.

    // Remove rows with more entries than this.
    // Note: if knobs contains a NaN, this is undefined:
    let dense_row_count = if knobs[DENSE_ROW] < 0.0 {
        // Only remove completely dense rows.
        ncol - 1
    } else {
        dense_degree(knobs[DENSE_ROW], ncol)
    };
    // Remove cols with more entries than this.
    let dense_col_count = if knobs[DENSE_COL] < 0.0 {
        // Only remove completely dense columns.
        nrow - 1
    } else {
        dense_degree(knobs[DENSE_COL], Int::min(nrow, ncol))
    };

    debug1!(
        "colamd: densecount: {} {}",
        dense_row_count,
        dense_col_count
    );
    let mut max_deg = 0; // Maximum row degree.
    let mut ncol2 = ncol; // Number of non-dense, non-empty columns.
    let mut nrow2 = nrow; // Number of non-dense, non-empty rows.

    // Kill empty columns.

    // Put the empty columns at the end in their natural order, so that LU
    // factorization can proceed as far as possible.
    // for c = ncol - 1; c >= 0; c-- {  TODO: check
    for c in (0..ncol as usize).rev() {
        let deg = cols[c].length; // Degree of column.
        if deg == 0 {
            // This is a empty column, kill and order it last.
            ncol2 -= 1;
            cols[c].shared2 = ColShared2::Order(ncol2);
            kill_principal_col(cols, c);
        }
    }
    debug1!("colamd: null columns killed: {}", ncol - ncol2);

    // Kill dense columns.

    // Put the dense columns at the end, in their natural order.
    // for c = ncol - 1; c >= 0; c-- {  TODO: check
    for c in (0..ncol as usize).rev() {
        // Skip any dead columns.
        if col_is_dead(cols, c) {
            continue;
        }
        let deg = cols[c].length; // Degree of column.
        if deg > dense_col_count {
            // This is a dense column, kill and order it last.
            ncol2 -= 1;
            cols[c].shared2 = ColShared2::Order(ncol2);
            // Decrement the row degrees.
            let mut cp = cols[c].start as usize; // Column pointer.
            let cp_end = cp + cols[c].length as usize; // Pointer to the end of the column.
            while cp < cp_end {
                rows[a_i[cp] as usize].shared1 =
                    RowShared1::Degree(rows[a_i[cp] as usize].shared1.degree() - 1);
                cp += 1;
            }
            kill_principal_col(cols, c);
        }
    }
    debug1!("colamd: Dense and null columns killed: {}", ncol - ncol2);

    // Kill dense and empty rows.
    for r in 0..nrow as usize {
        let deg = rows[r].shared1.degree(); // Degree of row.
        assert_debug!(deg >= 0 && deg <= ncol);
        if deg > dense_row_count || deg == 0 {
            // Kill a dense or empty row.
            kill_row(rows, r);
            nrow2 -= 1;
        } else {
            // keep track of max degree of remaining rows
            max_deg = Int::max(max_deg, deg);
        }
    }
    debug1!("colamd: Dense and null rows killed: {}", nrow - nrow2);

    // Compute initial column scores.

    // At this point the row degrees are accurate. They reflect the number
    // of "live" (non-dense) columns in each row. No empty rows exist.
    // Some "live" columns may contain only dead rows, however. These are
    // pruned in the code below.

    // Now find the initial matlab score for each column.
    // for c = ncol - 1; c >= 0; c-- { TODO: check
    for c in (0..ncol as usize).rev() {
        // Skip dead column.
        if col_is_dead(cols, c) {
            continue;
        }
        let mut score = 0; // Current column score.
        let mut cp = cols[c].start as usize;
        let mut new_cp = cp; // New column pointer.
        let cp_end = cp + cols[c].length as usize;
        while cp < cp_end {
            // Get a row.
            let row = a_i[cp] as usize; // row index
            cp += 1;
            // Skip if dead.
            if row_is_dead(rows, row) {
                continue;
            }
            // Compact the column.
            a_i[new_cp] = row as Int;
            new_cp += 1;
            // Add row's external degree.
            score += rows[row].shared1.degree() - 1;
            // Guard against integer overflow.
            score = Int::min(score, ncol);
        }
        // Determine pruned column length.
        let col_length = new_cp as Int - cols[c].start; // Length of pruned column.
        if col_length == 0 {
            // A newly-made null column (all rows in this col are "dense"
            // and have already been killed).
            debug2!("Newly null killed: {}", c);
            ncol2 -= 1;
            cols[c].shared2 = ColShared2::Order(ncol2);
            kill_principal_col(cols, c);
        } else {
            // Set column length and set score.
            assert_debug!(score >= 0);
            assert_debug!(score <= ncol);
            cols[c].length = col_length;
            cols[c].shared2 = ColShared2::Score(score);
        }
    }
    debug1!(
        "colamd: Dense, null, and newly-null columns killed: {}",
        ncol - ncol2
    );

    // At this point, all empty rows and columns are dead. All live columns
    // are "clean" (containing no dead rows) and simplicial (no supercolumns
    // yet). Rows may contain dead columns, but all live rows contain at
    // least one live column.

    #[cfg(feature = "debug")]
    debug_structures(nrow, ncol, rows, cols, a_i, ncol2);

    // Initialize degree lists.

    #[cfg(feature = "debug")]
    let mut debug_count = 0;

    // Clear the hash buckets.
    for c in 0..=ncol as usize {
        head[c] = EMPTY;
    }
    let mut min_score = ncol; // Smallest column score.

    // Place in reverse order, so low column indices are at the front
    // of the lists. This is to encourage natural tie-breaking.
    // for c = ncol - 1; c >= 0; c-- {
    for c in (0..ncol as usize).rev() {
        // Only add principal columns to degree lists.
        if col_is_alive(cols, c) {
            debug4!(
                "place {} score {} minscore {} ncol {}",
                c,
                cols[c].shared2.score(),
                min_score,
                ncol
            );

            // Add columns score to DList.

            let score = cols[c].shared2.score();

            assert_debug!(min_score >= 0);
            assert_debug!(min_score <= ncol);
            assert_debug!(score >= 0);
            assert_debug!(score <= ncol);
            assert_debug!(head[score as usize] >= EMPTY);

            // Now add this column to dList at proper score location.
            let next_col = head[score as usize]; // Used to add to degree list.
            cols[c].shared3 = ColShared3::Prev(EMPTY);
            cols[c].shared4 = ColShared4::DegreeNext(next_col);

            // if there already was a column with the same score, set its
            // previous pointer to this new column.
            if next_col != EMPTY {
                cols[next_col as usize].shared3 = ColShared3::Prev(c as Int);
            }
            head[score as usize] = c as Int;

            // See if this score is less than current min.
            min_score = Int::min(min_score, score);

            #[cfg(feature = "debug")]
            {
                debug_count += 1;
            }
        }
    }

    debug1!(
        "colamd: Live cols {} out of {}, non-princ: {}",
        debug_count,
        ncol,
        ncol - debug_count
    );

    assert_debug!(debug_count == ncol2);
    #[cfg(feature = "debug")]
    debug_deg_lists(nrow, ncol, rows, cols, head, min_score, ncol2, max_deg);

    // Return number of remaining columns, and max row degree.

    *pncol2 = ncol2;
    *pnrow2 = nrow2;
    *pmax_deg = max_deg;
}

// Order the principal columns of the supercolumn form of the matrix
// (no supercolumns on input). Uses a minimum approximate column minimum
// degree ordering method.
//
// `a_len` is the size of `A` and must be `2*nnz + n_col` or larger.
// `n_col2` is the number of remaining columns to order.
// `max_deg` is the maximum row degree.
// `pfree` is the index of first free slot (`2*nnz` on entry).
// Returns the number of garbage collections.
pub(crate) fn find_ordering(
    nrow: Int,
    ncol: Int,
    a_len: Int,
    rows: &mut [Row],
    cols: &mut [Col],
    a_i: &mut [Int],
    head: &mut [Int],
    n_col2: Int,
    mut maxdeg: Int,
    mut pfree: Int,
    aggressive: bool,
) -> Int {
    #[cfg(feature = "debug")]
    let mut debug_step = 0; // Debug loop counter.

    // Initialization and clear mark.

    let max_mark = Int::MAX - ncol; // Maximum value of tag_mark.
    let mut tag_mark = clear_mark(0, max_mark, nrow, rows); // Marker value for mark array.
    let mut min_score = 0; // Smallest column score.
    let mut ngarbage = 0; // Number of garbage collections performed.
    debug1!("colamd: Ordering, n_col2={}", n_col2);

    // Order the columns.

    let mut k = 0; // 'k' is incremented below
    while k < n_col2 {
        #[cfg(feature = "debug")]
        {
            if debug_step % 100 == 0 {
                debug2!("\n...      Step k: {} out of n_col2: {}", k, n_col2);
            } else {
                debug3!("\n----------Step k: {} out of n_col2: {}", k, n_col2);
            }
            debug_step += 1;
            debug_deg_lists(nrow, ncol, rows, cols, head, min_score, n_col2 - k, maxdeg);
            debug_matrix(nrow, ncol, rows, cols, a_i);
        }
        // ndebug

        // Select pivot column, and order it.

        #[cfg(feature = "debug")]
        {
            // Make sure degree list isn't empty.
            assert_debug!(min_score >= 0);
            assert_debug!(min_score <= ncol);
            assert_debug!(head[min_score as usize] >= EMPTY);

            for debugd in 0..min_score as usize {
                assert_debug!(head[debugd] == EMPTY);
            }
        }

        // Get pivot column from head of minimum degree list.
        while head[min_score as usize] == EMPTY && min_score < ncol {
            min_score += 1;
        }
        let pivot_col = head[min_score as usize] as usize; // Current pivot column.
        assert_debug!(/*pivot_col >= 0 &&*/ pivot_col <= ncol as usize);

        let next_col = cols[pivot_col].shared4.degree_next();
        head[min_score as usize] = next_col;
        if next_col != EMPTY {
            cols[next_col as usize].shared3 = ColShared3::Prev(EMPTY);
        }

        assert_debug!(col_is_alive(cols, pivot_col));
        debug3!("Pivot col: {}\n", pivot_col);

        // Remember score for defrag check.
        let pivot_col_score = cols[pivot_col].shared2.score(); // Score of pivot column.

        // The pivot column is the kth column in the pivot order.
        cols[pivot_col].shared2 = ColShared2::Order(k);

        // Increment order count by column thickness.
        let pivot_col_thickness = cols[pivot_col].shared1.thickness(); // Number of columns represented by pivot col.
        k += pivot_col_thickness;

        assert_debug!(pivot_col_thickness > 0);
        debug3!("Pivot col: {} thick {}", pivot_col, pivot_col_thickness);

        // Garbage_collection, if necessary.

        let needed_memory = Int::min(pivot_col_score, ncol - k); // Free space needed for pivot row.
        if pfree + needed_memory >= a_len {
            pfree = garbage_collection(nrow, ncol, rows, cols, a_i, pfree);
            ngarbage += 1;
            // After garbage collection we will have enough.
            assert_debug!(pfree + needed_memory < a_len);

            // Garbage collection has wiped out the Row[].shared2.mark array.
            tag_mark = clear_mark(0, max_mark, nrow, rows);

            #[cfg(feature = "debug")]
            debug_matrix(nrow, ncol, rows, cols, a_i);
        }

        // Compute pivot row pattern.

        // Get starting location for this new merged row.
        let pivot_row_start = pfree; // Pointer to start of pivot row.

        // Initialize new row counts to zero.
        let mut pivot_row_degree = 0; // Number of columns in pivot row.

        // Tag pivot column as having been visited so it isn't included
        // in merged pivot row.
        cols[pivot_col].shared1 = ColShared1::Thickness(-pivot_col_thickness);

        // Pivot row is the union of all rows in the pivot column pattern.
        let mut cp = cols[pivot_col].start as usize;
        let cp_end = cp + cols[pivot_col].length as usize; // Pointer to the end of the column.
        while cp < cp_end {
            // Get a row.
            let row = a_i[cp] as usize;
            cp += 1;
            debug4!("Pivot col pattern {} {}", row_is_alive(rows, row), row);

            // Skip if row is dead.
            if row_is_alive(rows, row) {
                let mut rp = rows[row].start as usize;
                let rp_end = rp + rows[row].length as usize; // Pointer to the end of the row.
                while rp < rp_end {
                    // Get a column.
                    let col = a_i[rp] as usize;
                    rp += 1;
                    // Add the column, if alive and untagged.
                    let col_thickness = cols[col].shared1.thickness(); // "thickness" (no. of columns in a supercol)
                    if col_thickness > 0 && col_is_alive(cols, col) {
                        // Tag column in pivot row.
                        cols[col].shared1 = ColShared1::Thickness(-col_thickness);
                        assert_debug!(pfree < a_len);

                        // Place column in pivot row.
                        a_i[pfree as usize] = col as Int;
                        pfree += 1;
                        pivot_row_degree += col_thickness;
                    }
                }
            }
        }

        // Clear tag on pivot column.
        cols[pivot_col].shared1 = ColShared1::Thickness(pivot_col_thickness);
        maxdeg = Int::max(maxdeg, pivot_row_degree);

        #[cfg(feature = "debug")]
        {
            debug3!("check2");
            debug_mark(nrow, rows, tag_mark, max_mark);
        }

        // Kill all rows used to construct pivot row.

        // Also kill pivot row, temporarily.
        let mut cp = cols[pivot_col].start as usize;
        let cp_end = cp + cols[pivot_col].length as usize; // Pointer to the end of the column.
        while cp < cp_end {
            // May be killing an already dead row.
            let row = a_i[cp] as usize;
            cp += 1;
            debug3!("Kill row in pivot col: {}", row);

            kill_row(rows, row)
        }

        // Select a row index to use as the new pivot row.

        let pivot_row_length = pfree - pivot_row_start; // Number of supercolumns in pivot row.
        let pivot_row = if pivot_row_length > 0 {
            // Pick the "pivot" row arbitrarily (first row in col).
            let pivot_row = a_i[cols[pivot_col].start as usize];
            debug3!("Pivotal row is {}", pivot_row);
            pivot_row
        } else {
            // There is no pivot row, since it is of zero length.
            let pivot_row = EMPTY;
            assert_debug!(pivot_row_length == 0);
            pivot_row
        };
        assert_debug!(cols[pivot_col].length > 0 || pivot_row_length == 0);

        // Approximate degree computation.

        // Here begins the computation of the approximate degree. The column
        // score is the sum of the pivot row "length", plus the size of the
        // set differences of each row in the column minus the pattern of the
        // pivot row itself. The column ("thickness") itself is also
        // excluded from the column score (we thus use an approximate
        // external degree).

        // The time taken by the following code (compute set differences, and
        // add them up) is proportional to the size of the data structure
        // being scanned - that is, the sum of the sizes of each column in
        // the pivot row. Thus, the amortized time to compute a column score
        // is proportional to the size of that column (where size, in this
        // context, is the column "length", or the number of row indices
        // in that column). The number of row indices in a column is
        // monotonically non-decreasing, from the length of the original
        // column on input to colamd.

        // Compute set differences.

        debug3!("** Computing set differences phase. **");

        // Pivot row is currently dead - it will be revived later.

        debug3!("Pivot row: ");

        // For each column in pivot row.
        let mut rp = pivot_row_start as usize;
        let rp_end = rp + pivot_row_length as usize;
        while rp < rp_end {
            let col = a_i[rp] as usize;
            rp += 1;

            assert_debug!(col_is_alive(cols, col) && col != pivot_col);
            debug3!("Col: {}", col);

            // Clear tags used to construct pivot row pattern.
            let col_thickness = -cols[col].shared1.thickness();
            assert_debug!(col_thickness > 0);

            cols[col].shared1 = ColShared1::Thickness(col_thickness);

            // Remove column from degree list.

            let cur_score = cols[col].shared2.score(); // Score of current column.
            let prev_col = cols[col].shared3.unwrap(); //.prev();
            let next_col = cols[col].shared4.degree_next();

            assert_debug!(cur_score >= 0);
            assert_debug!(cur_score <= ncol);
            assert_debug!(cur_score >= EMPTY);

            if prev_col == EMPTY {
                head[cur_score as usize] = next_col;
            } else {
                cols[prev_col as usize].shared4 = ColShared4::DegreeNext(next_col);
            }
            if next_col != EMPTY {
                cols[next_col as usize].shared3 = ColShared3::Prev(prev_col);
            }

            // Scan the column.

            let mut cp = cols[col].start as usize;
            let cp_end = cp + cols[col].length as usize; // Pointer to the end of a column.
            while cp < cp_end {
                // Get a row.
                let row = a_i[cp] as usize;
                cp += 1;
                let row_mark = rows[row].shared2.mark();
                // Skip if dead.
                if row_is_marked_dead(row_mark) {
                    continue;
                }
                assert_debug!(row != pivot_row as usize);

                // Set difference size of row with pivot row.
                let mut set_difference = row_mark - tag_mark;
                // Check if the row has been seen yet.
                if set_difference < 0 {
                    assert_debug!(rows[row].shared1.degree() <= maxdeg);

                    set_difference = rows[row].shared1.degree();
                }
                // Subtract column thickness from this row's set difference.
                set_difference -= col_thickness;
                assert_debug!(set_difference >= 0);

                // Absorb this row if the set difference becomes zero.
                if set_difference == 0 && aggressive {
                    debug3!("aggressive absorption. Row: {}", row);

                    kill_row(rows, row);
                } else {
                    // Save the new mark.
                    rows[row].shared2 = RowShared2::Mark(set_difference + tag_mark);
                }
            }
        }

        #[cfg(feature = "debug")]
        debug_deg_lists(
            nrow,
            ncol,
            rows,
            cols,
            head,
            min_score,
            n_col2 - k - pivot_row_degree,
            maxdeg,
        );

        // Add up set differences for each column.

        debug3!("** Adding set differences phase. **");

        // For each column in pivot row.
        let mut rp = pivot_row_start as usize;
        let rp_end = rp + pivot_row_length as usize;
        while rp < rp_end {
            // Get a column.
            let col = a_i[rp] as usize;
            rp += 1;
            assert_debug!(col_is_alive(cols, col) && col != pivot_col);

            let mut hash: usize = 0; // Hash value for supernode detection.
            let mut cur_score = 0; // Score of current column.
            let mut cp = cols[col].start as usize;
            // Compact the column.
            let mut new_cp = cp as Int; // Modified column pointer.
            let cp_end = cp + cols[col].length as usize;

            debug4!("Adding set diffs for Col: {}.", col);

            while cp < cp_end {
                // Get a row.
                let row = a_i[cp] as usize;
                cp += 1;
                assert_debug!(/*row >= 0 &&*/ row < nrow as usize);

                let row_mark = rows[row].shared2.mark();
                // Skip if dead.
                if row_is_marked_dead(row_mark) {
                    debug4!(" Row {}, dead", row);
                    continue;
                }
                debug4!(" Row {}, set diff {}", row, row_mark - tag_mark);
                assert_debug!(row_mark >= tag_mark);

                // Compact the column.
                a_i[new_cp as usize] = row as Int;
                new_cp += 1;
                // Compute hash function.
                hash += row as usize;
                // Add set difference.
                cur_score += row_mark - tag_mark;
                // Integer overflow...
                cur_score = Int::min(cur_score, ncol);
            }

            // Recompute the column's length.
            cols[col].length = new_cp - cols[col].start;

            // Further mass elimination.

            if cols[col].length == 0 {
                debug4!("further mass elimination. Col: {}", col);

                // Nothing left but the pivot row in this column.
                kill_principal_col(cols, col);
                pivot_row_degree -= cols[col].shared1.thickness();
                assert_debug!(pivot_row_degree >= 0);

                // Order it.
                cols[col].shared2 = ColShared2::Order(k);
                // Increment order count by column thickness.
                k += cols[col].shared1.thickness();
            } else {
                // Prepare for supercolumn detection.

                debug4!("Preparing supercol detection for Col: {}.", col);

                // Save score so far.
                cols[col].shared2 = ColShared2::Score(cur_score);

                // Add column to hash table, for supercolumn detection.
                hash %= (ncol as usize) + 1;

                debug4!(" Hash = {}, n_col = {}.", hash, ncol);
                assert_debug!((hash as Int) <= ncol);

                let head_column = head[hash]; // Head of hash bucket.
                                              // First column in hash bucket.
                let first_col = if head_column > EMPTY {
                    // Degree list "hash" is non-empty, use prev (shared3) of
                    // first column in degree list as head of hash bucket.
                    let first_col = cols[head_column as usize].shared3.prev(); //.head_hash();
                    cols[head_column as usize].shared3 = ColShared3::HeadHash(col as Int);
                    first_col
                } else {
                    // degree list "hash" is empty, use head as hash bucket.
                    let first_col = -(head_column + 2);
                    head[hash] = -1 * (col + 2) as Int;
                    first_col
                };
                cols[col].shared4 = ColShared4::HashNext(first_col);

                // Save hash function in Col [col].shared3.hash
                cols[col].shared3 = ColShared3::Hash(hash as Int);
                assert_debug!(col_is_alive(cols, col));
            }
        }

        // The approximate external column degree is now computed.

        // Supercolumn detection.

        debug3!("** Supercolumn detection phase. **");

        #[cfg(feature = "debug")]
        detect_super_cols(
            ncol,
            Some(rows),
            cols,
            a_i,
            head,
            pivot_row_start,
            pivot_row_length,
        );
        #[cfg(not(feature = "debug"))]
        detect_super_cols(0, None, cols, a_i, head, pivot_row_start, pivot_row_length);

        // Kill the pivotal column.

        kill_principal_col(cols, pivot_col);

        // Clear mark.

        tag_mark = clear_mark(tag_mark + maxdeg + 1, max_mark, nrow, rows);

        #[cfg(feature = "debug")]
        {
            debug3!("check3");
            debug_mark(nrow, rows, tag_mark, max_mark);
        }

        // Finalize the new pivot row, and column scores.

        debug3!("** Finalize scores phase. **");

        // For each column in pivot row.
        let mut rp = pivot_row_start as usize;
        // Compact the pivot row.
        let mut newrp = rp as Int; // Modified row pointer.
        let rp_end = rp + pivot_row_length as usize;
        while rp < rp_end {
            let col = a_i[rp] as usize;
            rp += 1;
            // Skip dead columns.
            if col_is_dead(cols, col) {
                continue;
            }
            a_i[newrp as usize] = col as Int;
            newrp += 1;
            // Add new pivot row to column.
            a_i[(cols[col].start + cols[col].length) as usize] = pivot_row;
            cols[col].length += 1;

            // Retrieve score so far and add on pivot row's degree.
            // (we wait until here for this in case the pivot
            // row's degree was reduced due to mass elimination).
            let mut cur_score = cols[col].shared2.score() + pivot_row_degree;

            // Calculate the max possible score as the number of
            // external columns minus the 'k' value minus the
            // columns thickness.
            let max_score = ncol - k - cols[col].shared1.thickness();

            // Make the score the external degree of the union-of-rows.
            cur_score -= cols[col].shared1.thickness();

            // Make sure score is less or equal than the max score.
            cur_score = Int::min(cur_score, max_score);
            assert_debug!(cur_score >= 0);

            // Store updated score.
            cols[col].shared2 = ColShared2::Score(cur_score);

            // Place column back in degree list.

            assert_debug!(min_score >= 0);
            assert_debug!(min_score <= ncol);
            assert_debug!(cur_score >= 0);
            assert_debug!(cur_score <= ncol);
            assert_debug!(head[cur_score as usize] >= EMPTY);

            let next_col = head[cur_score as usize];
            cols[col].shared4 = ColShared4::DegreeNext(next_col);
            cols[col].shared3 = ColShared3::Prev(EMPTY);
            if next_col != EMPTY {
                cols[next_col as usize].shared3 = ColShared3::Prev(col as Int);
            }
            head[cur_score as usize] = col as Int;

            // See if this score is less than current min.
            min_score = Int::min(min_score, cur_score);
        }

        #[cfg(feature = "debug")]
        debug_deg_lists(nrow, ncol, rows, cols, head, min_score, n_col2 - k, maxdeg);

        // Resurrect the new pivot row.

        if pivot_row_degree > 0 {
            // Update pivot row length to reflect any cols that were killed
            // during super-col detection and mass elimination.
            rows[pivot_row as usize].start = pivot_row_start;
            rows[pivot_row as usize].length = newrp - pivot_row_start;
            assert_debug!(rows[pivot_row as usize].length > 0);

            rows[pivot_row as usize].shared1 = RowShared1::Degree(pivot_row_degree);
            rows[pivot_row as usize].shared2 = RowShared2::Mark(0);
            // Pivot row is no longer dead.

            debug1!(
                "Resurrect Pivot_row {} deg: {}",
                pivot_row,
                pivot_row_degree
            );
        }
    }

    // All principal columns have now been ordered.

    ngarbage
}

// The `find_ordering` routine has ordered all of the principal columns (the
// representatives of the supercolumns). The non-principal columns have not
// yet been ordered. This routine orders those columns by walking up the
// parent tree (a column is a child of the column which absorbed it). The
// final permutation vector is then placed in `p[0...n_col-1]`, with `p[0]`
// being the first column, and `p[n_col-1]` being the last. It doesn't look
// like it at first glance, but be assured that this routine takes time linear
// in the number of columns. Although not immediately obvious, the time
// taken by this routine is `O(n_col)`, that is, linear in the number of
// columns.
pub(crate) fn order_children(ncol: Int, cols: &mut [Col], p: &mut [Int]) {
    // Order each non-principal column.

    for i in 0..ncol as usize {
        // Find an un-ordered non-principal column.
        assert_debug!(col_is_dead(cols, i));

        if !col_is_dead_principal(cols, i) && cols[i].shared2.order() == EMPTY {
            let mut parent = i; // Index of column's parent.

            // Once found, find its principal parent.
            // do {
            //     parent = Col [parent].shared1.parent ;
            // } while (!COL_IS_DEAD_PRINCIPAL (parent)) ;
            loop {
                parent = cols[parent].shared1.parent() as usize;
                if col_is_dead_principal(cols, parent) {
                    break; // TODO: check
                }
            }

            // Now, order all un-ordered non-principal columns along path
            // to this parent. Collapse tree at the same time.
            let mut c = i;
            // Get order of parent.
            let mut order = cols[parent].shared2.order();

            // for ok := true; ok; ok = (cols[c].order() == empty) {
            loop {
                assert_debug!(cols[c].shared2.order() == EMPTY);

                // Order this column.
                cols[c].shared2 = ColShared2::Order(order);
                order += 1;
                // Collapse tree.
                cols[c].shared1 = ColShared1::Parent(parent as Int);

                // Get immediate parent of this column.
                c = cols[c].shared1.parent() as usize;

                // Continue until we hit an ordered column. There are
                // guaranteed not to be anymore unordered columns
                // above an ordered column.

                if cols[c].shared2.order() != EMPTY {
                    break; // TODO: check
                }
            }

            // Re-order the super_col parent to largest order for this group.
            cols[parent].shared2 = ColShared2::Order(order);
        }
    }

    // Generate the permutation.
    for c in 0..ncol {
        p[cols[c as usize].shared2.order() as usize] = c;
    }
}

// Detects supercolumns by finding matches between columns in the hash buckets.
// Check amongst columns in the set `A[row_start ... row_start + row_length-1]`.
// The columns under consideration are currently *not* in the degree lists,
// and have already been placed in the hash buckets.
//
// The hash bucket for columns whose hash function is equal to h is stored
// as follows:
//
// if `head[h]` is >= 0, then `head[h]` contains a degree list, so:
//
//     `head[h]` is the first column in degree bucket `h`.
//     `Col[head [h]].headhash` gives the first column in hash bucket `h`.
//
// otherwise, the degree list is empty, and:
//
//     `-(head [h] + 2)` is the first column in hash bucket `h`.
//
// For a column `c` in a hash bucket, `Col[c].shared3.prev` is NOT a "previous
// column" pointer. `Col[c].shared3.hash` is used instead as the hash number
// for that column. The value of `Col[c].shared4.hash_next` is the next column
// in the same hash bucket.
//
// Assuming no, or "few" hash collisions, the time taken by this routine is
// linear in the sum of the sizes (lengths) of each column whose score has
// just been computed in the approximate degree computation.
pub(crate) fn detect_super_cols(
    _ncol: Int,
    _rows: Option<&[Row]>,
    cols: &mut [Col],
    a_i: &[Int],
    head: &mut [Int],
    row_start: Int,
    row_length: Int,
) {
    // Consider each column in the row.

    let mut rp = row_start as usize;
    let rp_end = rp + row_length as usize;
    while rp < rp_end {
        let col = a_i[rp] as usize; // Column index in the row to check.
        rp += 1;
        if col_is_dead(cols, col) {
            continue;
        }

        // Get hash number for this column.
        let hash = cols[col].shared3.hash();
        assert_debug!(hash <= _ncol);

        // Get the first column in this hash bucket.

        let head_column = head[hash as usize]; // First column in hash bucket or degree list.

        // First column in hash bucket.
        let first_col = if head_column > EMPTY {
            cols[head_column as usize].shared3.head_hash()
        } else {
            -(head_column + 2)
        };

        // Consider each column in the hash bucket.

        // for super_c = first_col; super_c != empty; super_c = cols[super_c].hashNext() {  // TODO: check
        let mut super_c = first_col; // Column index of the column to absorb into.
        while super_c != EMPTY {
            assert_debug!(col_is_alive(cols, super_c as usize));
            assert_debug!(cols[super_c as usize].shared3.hash() == hash);

            let length = cols[super_c as usize].length; // Length of column super_c.

            // prev_c is the column preceding column c in the hash bucket.
            let mut prev_c = super_c;

            // Compare super_c with all columns after it.

            // for c = cols[super_c].hashNext(); c != empty; c = cols[c].hashNext() {  TODO: check
            let mut c = cols[super_c as usize].shared4.hash_next();
            while c != EMPTY {
                assert_debug!(c != super_c);
                assert_debug!(col_is_alive(cols, c as usize));
                assert_debug!(cols[c as usize].shared3.hash() == hash);

                // Not identical if lengths or scores are different.
                if cols[c as usize].length != length
                    || cols[c as usize].shared2.score() != cols[super_c as usize].shared2.score()
                {
                    prev_c = c;

                    c = cols[c as usize].shared4.hash_next();
                    continue;
                }

                // Compare the two columns.
                let mut cp1 = cols[super_c as usize].start as usize; // Column pointer for column super_c.
                let mut cp2 = cols[c as usize].start as usize; // Column pointer for column c.

                let mut i = 0;
                while i < length {
                    // The columns are "clean" (no dead rows).
                    #[cfg(feature = "debug")]
                    if let Some(rows) = _rows {
                        assert_debug!(row_is_alive(rows, a_i[cp1] as usize));
                        assert_debug!(row_is_alive(rows, a_i[cp2] as usize));
                    }
                    // Row indices will same order for both supercols,
                    // no gather scatter nessasary.
                    let _cp1 = cp1;
                    cp1 += 1;
                    let _cp2 = cp2;
                    cp2 += 1;
                    if a_i[_cp1] != a_i[_cp2] {
                        break;
                    }
                    i += 1;
                }

                // The two columns are different if the for-loop "broke".
                if i != length {
                    prev_c = c;

                    c = cols[c as usize].shared4.hash_next();
                    continue;
                }

                // Got it! Two columns are identical.

                assert_debug!(
                    cols[c as usize].shared2.score() == cols[super_c as usize].shared2.score()
                );

                cols[super_c as usize].shared1 = ColShared1::Thickness(
                    cols[super_c as usize].shared1.thickness()
                        + cols[c as usize].shared1.thickness(),
                );
                cols[c as usize].shared1 = ColShared1::Parent(super_c);
                kill_non_principal_col(cols, c as usize);
                // Order c later, in orderChildren().
                cols[c as usize].shared2 = ColShared2::Order(EMPTY);
                // Remove c from hash bucket.
                cols[prev_c as usize].shared4 =
                    ColShared4::HashNext(cols[c as usize].shared4.hash_next());

                c = cols[c as usize].shared4.hash_next();
            }
            super_c = cols[super_c as usize].shared4.hash_next();
        }

        // Empty this hash bucket.
        if head_column > EMPTY {
            // Corresponding degree list "hash" is not empty.
            cols[head_column as usize].shared3 = ColShared3::HeadHash(EMPTY);
        } else {
            // Corresponding degree list "hash" is empty.
            head[hash as usize] = EMPTY;
        }
    }
}

// Defragments and compacts columns and rows in the workspace A. Used when
// all avaliable memory has been used while performing row merging. Returns
// the index of the first free position in A, after garbage collection. The
// time taken by this routine is linear in the size of the array A, which is
// itself linear in the number of nonzeros in the input matrix.
pub(crate) fn garbage_collection(
    n_row: Int,
    n_col: Int,
    rows: &mut [Row],
    cols: &mut [Col],
    a_i: &mut [Int],
    pfree: Int,
) -> Int {
    #[cfg(feature = "debug")]
    let mut debug_rows = {
        debug2!("Defrag..");
        for psrc in 0..pfree as usize {
            assert_debug!(a_i[psrc] >= 0);
        }
        0
    };

    // Defragment the columns.

    let mut pdest = 0;
    for c in 0..n_col as usize {
        if col_is_alive(cols, c) {
            let mut psrc = cols[c].start;

            // Move and compact the column.
            assert_debug!(pdest <= psrc);

            cols[c].start = pdest - 0;
            let length = cols[c].length;
            for _j in 0..length {
                let r = a_i[psrc as usize] as usize;
                psrc += 1;
                if row_is_alive(rows, r) {
                    a_i[pdest as usize] = r as Int;
                }
            }
            cols[c].length = pdest - cols[c].start;
        }
    }

    // Prepare to defragment the rows.

    for r in 0..n_row as usize {
        if row_is_dead(rows, r) || (rows[r].length == 0) {
            // This row is already dead, or is of zero length. Cannot compact
            // a row of zero length, so kill it. NOTE: in the current version,
            // there are no zero-length live rows. Kill the row (for the first
            // time, or again) just to be safe.
            kill_row(rows, r);
        } else {
            // Save first column index in Row [r].shared2.first_column.
            let psrc = rows[r].start;
            rows[r].shared2 = RowShared2::FirstColumn(a_i[psrc as usize]);
            assert_debug!(row_is_alive(rows, r));

            // Flag the start of the row with the one's complement of row.
            a_i[psrc as usize] = ones_complement(r as Int);

            #[cfg(feature = "debug")]
            {
                debug_rows += 1;
            }
        }
    }

    // Defragment the rows.

    let mut psrc = pdest;
    while psrc < pfree {
        // Find a negative number ... the start of a row.
        let _psrc = psrc;
        psrc += 1;
        if a_i[_psrc as usize] < 0 {
            psrc -= 1;
            // Get the row index.
            let r = ones_complement(a_i[psrc as usize]) as usize;
            assert_debug!(/*r >= 0 &&*/ r < n_row as usize);

            // Restore first column index.
            a_i[psrc as usize] = rows[r].shared2.first_column();
            assert_debug!(row_is_alive(rows, r));
            assert_debug!(rows[r].length > 0);

            // Move and compact the row.
            assert_debug!(pdest <= psrc);

            rows[r].start = pdest - 0;
            let length = rows[r].length;
            for _j in 0..length {
                let c = a_i[psrc as usize] as usize;
                psrc += 1;
                if col_is_alive(cols, c) {
                    a_i[pdest as usize] = c as Int;
                    pdest += 1;
                }
            }
            rows[r].length = pdest - rows[r].start;
            #[cfg(feature = "debug")]
            {
                assert_debug!(rows[r].length > 0);
                debug_rows -= 1;
            }
        }
    }

    // Ensure we found all the rows.
    assert_debug!(debug_rows == 0);

    // Return the new value of pfree.

    pdest - 0
}

// Clears the Row[].shared2.mark array, and returns the new tag_mark.
fn clear_mark(mut tag_mark: Int, max_mark: Int, nrow: Int, row: &mut [Row]) -> Int {
    if tag_mark <= 0 || tag_mark >= max_mark {
        for r in 0..nrow as usize {
            if row_is_alive(row, r) {
                row[r].shared2 = RowShared2::Mark(0)
            }
        }
        tag_mark = 1;
    }
    tag_mark
}

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
    A: &mut [i32],
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

    for col in 0..ncol {
        let mut lastRow = -1;

        let mut cp = p[col as usize] as usize;
        let cpend = p[(col + 1) as usize] as usize;

        while cp < cpend {
            let row = A[cp];
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

            if row <= lastRow || rows[row as usize].mark() == col {
                // Row index are unsorted or repeated (or both), thus col
                // is jumbled. This is a notice, not an error condition.
                stats[STATUS] = OK_BUT_JUMBLED;
                stats[INFO1] = col;
                stats[INFO2] = row;
                stats[INFO3] += 1;
                debug1!("colamd: row {} col {} unsorted/duplicate", row, col);
            }

            if rows[row as usize].mark() != col {
                rows[row as usize].length += 1;
            } else {
                // This is a repeated entry in the column, it will be removed.
                cols[col as usize].length -= 1;
            }

            // Mark the row as having been seen in this column.
            rows[row as usize].set_mark(col);

            lastRow = row;
        }
    }

    // Compute row pointers.

    // Row form of the matrix starts directly after the column
    // form of matrix in A.
    rows[0].start = p[ncol as usize];
    rows[0].set_p(rows[0].start);
    rows[0].set_mark(-1);
    for row in 1..nrow as usize {
        rows[row].start = rows[row - 1].start + rows[row - 1].length;
        rows[row].set_p(rows[row].start);
        rows[row].set_mark(-1);
    }

    // Create row form.

    if stats[STATUS] == OK_BUT_JUMBLED {
        // If cols jumbled, watch for repeated row indices.
        for col in 0..ncol {
            let mut cp = p[col as usize];
            let cpend = p[(col + 1) as usize];

            while cp < cpend {
                let row = A[cp as usize] as usize;
                cp += 1;

                if rows[row].mark() != col {
                    A[rows[row].p() as usize] = col;
                    rows[row].set_p(rows[row].p() + 1);
                    rows[row].set_mark(col);
                }
            }
        }
    } else {
        // If cols not jumbled, we don't need the mark (this is faster).
        for col in 0..ncol {
            let mut cp = p[col as usize] as usize;
            let cpend = p[(col + 1) as usize] as usize;
            while cp < cpend {
                A[rows[A[cp] as usize].p() as usize] = col;
                rows[A[cp] as usize].set_p(rows[A[cp] as usize].p() + 1);
                cp += 1;
            }
        }
    }

    // Clear the row marks and set row degrees.

    for row in 0..nrow as usize {
        rows[row].set_mark(0);
        rows[row].set_degree(rows[row].length);
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
        for col in 1..ncol as usize {
            // Note that the lengths here are for pruned columns, i.e.
            // no duplicate row indices will exist for these columns.
            cols[col].start = cols[col - 1].start + cols[col - 1].length;
            p[col] = cols[col].start;
        }

        // Re-create col form.

        for row in 0..nrow {
            let mut rp = rows[row as usize].start as usize;
            let rpend = rp + rows[row as usize].length as usize;
            while rp < rpend {
                A[p[A[rp] as usize] as usize] = row;
                p[A[rp] as usize] += 1;
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
    rows: &mut [Row],
    cols: &mut [Col],
    A: &mut [i32],
    head: &mut [i32],
    knobs: [f64; KNOBS],
    pnrow2: &mut i32,
    pncol2: &mut i32,
    pmaxDeg: &mut i32,
) {
    let c: i32; // A column index.
    let r: i32; // A row index.
                // let row: usize; // A row index.
                // let cp: usize; // A column pointer.
                // let cpend: usize; // A pointer to the end of a column.
                // let newcp: usize; // New column pointer.
                // let colLength: i32; // Length of pruned column.
                // let score: i32; // Current column score.
    let denseRowCount: i32; // Remove rows with more entries than this.
    let denseColCount: i32; // Remove cols with more entries than this.
                            // let minScore: i32; // Smallest column score.
                            // let maxDeg: i32; // Maximum row degree.
                            // let nextCol: i32; // Used to add to degree list.

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
    let mut maxDeg = 0; // Maximum row degree.
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
            cols[c].set_order(ncol2);
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
        if deg > denseColCount {
            // This is a dense column, kill and order it last.
            ncol2 -= 1;
            cols[c].set_order(ncol2);
            // Decrement the row degrees.
            let mut cp = cols[c].start as usize; // Column pointer.
            let cpend = cp + cols[c].length as usize; // Pointer to the end of the column.
            while cp < cpend {
                rows[A[cp] as usize].set_degree(rows[A[cp] as usize].degree() - 1);
                cp += 1;
            }
            kill_principal_col(cols, c);
        }
    }
    debug1!("colamd: Dense and null columns killed: {}", ncol - ncol2);

    // Kill dense and empty rows.
    for r in 0..nrow as usize {
        let deg = rows[r].degree(); // Degree of row.
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

    // Compute initial column scores.

    // At this point the row degrees are accurate. They reflect the number
    // of "live" (non-dense) columns in each row. No empty rows exist.
    // Some "live" columns may contain only dead rows, however. These are
    // pruned in the code below.

    // Now find the initial matlab score for each column.
    // for c = ncol - 1; c >= 0; c-- {
    for c in (0..ncol as usize).rev() {
        // TODO: check
        // Skip dead column.
        if col_is_dead(cols, c) {
            continue;
        }
        let mut score = 0; // Current column score.
        let mut cp = cols[c].start as usize;
        let mut newcp = cp; // New column pointer.
        let cpend = cp + cols[c].length as usize;
        while cp < cpend {
            // Get a row.
            let row = A[cp] as usize; // row index
            cp += 1;
            // Skip if dead.
            if row_is_dead(rows, row) {
                continue;
            }
            // Compact the column.
            A[newcp] = row as i32;
            newcp += 1;
            // Add row's external degree.
            score += rows[row].degree() - 1;
            // Guard against integer overflow.
            score = i32::min(score, ncol);
        }
        // Determine pruned column length.
        let colLength = (newcp as i32 - cols[c].start); // Length of pruned column.
        if colLength == 0 {
            // A newly-made null column (all rows in this col are "dense"
            // and have already been killed).
            debug2!("Newly null killed: {}", c);
            ncol2 -= 1;
            cols[c].set_order(ncol2);
            kill_principal_col(cols, c);
        } else {
            // Set column length and set score.
            debug_assert!(score >= 0);
            debug_assert!(score <= ncol);
            cols[c].length = colLength;
            cols[c].set_score(score);
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
    debug_structures(nrow, ncol, rows, cols, A, ncol2);

    // Initialize degree lists.

    #[cfg(feature = "debug")]
    let mut debugCount = 0;

    // Clear the hash buckets.
    for c in 0..=ncol as usize {
        head[c] = empty;
    }
    let mut minScore = ncol; // Smallest column score.

    // Place in reverse order, so low column indices are at the front
    // of the lists. This is to encourage natural tie-breaking.
    // for c = ncol - 1; c >= 0; c-- {
    for c in (0..ncol as usize).rev() {
        // Only add principal columns to degree lists.
        if col_is_alive(cols, c) {
            debug4!(
                "place {} score {} minscore {} ncol {}",
                c,
                cols[c].score(),
                minScore,
                ncol
            );

            // Add columns score to DList.

            let score = cols[c].score();

            debug_assert!(minScore >= 0);
            debug_assert!(minScore <= ncol);
            debug_assert!(score >= 0);
            debug_assert!(score <= ncol);
            debug_assert!(head[score as usize] >= empty);

            // Now add this column to dList at proper score location.
            let nextCol = head[score as usize]; // Used to add to degree list.
            cols[c].set_prev(empty);
            cols[c].set_degree_next(nextCol);

            // if there already was a column with the same score, set its
            // previous pointer to this new column.
            if nextCol != empty {
                cols[nextCol as usize].set_prev(c as i32);
            }
            head[score as usize] = c as i32;

            // See if this score is less than current min.
            minScore = i32::min(minScore, score);

            #[cfg(feature = "debug")]
            {
                debugCount += 1;
            }
        }
    }

    debug1!(
        "colamd: Live cols {} out of {}, non-princ: {}",
        debugCount,
        ncol,
        ncol - debugCount
    );
    #[cfg(feature = "debug")]
    {
        assert!(debugCount == ncol2);
        debugDegLists(nrow, ncol, rows, cols, head, minScore, ncol2, maxDeg);
    }

    // Return number of remaining columns, and max row degree.

    *pncol2 = ncol2;
    *pnrow2 = nrow2;
    *pmaxDeg = maxDeg;
}

// Order the principal columns of the supercolumn form of the matrix
// (no supercolumns on input). Uses a minimum approximate column minimum
// degree ordering method.
//
// `Alen` is the size of `A` and must be `2*nnz + n_col` or larger.
// `n_col2` is the number of remaining columns to order.
// `max_deg` is the maximum row degree.
// `pfree` is the index of first free slot (`2*nnz` on entry).
// Returns the number of garbage collections.
pub(crate) fn find_ordering(
    nrow: i32,
    ncol: i32,
    Alen: i32,
    rows: &mut [Row],
    cols: &mut [Col],
    A: &mut [i32],
    head: &mut [i32],
    n_col2: i32,
    mut maxdeg: i32,
    mut pfree: i32,
    aggressive: bool,
) -> i32 {
    let k: i32; // Current pivot ordering step.
                // let pivotCol: i32; // Current pivot column.
                // let cp: usize; // A column pointer.
                // let rp: usize; // A row pointer.
                // let pivotRow: i32; // Current pivot row.
                // let newcp: i32; // Modified column pointer.
                // let newrp: i32; // Modified row pointer.
                // let pivotRowStart: i32; // Pointer to start of pivot row.
                // let pivotRowDegree: i32; // Number of columns in pivot row.
                // let pivotRowLength: i32; // Number of supercolumns in pivot row.
                // let pivotColScore: i32; // Score of pivot column.
                // let neededMemory: i32; // Free space needed for pivot row.
                // let cpend: usize; // Pointer to the end of a column.
                // let rpend: usize; // Pointer to the end of a row.
                // let row: usize; // A row index.
                // let col: usize; // A column index.
                // let maxScore: i32; // Maximum possible score.
                // let curScore: i32; // Score of current column.
                // let hash: usize; // Hash value for supernode detection.
                // let headColumn: i32; // Head of hash bucket.
                // let firstCol: i32; // First column in hash bucket.
                // let tagMark: i32; // Marker value for mark array.
                // let rowMark: i32; // Row [row].shared2.mark
                // let setDifference: i32; /// Set difference size of row with pivot row.
                // let minScore: i32; // Smallest column score.
                // let colThickness: i32; // "thickness" (no. of columns in a supercol)
    let maxMark: i32; // Maximum value of tag_mark.
                      // let pivotColThickness: i32; // Number of columns represented by pivot col.
                      // let prevCol: i32; // Used by Dlist operations.
                      // let nextCol: i32; // Used by Dlist operations.
                      // let ngarbage: i32; // Number of garbage collections performed.

    #[cfg(feature = "debug")]
    let debugd: i32; // Debug loop counter.
    #[cfg(feature = "debug")]
    let mut debugStep = 0; // Debug loop counter.

    // Initialization and clear mark.

    let maxMark = i32::MAX - ncol;
    let mut tagMark = clear_mark(0, maxMark, nrow, rows); // Marker value for mark array.
    let mut minScore = 0; // Smallest column score.
    let mut ngarbage = 0; // Number of garbage collections performed.
    debug1!("colamd: Ordering, n_col2={}", n_col2);

    // Order the columns.

    let mut k = 0; // 'k' is incremented below
    while k < n_col2 {
        #[cfg(feature = "debug")]
        {
            if debugStep % 100 == 0 {
                debug2!("\n...      Step k: {} out of n_col2: {}", k, n_col2);
            } else {
                debug3!("\n----------Step k: {} out of n_col2: {}", k, n_col2);
            }
            debugStep += 1;
            debugDegLists(nrow, ncol, rows, cols, head, minScore, n_col2 - k, maxdeg);
            debugMatrix(nrow, ncol, rows, cols, A);
        }
        // ndebug

        // Select pivot column, and order it.

        #[cfg(feature = "debug")]
        {
            // Make sure degree list isn't empty.
            assert!(minScore >= 0);
            assert!(minScore <= ncol);
            assert!(head[minScore] >= empty);

            for debugd in 0..minScore {
                assert!(head[debugd] == empty);
            }
        }
        // ndebug

        // Get pivot column from head of minimum degree list.
        while head[minScore as usize] == empty && minScore < ncol {
            minScore += 1;
        }
        let pivotCol = head[minScore as usize] as usize; // Current pivot column.
        #[cfg(feature = "debug")]
        assert(pivotCol >= 0 && pivotCol <= ncol as usize);

        let nextCol = cols[pivotCol].degree_next();
        head[minScore as usize] = nextCol;
        if nextCol != empty {
            cols[nextCol as usize].set_prev(empty);
        }

        #[cfg(feature = "debug")]
        {
            assert!(colIsAlive(cols, pivotCol));
            debug3!("Pivot col: {}\n", pivotCol);
        }

        // Remember score for defrag check.
        let pivotColScore = cols[pivotCol].score(); // Score of pivot column.

        // The pivot column is the kth column in the pivot order.
        cols[pivotCol].set_order(k);

        // Increment order count by column thickness.
        let pivotColThickness = cols[pivotCol].thickness(); // Number of columns represented by pivot col.
        k += pivotColThickness;
        #[cfg(feature = "debug")]
        {
            assert!(pivotColThickness > 0);
            debug3!("Pivot col: {} thick {}", pivotCol, pivotColThickness);
        }

        // Garbage_collection, if necessary.

        let neededMemory = i32::min(pivotColScore, ncol - k); // Free space needed for pivot row.
        if pfree + neededMemory >= Alen {
            pfree = garbage_collection(nrow, ncol, rows, cols, A, pfree);
            ngarbage += 1;
            // After garbage collection we will have enough.
            #[cfg(feature = "debug")]
            assert!(pfree + neededMemory < Alen);

            // Garbage collection has wiped out the Row[].shared2.mark array.
            tagMark = clear_mark(0, maxMark, nrow, rows);

            #[cfg(feature = "debug")]
            debug_matrix(nrow, ncol, rows, cols, A);
        }

        // Compute pivot row pattern.

        // Get starting location for this new merged row.
        let pivotRowStart = pfree; // Pointer to start of pivot row.

        // Initialize new row counts to zero.
        let mut pivotRowDegree = 0; // Number of columns in pivot row.

        // Tag pivot column as having been visited so it isn't included
        // in merged pivot row.
        cols[pivotCol].set_thickness(-pivotColThickness);

        // Pivot row is the union of all rows in the pivot column pattern.
        let mut cp = cols[pivotCol].start as usize;
        let cpend = cp + cols[pivotCol].length as usize; // Pointer to the end of the column.
        while cp < cpend {
            // Get a row.
            let row = A[cp] as usize;
            cp += 1;
            debug4!("Pivot col pattern {] {}", rowIsAlive(rows, row), row);

            // Skip if row is dead.
            if row_is_alive(rows, row) {
                let mut rp = rows[row].start as usize;
                let rpend = rp + rows[row].length as usize; // Pointer to the end of the row.
                while rp < rpend {
                    // Get a column.
                    let col = A[rp] as usize;
                    rp += 1;
                    // Add the column, if alive and untagged.
                    let colThickness = cols[col].thickness(); // "thickness" (no. of columns in a supercol)
                    if colThickness > 0 && col_is_alive(cols, col) {
                        // Tag column in pivot row.
                        cols[col].set_thickness(-colThickness);
                        #[cfg(feature = "debug")]
                        assert!(pfree < Alen);

                        // Place column in pivot row.
                        A[pfree as usize] = col as i32;
                        pfree += 1;
                        pivotRowDegree += colThickness;
                    }
                }
            }
        }

        // Clear tag on pivot column.
        cols[pivotCol].set_thickness(pivotColThickness);
        maxdeg = i32::max(maxdeg, pivotRowDegree);

        #[cfg(feature = "debug")]
        {
            debug3!("check2");
            debug_mark(nrow, rows, tagMark, maxMark);
        }

        // Kill all rows used to construct pivot row.

        // Also kill pivot row, temporarily.
        let mut cp = cols[pivotCol].start as usize;
        let cpend = cp + cols[pivotCol].length as usize; // Pointer to the end of the column.
        while cp < cpend {
            // May be killing an already dead row.
            let row = A[cp] as usize;
            cp += 1;
            debug3!("Kill row in pivot col: {}", row);

            kill_row(rows, row)
        }

        // Select a row index to use as the new pivot row.

        let pivotRowLength = pfree - pivotRowStart; // Number of supercolumns in pivot row.
        let pivotRow = if pivotRowLength > 0 {
            // Pick the "pivot" row arbitrarily (first row in col).
            let pivotRow = A[cols[pivotCol].start as usize];
            debug3!("Pivotal row is {}", pivotRow);
            pivotRow
        } else {
            // There is no pivot row, since it is of zero length.
            let pivotRow = empty;
            #[cfg(feature = "debug")]
            assert!(pivotRowLength == 0);
            pivotRow
        };
        #[cfg(feature = "debug")]
        assert!(cols[pivotCol].length > 0 || pivotRowLength == 0);

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
        let mut rp = pivotRowStart as usize;
        let rpend = rp + pivotRowLength as usize;
        while rp < rpend {
            let col = A[rp] as usize;
            rp += 1;
            #[cfg(feature = "debug")]
            {
                assert!(colIsAlive(cols, col) && col != pivotCol);
                debug3!("Col: {}", col);
            }

            // Clear tags used to construct pivot row pattern.
            let colThickness = -cols[col].thickness();
            #[cfg(feature = "debug")]
            assert!(colThickness > 0);

            cols[col].set_thickness(colThickness);

            // Remove column from degree list.

            let curScore = cols[col].score(); // Score of current column.
            let prevCol = cols[col].prev();
            let nextCol = cols[col].degree_next();
            #[cfg(feature = "debug")]
            {
                assert(curScore >= 0);
                assert(curScore <= ncol);
                assert(curScore >= empty);
            }
            if prevCol == empty {
                head[curScore as usize] = nextCol;
            } else {
                cols[prevCol as usize].set_degree_next(nextCol);
            }
            if nextCol != empty {
                cols[nextCol as usize].set_prev(prevCol);
            }

            // Scan the column.

            let mut cp = cols[col].start as usize;
            let cpend = cp + cols[col].length as usize; // Pointer to the end of a column.
            while cp < cpend {
                // Get a row.
                let row = A[cp] as usize;
                cp += 1;
                let rowMark = rows[row].mark();
                // Skip if dead.
                if row_is_marked_dead(rowMark) {
                    continue;
                }
                #[cfg(feature = "debug")]
                assert!(row != pivotRow);

                let mut setDifference = rowMark - tagMark;
                // Check if the row has been seen yet.
                if setDifference < 0 {
                    #[cfg(feature = "debug")]
                    assert!(rows[row].degree() <= maxdeg);

                    setDifference = rows[row].degree();
                }
                // Subtract column thickness from this row's set difference.
                setDifference -= colThickness;
                #[cfg(feature = "debug")]
                assert!(setDifference >= 0);

                // Absorb this row if the set difference becomes zero.
                if setDifference == 0 && aggressive {
                    debug3!("aggressive absorption. Row: %d", row);

                    kill_row(rows, row);
                } else {
                    // Save the new mark.
                    rows[row].set_mark(setDifference + tagMark);
                }
            }
        }

        #[cfg(feature = "debug")]
        debugDegLists(
            nrow,
            ncol,
            rows,
            cols,
            head,
            minScore,
            n_col2 - k - pivotRowDegree,
            maxdeg,
        );

        // Add up set differences for each column.

        debug3!("** Adding set differences phase. **");

        // For each column in pivot row.
        let mut rp = pivotRowStart as usize;
        let rpend = rp + pivotRowLength as usize;
        while rp < rpend {
            // Get a column.
            let col = A[rp] as usize;
            rp += 1;
            #[cfg(feature = "debug")]
            assert(colIsAlive(cols, col) && col != pivotCol);

            let mut hash: usize = 0; // Hash value for supernode detection.
            let mut curScore = 0; // Score of current column.
            let mut cp = cols[col].start as usize;
            // Compact the column.
            let mut newcp = cp as i32; // Modified column pointer.
            let cpend = cp + cols[col].length as usize;

            debug4!("Adding set diffs for Col: {}.", col);

            while cp < cpend {
                // Get a row.
                let row = A[cp] as usize;
                cp += 1;
                #[cfg(feature = "debug")]
                assert!(row >= 0 && row < nrow);

                let rowMark = rows[row].mark();
                // Skip if dead.
                if row_is_marked_dead(rowMark) {
                    debug4!(" Row {}, dead", row);
                    continue;
                }
                debug4!(" Row {}, set diff {}", row, rowMark - tagMark);
                #[cfg(feature = "debug")]
                assert!(rowMark >= tagMark);

                // Compact the column.
                A[newcp as usize] = row as i32;
                newcp += 1;
                // Compute hash function.
                hash += row as usize;
                // Add set difference.
                curScore += rowMark - tagMark;
                // Integer overflow...
                curScore = i32::min(curScore, ncol);
            }

            // Recompute the column's length.
            cols[col].length = (newcp - cols[col].start);

            // Further mass elimination.

            if cols[col].length == 0 {
                debug4!("further mass elimination. Col: {}", col);

                // Nothing left but the pivot row in this column.
                kill_principal_col(cols, col);
                pivotRowDegree -= cols[col].thickness();
                #[cfg(feature = "debug")]
                assert(pivotRowDegree >= 0);

                // Order it.
                cols[col].set_order(k);
                // Increment order count by column thickness.
                k += cols[col].thickness();
            } else {
                // Prepare for supercolumn detection.

                debug4!("Preparing supercol detection for Col: {}.", col);

                // Save score so far.
                cols[col].set_score(curScore);

                // Add column to hash table, for supercolumn detection.
                hash %= (ncol as usize) + 1;

                debug4!(" Hash = {}, n_col = {}.", hash, ncol);
                #[cfg(feature = "debug")]
                assert!((hash as i32) <= ncol);

                let headColumn = head[hash]; // Head of hash bucket.
                                             // First column in hash bucket.
                let firstCol = if headColumn > empty {
                    // Degree list "hash" is non-empty, use prev (shared3) of
                    // first column in degree list as head of hash bucket.
                    let firstCol = cols[headColumn as usize].headhash();
                    cols[headColumn as usize].set_headhash(col as i32);
                    firstCol
                } else {
                    // degree list "hash" is empty, use head as hash bucket.
                    let firstCol = -(headColumn + 2);
                    head[hash] = -1 * (col + 2) as i32;
                    firstCol
                };
                cols[col].set_hash_next(firstCol);

                // Save hash function in Col [col].shared3.hash
                cols[col].set_hash(hash as i32);
                #[cfg(feature = "debug")]
                assert!(colIsAlive(cols, col));
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
            A,
            head,
            pivotRowStart,
            pivotRowLength,
        );
        #[cfg(not(feature = "debug"))]
        detect_super_cols(0, None, cols, A, head, pivotRowStart, pivotRowLength);

        // Kill the pivotal column.

        kill_principal_col(cols, pivotCol);

        // Clear mark.

        tagMark = clear_mark(tagMark + maxdeg + 1, maxMark, nrow, rows);

        #[cfg(feature = "debug")]
        {
            debug3!("check3");
            debugMark(nrow, rows, tagMark, maxMark);
        }

        // Finalize the new pivot row, and column scores.

        debug3!("** Finalize scores phase. **");

        // For each column in pivot row.
        let mut rp = pivotRowStart as usize;
        // Compact the pivot row.
        let mut newrp = rp as i32; // Modified row pointer.
        let rpend = rp + pivotRowLength as usize;
        while rp < rpend {
            let col = A[rp] as usize;
            rp += 1;
            // Skip dead columns.
            if col_is_dead(cols, col) {
                continue;
            }
            A[newrp as usize] = col as i32;
            newrp += 1;
            // Add new pivot row to column.
            A[(cols[col].start + cols[col].length) as usize] = pivotRow;
            cols[col].length += 1;

            // Retrieve score so far and add on pivot row's degree.
            // (we wait until here for this in case the pivot
            // row's degree was reduced due to mass elimination).
            let mut curScore = cols[col].score() + pivotRowDegree;

            // Calculate the max possible score as the number of
            // external columns minus the 'k' value minus the
            // columns thickness.
            let maxScore = ncol - k - cols[col].thickness();

            // Make the score the external degree of the union-of-rows.
            curScore -= cols[col].thickness();

            // Make sure score is less or equal than the max score.
            curScore = i32::min(curScore, maxScore);
            #[cfg(feature = "debug")]
            assert!(curScore >= 0);

            // Store updated score.
            cols[col].set_score(curScore);

            // Place column back in degree list.

            #[cfg(feature = "debug")]
            {
                assert!(minScore >= 0);
                assert!(minScore <= ncol);
                assert!(curScore >= 0);
                assert!(curScore <= ncol);
                assert!(head[curScore] >= empty);
            }
            let nextCol = head[curScore as usize];
            cols[col].set_degree_next(nextCol);
            cols[col].set_prev(empty);
            if nextCol != empty {
                cols[nextCol as usize].set_prev(col as i32);
            }
            head[curScore as usize] = col as i32;

            // See if this score is less than current min.
            minScore = i32::min(minScore, curScore);
        }

        #[cfg(feature = "debug")]
        debugDegLists(nrow, ncol, rows, cols, head, minScore, n_col2 - k, maxdeg);

        // Resurrect the new pivot row.

        if pivotRowDegree > 0 {
            // Update pivot row length to reflect any cols that were killed
            // during super-col detection and mass elimination.
            rows[pivotRow as usize].start = pivotRowStart;
            rows[pivotRow as usize].length = newrp - pivotRowStart;
            #[cfg(feature = "debug")]
            assert!(rows[pivotRow as usize].length > 0);

            rows[pivotRow as usize].set_degree(pivotRowDegree);
            rows[pivotRow as usize].set_mark(0);
            // Pivot row is no longer dead.

            debug1!("Resurrect Pivot_row %d deg: %d", pivotRow, pivotRowDegree);
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
pub(crate) fn order_children(ncol: i32, cols: &mut [Col], p: &mut [i32]) {
    // let parent: usize; // Index of column's parent.
    // let order: i32; // Column's order.

    // Order each non-principal column.

    for i in 0..ncol as usize {
        // Find an un-ordered non-principal column.
        #[cfg(feature = "debug")]
        assert!(colIsDead(cols, i));

        if !col_is_dead_principal(cols, i) && cols[i].order() == empty {
            let mut parent = i; // Index of column's parent.

            // Once found, find its principal parent.
            // do {
            //     parent = Col [parent].shared1.parent ;
            // } while (!COL_IS_DEAD_PRINCIPAL (parent)) ;
            loop {
                parent = cols[parent].parent() as usize;
                if col_is_dead_principal(cols, parent) {
                    break; // TODO: check
                }
            }

            // Now, order all un-ordered non-principal columns along path
            // to this parent. Collapse tree at the same time.
            let mut c = i;
            // Get order of parent.
            let mut order = cols[parent].order();

            // for ok := true; ok; ok = (cols[c].order() == empty) {
            loop {
                #[cfg(feature = "debug")]
                assert!(cols[c].order() == empty);

                // Order this column.
                cols[c].set_order(order);
                order += 1;
                // Collapse tree.
                cols[c].set_parent(parent as i32);

                // Get immediate parent of this column.
                c = cols[c].parent() as usize;

                // Continue until we hit an ordered column. There are
                // guaranteed not to be anymore unordered columns
                // above an ordered column.

                if cols[c].order() != empty {
                    break; // TODO: check
                }
            }

            // Re-order the super_col parent to largest order for this group.
            cols[parent].set_order(order);
        }
    }

    // Generate the permutation.
    for c in 0..ncol {
        p[cols[c as usize].order() as usize] = c;
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
    ncol: i32,
    rows: Option<&[Row]>,
    cols: &mut [Col],
    A: &[i32],
    head: &mut [i32],
    rowStart: i32,
    rowLength: i32,
) {
    // let hash: i32; // Hash value for a column.
    let rp: usize; // Pointer to a row.
                   // let c: i32; // A column index.
                   // let superc: i32; // Column index of the column to absorb into.
                   // let cp1: usize; // Column pointer for column super_c.
                   // let cp2: usize; // Column pointer for column c.
                   // let length: i32; // Length of column super_c.
                   // let prevc: i32; // Column preceding c in hash bucket.
                   // let i: i32; // Loop counter.
    let rpend: usize; // Pointer to the end of the row.
                      // let col: usize; // A column index in the row to check.
                      // let headColumn: i32; // First column in hash bucket or degree list.
                      // let firstCol: i32; // First column in hash bucket.

    // Consider each column in the row.

    let mut rp = rowStart as usize;
    let rpend = rp + rowLength as usize;
    while rp < rpend {
        let col = A[rp] as usize; // Column index in the row to check.
        rp += 1;
        if col_is_dead(cols, col) {
            continue;
        }

        // Get hash number for this column.
        let hash = cols[col].hash();
        #[cfg(feature = "debug")]
        assert!(hash <= ncol);

        // Get the first column in this hash bucket.

        let headColumn = head[hash as usize]; // First column in hash bucket or degree list.

        // First column in hash bucket.
        let firstCol = if headColumn > empty {
            cols[headColumn as usize].headhash()
        } else {
            -(headColumn + 2)
        };

        // Consider each column in the hash bucket.

        // for superc = firstCol; superc != empty; superc = cols[superc].hashNext() {  // TODO: check
        let mut superc = firstCol; // Column index of the column to absorb into.
        while superc != empty {
            #[cfg(feature = "debug")]
            {
                assert!(colIsAlive(cols, superc));
                assert!(cols[superc].hash() == hash);
            }
            let length = cols[superc as usize].length; // Length of column super_c.

            // prev_c is the column preceding column c in the hash bucket.
            let mut prevc = superc;

            // Compare super_c with all columns after it.

            // for c = cols[superc].hashNext(); c != empty; c = cols[c].hashNext() {  TODO: check
            let mut c = cols[superc as usize].hash_next();
            while c != empty {
                #[cfg(feature = "debug")]
                {
                    assert!(c != superc);
                    assert!(colIsAlive(cols, c));
                    assert!(cols[c].hash() == hash);
                }

                // Not identical if lengths or scores are different.
                if cols[c as usize].length != length
                    || cols[c as usize].score() != cols[superc as usize].score()
                {
                    prevc = c;

                    c = cols[c as usize].hash_next();
                    continue;
                }

                // Compare the two columns.
                let mut cp1 = cols[superc as usize].start as usize; // Column pointer for column super_c.
                let mut cp2 = cols[c as usize].start as usize; // Column pointer for column c.

                let mut i = 0;
                while i < length {
                    // The columns are "clean" (no dead rows).
                    #[cfg(feature = "debug")]
                    {
                        assert!(rowIsAlive(rows, A[cp1]));
                        assert!(rowIsAlive(rows, A[cp2]));
                    }
                    // Row indices will same order for both supercols,
                    // no gather scatter nessasary.
                    let _cp1 = cp1;
                    cp1 += 1;
                    let _cp2 = cp2;
                    cp2 += 1;
                    if A[_cp1] != A[_cp2] {
                        break;
                    }
                    i += 1;
                }

                // The two columns are different if the for-loop "broke".
                if i != length {
                    prevc = c;

                    c = cols[c as usize].hash_next();
                    continue;
                }

                // Got it! Two columns are identical.

                #[cfg(feature = "debug")]
                assert!(cols[c].score() == cols[superc].score());

                cols[superc as usize].set_thickness(
                    cols[superc as usize].thickness() + cols[c as usize].thickness(),
                );
                cols[c as usize].set_parent(superc);
                kill_non_principal_col(cols, c as usize);
                // Order c later, in orderChildren().
                cols[c as usize].set_order(empty);
                // Remove c from hash bucket.
                cols[prevc as usize].set_hash_next(cols[c as usize].hash_next());

                c = cols[c as usize].hash_next();
            }
            superc = cols[superc as usize].hash_next();
        }

        // Empty this hash bucket.
        if headColumn > empty {
            // Corresponding degree list "hash" is not empty.
            cols[headColumn as usize].set_headhash(empty);
        } else {
            // Corresponding degree list "hash" is empty.
            head[hash as usize] = empty;
        }
    }
}

// Defragments and compacts columns and rows in the workspace A. Used when
// all avaliable memory has been used while performing row merging. Returns
// the index of the first free position in A, after garbage collection. The
// time taken by this routine is linear in the size of the array A, which is
// itself linear in the number of nonzeros in the input matrix.
pub(crate) fn garbage_collection(
    n_row: i32,
    n_col: i32,
    rows: &mut [Row],
    cols: &mut [Col],
    A: &mut [i32],
    pfree: i32,
) -> i32 {
    // let psrc: i32; // Source pointer.
    // let pdest: i32; // Destination pointer.
    // let j: i32; // Counter
    // let r: usize; // A row index.
    // let c: usize; // A column index.
    // let length: i32; // Length of a row or column.

    let mut debugRows = 0;
    #[cfg(feature = "debug")]
    {
        debug2!("Defrag..");
        for psrc in 0..pfree {
            assert!(A[psrc] >= 0);
        }
        debugRows = 0
    }

    // Defragment the columns.

    let mut pdest = 0;
    for c in 0..n_col as usize {
        if col_is_alive(cols, c) {
            let mut psrc = cols[c].start;

            // Move and compact the column.
            #[cfg(feature = "debug")]
            assert!(pdest <= psrc);

            cols[c].start = pdest - 0;
            let length = cols[c].length;
            for j in 0..length {
                let r = A[psrc as usize] as usize;
                psrc += 1;
                if row_is_alive(rows, r) {
                    A[pdest as usize] = r as i32;
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
            rows[r].set_first_column(A[psrc as usize]);
            #[cfg(feature = "debug")]
            assert!(rowIsAlive(rows, r));

            // Flag the start of the row with the one's complement of row.
            A[psrc as usize] = ones_complement(r as i32);

            #[cfg(feature = "debug")]
            {
                debugRows += 1;
            }
        }
    }

    // Defragment the rows.

    let mut psrc = pdest;
    while psrc < pfree {
        // Find a negative number ... the start of a row.
        let _psrc = psrc;
        psrc += 1;
        if A[_psrc as usize] < 0 {
            psrc -= 1;
            // Get the row index.
            let r = ones_complement(A[psrc as usize]) as usize;
            #[cfg(feature = "debug")]
            assert!(r >= 0 && r < n_row);

            // Restore first column index.
            A[psrc as usize] = rows[r].first_column();
            #[cfg(feature = "debug")]
            {
                assert!(rowIsAlive(rows, r));
                assert!(rows[r].length > 0);
            }
            // Move and compact the row.
            #[cfg(feature = "debug")]
            assert!(pdest <= psrc);

            rows[r].start = pdest - 0;
            let length = rows[r].length;
            for j in 0..length {
                let c = A[psrc as usize] as usize;
                psrc += 1;
                if col_is_alive(cols, c) {
                    A[pdest as usize] = c as i32;
                    pdest += 1;
                }
            }
            rows[r].length = pdest - rows[r].start;
            #[cfg(feature = "debug")]
            {
                assert!(rows[r].length > 0);
                debugRows -= 1;
            }
        }
    }

    // Ensure we found all the rows.
    #[cfg(feature = "debug")]
    assert!(debugRows == 0);

    // Return the new value of pfree.

    pdest - 0
}

// Clears the Row[].shared2.mark array, and returns the new tagMark.
fn clear_mark(mut tagMark: i32, maxMark: i32, nrow: i32, row: &mut [Row]) -> i32 {
    if tagMark <= 0 || tagMark >= maxMark {
        for r in 0..nrow as usize {
            if row_is_alive(row, r) {
                row[r].set_mark(0)
            }
        }
        tagMark = 1;
    }
    tagMark
}

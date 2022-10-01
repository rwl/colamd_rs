#[cfg(feature = "debug")]
use crate::col::Col;
#[cfg(feature = "debug")]
use crate::internal::*;
#[cfg(feature = "debug")]
use crate::row::Row;

// At this point, all empty rows and columns are dead. All live columns
// are "clean" (containing no dead rows) and simplicial (no supercolumns
// yet). Rows may contain dead columns, but all live rows contain at
// least one live column.
#[cfg(feature = "debug")]
pub(crate) fn debug_structures(
    n_row: Int,
    n_col: Int,
    rows: &[Row],
    cols: &[Col],
    a_i: &[Int],
    n_col2: Int,
) {
    // Check A, Row, and Col.

    for c in 0..n_col as usize {
        if col_is_alive(cols, c) {
            let len = cols[c].length as usize;
            let score = cols[c].shared2.score();
            debug4!("initial live col {:5} {:5} {:5}\n", c, len, score);
            assert_debug!(len > 0);
            assert_debug!(score >= 0);
            assert_debug!(cols[c].shared1.thickness() == 1);
            let mut cp = cols[c].start as usize;
            let cp_end = cp + len;
            while cp < cp_end {
                let r = a_i[cp] as usize;
                cp += 1;
                assert_debug!(row_is_alive(rows, r));
            }
        } else {
            let i = cols[c].shared2.order();
            assert_debug!(i >= n_col2 && i < n_col);
        }
    }

    for r in 0..n_row as usize {
        if row_is_alive(rows, r) {
            let mut i = 0;
            let len = rows[r].length as usize;
            let deg = rows[r].shared1.degree();
            assert_debug!(len > 0);
            assert_debug!(deg > 0);
            let mut rp = rows[r].start as usize;
            let rp_end = rp + len;
            while rp < rp_end {
                let c = a_i[rp] as usize;
                rp += 1;
                if col_is_alive(cols, c) {
                    i += 1;
                }
            }
            assert_debug!(i > 0)
        }
    }
}

// Prints the contents of the degree lists. Counts the number of columns
// in the degree list and compares it to the total it should have. Also
// checks the row degrees.
#[cfg(feature = "debug")]
pub(crate) fn debug_deg_lists(
    n_row: Int,
    n_col: Int,
    rows: &[Row],
    cols: &[Col],
    head: &[Int],
    min_score: Int,
    should: Int,
    max_deg: Int,
) {
    // Check the degree lists.

    #[cfg(not(feature = "debug1"))]
    if n_col > 10000 && debugLevel <= 0 {
        return;
    }
    let mut have = 0;
    let mut d = format!("Degree lists: {}\n", min_score);
    for deg in 0..=n_col as usize {
        let mut col = head[deg];
        if col == EMPTY {
            continue;
        }
        d.push_str(&format!("{}:", deg));
        while col != EMPTY {
            d.push_str(&format!(" {}", col));
            have += cols[col as usize].shared1.thickness();
            assert_debug!(col_is_alive(cols, col as usize));
            col = cols[col as usize].shared4.degree_next();
        }
        debug4!("{}", d);
    }
    debug4!("should {} have {}", should, have);
    assert_debug!(should == have);

    // Check the row degrees.

    #[cfg(not(feature = "debug1"))]
    if n_row > 10000 && debugLevel <= 0 {
        return;
    }
    for row in 0..n_row as usize {
        if row_is_alive(rows, row) {
            assert_debug!(rows[row].shared1.degree() <= max_deg);
        }
    }
}

// Ensures that the tag_mark is less that the maximum and also ensures that
// each entry in the mark array is less than the tag mark.
#[cfg(feature = "debug")]
pub(crate) fn debug_mark(n_row: Int, rows: &[Row], tag_mark: Int, max_mark: Int) {
    // Check the Row marks.
    assert_debug!(tag_mark > 0 && tag_mark <= max_mark);
    #[cfg(not(feature = "debug1"))]
    if n_row > 10000 && debugLevel <= 0 {
        return;
    }
    for r in 0..n_row as usize {
        assert_debug!(rows[r].shared2.mark() < tag_mark)
    }
}

// Prints out the contents of the columns and the rows.
#[cfg(feature = "debug3")]
pub(crate) fn debug_matrix(n_row: Int, n_col: Int, rows: &[Row], cols: &[Col], a_i: &[Int]) {
    // Dump the rows and columns of the matrix.
    debug3!("DUMP MATRIX:");
    for r in 0..n_row as usize {
        debug3!("Row {} alive? {}", r, row_is_alive(rows, r));
        if row_is_dead(rows, r) {
            continue;
        }
        debug3!(
            "start {} length {} degree {}",
            rows[r].start,
            rows[r].length,
            rows[r].shared1.degree()
        );
        let mut rp = rows[r].start as usize;
        let rp_end = rp + rows[r].length as usize;
        while rp < rp_end {
            let c = a_i[rp] as usize;
            rp += 1;
            debug4!("	{} col {}\n", col_is_alive(cols, c), c);
            //if colIsAlive(cols, c) {
            //	debug4!("  1 col {}", c);
            //} else {
            //	debug4!("  0 col {}", c);
            //}
        }
    }

    for c in 0..n_col as usize {
        if col_is_alive(cols, c) {
            debug3!("Col {} alive? 1", c);
        } else {
            debug3!("Col {} alive? 0", c);
        }
        if col_is_dead(cols, c) {
            continue;
        }
        debug3!(
            "start {} length {} shared1 {} shared2 {}",
            cols[c].start,
            cols[c].length,
            cols[c].shared1.thickness(),
            cols[c].shared2.score()
        );
        let mut cp = cols[c].start as usize;
        let cp_end = cp + cols[c].length as usize;
        while cp < cp_end {
            let r = a_i[cp] as usize;
            cp += 1;
            if row_is_alive(rows, r) {
                debug4!("  1 row {}", r);
            } else {
                debug4!("  0 row {}", r);
            }
        }
    }
}

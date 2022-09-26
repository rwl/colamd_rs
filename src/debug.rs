use crate::col::Col;
use crate::internal::{col_is_alive, col_is_dead, empty, row_is_alive};
use crate::row::Row;

// At this point, all empty rows and columns are dead. All live columns
// are "clean" (containing no dead rows) and simplicial (no supercolumns
// yet). Rows may contain dead columns, but all live rows contain at
// least one live column.
fn debug_structures(n_row: i32, n_col: i32, rows: &[Row], cols: &[Col], A: &[i32], n_col2: i32) {
    #[cfg(feature = "debug")]
    {
        // Local variables.

        let (i, c, cp, cp_end, len, score, r, rp, rp_end, deg): (
            i32,
            i32,
            i32,
            i32,
            i32,
            i32,
            i32,
            i32,
            i32,
            i32,
        );

        // Check A, Row, and Col.

        for c in 0..n_col {
            if col_is_alive(cols, c) {
                len = cols[c].length;
                score = cols[c].score();
                debug4!("initial live col %5d %5d %5d\n", c, len, score);
                assert!(len > 0);
                assert!(score >= 0);
                assert!(cols[c].thickness() == 1);
                cp = cols[c].start;
                cp_end = cp + len;
                while cp < cp_end {
                    r = A[cp];
                    cp += 1;
                    assert!(rowIsAlive(rows, r));
                }
            } else {
                i = cols[c].order();
                assert!(i >= n_col2 && i < n_col);
            }
        }

        for r in 0..n_row {
            if rowIsAlive(rows, r) {
                i = 0;
                len = rows[r].length;
                deg = rows[r].degree();
                assert!(len > 0);
                assert!(deg > 0);
                rp = rows[r].start;
                rp_end = rp + len;
                while rp < rp_end {
                    c = A[rp];
                    rp += 1;
                    if colIsAlive(cols, c) {
                        i += 1;
                    }
                }
                assert!(i > 0)
            }
        }
    }
}

// Prints the contents of the degree lists. Counts the number of columns
// in the degree list and compares it to the total it should have. Also
// checks the row degrees.
fn debug_deg_lists(
    n_row: i32,
    n_col: i32,
    rows: &[Row],
    cols: &[Col],
    head: &[i32],
    min_score: i32,
    should: i32,
    max_deg: i32,
) {
    #[cfg(feature = "debug")]
    {
        // Local variables.

        let (deg, col, have, row): (i32, i32, i32, i32);

        // Check the degree lists.

        if n_col > 10000 && debugLevel <= 0 {
            return;
        }
        have = 0;
        let mut d = format!("Degree lists: {}\n", min_score);
        for deg in 0..=n_col {
            col = head[deg];
            if col == empty {
                continue;
            }
            d.push_str(format!("{}:", deg));
            while col != empty {
                d.push_str(format!(" {}", col));
                have += cols[col].thickness();
                assert!(colIsAlive(cols, col));
                col = cols[col].degreeNext();
            }
            debug4!(d);
        }
        debug4!("should {} have {}", should, have);
        assert!(should == have);

        // Check the row degrees.

        if n_row > 10000 && debugLevel <= 0 {
            return;
        }
        for row in 0..n_row {
            if row_is_alive(rows, row) {
                assert!(rows[row].degree() <= max_deg);
            }
        }
    }
}

// Ensures that the tag_mark is less that the maximum and also ensures that
// each entry in the mark array is less than the tag mark.
fn debug_mark(n_row: i32, rows: &[Row], tag_mark: i32, max_mark: i32) {
    #[cfg(feature = "debug")]
    {
        // Check the Row marks.
        assert!(tag_mark > 0 && tag_mark <= max_mark);
        if n_row > 10000 && debugLevel <= 0 {
            return;
        }
        for r in 0..n_row {
            assert!(rows[r].mark() < tag_mark)
        }
    }
}

// Prints out the contents of the columns and the rows.
fn debug_matrix(n_row: i32, n_col: i32, rows: &[Row], cols: &[Col], A: &[i32]) {
    #[cfg(feature = "debug3")]
    {
        // Dump the rows and columns of the matrix.
        // if debugLevel < 3 {
        // 	return
        // }
        debug3!("DUMP MATRIX:");
        for r in 0..n_row {
            debug3!("Row {} alive? {}", r, rowIsAlive(rows, r));
            if rowIsDead(rows, r) {
                continue;
            }
            debug3!(
                "start {} length {} degree {}",
                rows[r].start,
                rows[r].length,
                rows[r].degree()
            );
            let mut rp = rows[r].start;
            let rp_end = rp + rows[r].length;
            while rp < rp_end {
                let c = A[rp];
                rp += 1;
                debug4!("	{} col {}\n", colIsAlive(cols, c), c);
                //if colIsAlive(cols, c) {
                //	debug4!("  1 col {}", c);
                //} else {
                //	debug4!("  0 col {}", c);
                //}
            }
        }

        for c in 0..n_col {
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
                cols[c].thickness(),
                cols[c].score()
            );
            let mut cp = cols[c].start;
            let cp_end = cp + cols[c].length;
            while cp < cp_end {
                let r = A[cp];
                cp += 1;
                if row_is_alive(rows, r) {
                    debug4!("  1 row {}", r);
                } else {
                    debug4!("  0 row {}", r);
                }
            }
        }
    }
}

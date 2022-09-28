use crate::codes::*;
use crate::stats::*;

pub fn colamd_report(stats: &[i32; STATS]) {
    print_report("colamd", stats)
}

pub fn symamd_report(stats: &[i32; STATS]) {
    print_report("symamd", stats)
}

fn print_report(method: &str, stats: &[i32; STATS]) {
    #[cfg(feature = "nprint")]
    {
        return;
    }
    let s = format!("\n{}: ", method);

    // if stats == nil {
    // 	println!(s + "No statistics available.")
    // 	return
    // }

    let i1 = stats[INFO1].clone();
    let i2 = stats[INFO2].clone();
    let i3 = stats[INFO3].clone();

    if stats[STATUS] >= 0 {
        println!("{} OK.  ", s)
    } else {
        println!("{} ERROR.  ", s)
    }

    match stats[STATUS].clone() {
        OK_BUT_JUMBLED => {
            println!("Matrix has unsorted or duplicate row indices.");
            println!(
                "{}: number of duplicate or out-of-order row indices: {}",
                method, i3
            );
            println!(
                "{}: last seen duplicate or out-of-order row index:   {}",
                method,
                index(i2)
            );
            println!(
                "{}: last seen in column:                             {}",
                method,
                index(i1)
            );
            println!("");
            println!(
                "{}: number of dense or empty rows ignored:           {}",
                method, stats[DENSE_ROW]
            );
            println!(
                "{}: number of dense or empty columns ignored:        {}",
                method, stats[DENSE_COL]
            );
            println!(
                "{}: number of garbage collections performed:         {}",
                method, stats[DEFRAG_COUNT]
            );
        }
        OK => {
            println!();
            println!(
                "{}: number of dense or empty rows ignored:           {}",
                method, stats[DENSE_ROW]
            );
            println!(
                "{}: number of dense or empty columns ignored:        {}",
                method, stats[DENSE_COL]
            );
            println!(
                "{}: number of garbage collections performed:         {}",
                method, stats[DEFRAG_COUNT]
            );
        }
        // ERROR_A_NOT_PRESENT => {
        //     println!("Array A (row indices of matrix) not present.");
        // }
        // ERROR_P_NOT_PRESENT => {
        //     println!("Array p (column pointers for matrix) not present.");
        // }
        ERROR_NROW_NEGATIVE => {
            println!("Invalid number of rows ({}).", i1);
        }
        ERROR_NCOL_NEGATIVE => {
            println!("Invalid number of columns ({}).", i1);
        }
        ERROR_NNZ_NEGATIVE => {
            println!("Invalid number of nonzero entries ({}).", i1);
        }
        ERROR_P0_NONZERO => {
            println!("Invalid column pointer, p[0] = {}, must be zero.", i1);
        }
        ERROR_A_TOO_SMALL => {
            println!("Array A too small.");
            println!("        Need Alen >= {}, but given only Alen = {}.", i1, i2);
        }
        ERROR_COL_LENGTH_NEGATIVE => {
            println!(
                "Column {} has a negative number of nonzero entries ({}).",
                index(i1),
                i2
            );
        }
        ERROR_ROW_INDEX_OUT_OF_BOUNDS => {
            println!(
                "Row index (row {}) out of bounds ({} to {}) in column {}.",
                index(i2),
                index(0),
                index(i3 - 1),
                index(i1)
            );
        }
        ERROR_OUT_OF_MEMORY => {
            println!("Out of memory.");
        }
        _ => {}
    }
}

// Matrixes are 0-based and indices are reported as such.
fn index(i: i32) -> i32 {
    i
}

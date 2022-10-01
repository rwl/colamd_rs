use colamd::{colamd, colamd_report, symamd, symamd_report, Int, STATS};

const A_NNZ: Int = 11;
const A_NROW: Int = 5;
const A_NCOL: Int = 4;
const ALEN: Int = 150;

const B_NNZ: Int = 4;
const B_N: Int = 5;

/// COLAMD / SYMAMD example
///
/// `colamd` example of use, to order the columns of a 5-by-4 matrix with
/// 11 nonzero entries in the following nonzero pattern, with default knobs.
///
///     x 0 x 0
///     x 0 x x
///     0 x x 0
///     0 0 x x
///     x x 0 0
///
/// `symamd` example of use, to order the rows and columns of a 5-by-5
/// matrix with 13 nonzero entries in the following nonzero pattern,
/// with default knobs.
///
///     x x 0 0 0
///     x x x x 0
///     0 x x 0 0
///     0 x 0 x x
///     0 0 0 x x
///
/// (where x denotes a nonzero value).
fn main() -> Result<(), String> {
    let mut a_i = vec![0; ALEN as usize];
    let aa_i = vec![
        0, 1, 4, // Row indices of nonzeros in column 0.
        2, 4, // Row indices of nonzeros in column 1.
        0, 1, 2, 3, // Row indices of nonzeros in column 2.
        1, 3, // Row indices of nonzeros in column 3.
    ];
    for (i, a) in aa_i.into_iter().enumerate() {
        a_i[i] = a;
    }

    let mut p = vec![
        0,     // Column 0 is in A[0..2].
        3,     // Column 1 is in A[3..4].
        5,     // Column 2 is in A[5..8].
        9,     // Column 3 is in A[9..10].
        A_NNZ, // Number of nonzeros in A.
    ];

    // Input matrix B definition.

    // Note: Only strictly lower triangular part
    // is included, since symamd ignores the
    // diagonal and upper triangular part of B.
    let b_i = vec![
        1, // Row indices of nonzeros in column 0.
        2, 3, // Row indices of nonzeros in column 1.
        // Row indices of nonzeros in column 2 (none).
        4, // Row indices of nonzeros in column 3.
           // Row indices of nonzeros in column 4 (none).
    ];

    let q = vec![
        0,     // Column 0 is in B[0].
        1,     // Column 1 is in B[1..2].
        3,     // Column 2 is empty.
        3,     // Column 3 is in B[3].
        4,     // Column 4 is empty.
        B_NNZ, // Number of nonzeros in strictly lower B.
    ];

    // Other variable definitions.

    let mut perm = vec![0; (B_N + 1) as usize]; // Note the size is N+1.
    let mut stats = [0; STATS]; // For colamd and symamd output statistics.

    // Dump the input matrix A.

    println!("colamd {}-by-{} input matrix:", A_NROW, A_NCOL);
    for col in 0..A_NCOL as usize {
        let length = p[col + 1] - p[col];
        println!("Column {}, with {} entries:", col, length);
        for pp in p[col]..p[col + 1] {
            let row = a_i[pp as usize];
            println!("    row {}", row)
        }
    }

    // Order the matrix. Note that this destroys A and overwrites p.

    let ok = colamd(A_NROW, A_NCOL, ALEN, &mut a_i, &mut p, None, &mut stats);

    colamd_report(&stats);

    if !ok {
        return Err("colamd error!\n".to_string());
    }

    // Write the column ordering.

    println!("colamd column ordering:");
    println!("1st column: {}", p[0]);
    println!("2nd column: {}", p[1]);
    println!("3rd column: {}", p[2]);
    println!("4th column: {}", p[3]);

    if p[0] != 1 {
        return Err(format!("expected {} actual {}", 1, p[0]).to_string());
    }
    if p[1] != 0 {
        return Err(format!("expected {} actual {}", 0, p[1]).to_string());
    }
    if p[2] != 2 {
        return Err(format!("expected {} actual {}", 2, p[2]).to_string());
    }
    if p[3] != 3 {
        return Err(format!("expected {} actual {}", 3, p[3]).to_string());
    }

    // Dump the strictly lower triangular part of symmetric input matrix B.

    println!("\n\nsymamd {}-by-{} input matrix:", B_N, B_N);
    println!("Entries in strictly lower triangular part:");
    for col in 0..B_N as usize {
        let length = q[col + 1] - q[col];
        println!("Column {}, with {} entries:", col, length);
        for pp in q[col]..q[col + 1] {
            let row = b_i[pp as usize];
            println!("    row {}", row);
        }
    }

    // Order the matrix B.  Note that this does not modify B or q.

    let ok = symamd(B_N, &b_i, &q, &mut perm, None, &mut stats);
    symamd_report(&stats);

    if !ok {
        return Err("symamd error!\n".to_string());
    }

    // Write the symmetric ordering.

    println!("symamd column ordering:");
    println!("1st row/column: {}", perm[0]);
    println!("2nd row/column: {}", perm[1]);
    println!("3rd row/column: {}", perm[2]);
    println!("4th row/column: {}", perm[3]);
    println!("5th row/column: {}", perm[4]);

    if perm[0] != 0 {
        return Err(format!("expected {} actual {}", 0, perm[0]).to_string());
    }
    if perm[1] != 2 {
        return Err(format!("expected {} actual {}", 2, perm[1]).to_string());
    }
    if perm[2] != 1 {
        return Err(format!("expected {} actual {}", 1, perm[2]).to_string());
    }
    if perm[3] != 3 {
        return Err(format!("expected {} actual {}", 3, perm[3]).to_string());
    }
    if perm[4] != 4 {
        return Err(format!("expected {} actual {}", 4, perm[4]).to_string());
    }
    Ok(())
}

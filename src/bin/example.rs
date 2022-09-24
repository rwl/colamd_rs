use colamd::colamd;

const A_NNZ: i32 = 11;
const A_NROW: i32 = 5;
const A_NCOL: i32 = 4;
const ALEN: i32 = 150;

const B_NNZ: i32 = 4;
const B_N: i32 = 5;

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
fn main() {
    let AA = vec![
        0, 1, 4, // Row indices of nonzeros in column 0.
        2, 4, // Row indices of nonzeros in column 1.
        0, 1, 2, 3, // Row indices of nonzeros in column 2.
        1, 3, // Row indices of nonzeros in column 3.
    ];
    let A = AA.clone();

    let p = vec![
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
    let B = vec![
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

    let perm = vec![0; (B_N + 1) as usize]; // Note the size is N+1.
    let stats = vec![0; STATS]; // For colamd and symamd output statistics.

    // Dump the input matrix A.

    println!("colamd {}-by-{} input matrix:", A_NROW, A_NCOL);
    for col in 0..A_NCOL {
        let length = p[col + 1] - p[col];
        println!("Column {}, with {} entries:", col, length);
        for pp in p[col]..p[col + 1] {
            row = A[pp];
            println!("    row {}", row)
        }
    }

    // Order the matrix.  Note that this destroys A and overwrites p.

    let ok = colamd(A_NROW, A_NCOL, ALEN, A, p, nil, stats);
}

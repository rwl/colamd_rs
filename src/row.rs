#[derive(Clone, Default)]
pub(crate) struct Row {
    // Index for A of first col in this row.
    pub(crate) start: i32,
    // Number of principal columns in this row.
    pub(crate) length: i32,

    shared1: i32,
    shared2: i32,
}

// pub(crate) enum RowShared1 {
//     Degree(i32),
//     P(i32),
// }
//
// impl RowShared1 {
//     // Number of principal & non-principal columns in row.
//     pub(crate) fn degree(self) -> i32 {
//         match self {
//             RowShared1::Degree(val) => val,
//             RowShared1::P(val) => panic!("called `RowShared1::degree()` on a `P` value: {}", val),
//         }
//     }
//
//     // Used as a row pointer in `init_rows_cols()`.
//     pub(crate) fn p(self) -> i32 {
//         match self {
//             RowShared1::P(val) => val,
//             RowShared1::Degree(val) => {
//                 panic!("called `RowShared1::p()` on a `Degree` value: {}", val)
//             }
//         }
//     }
// }

impl Row {
    // Number of principal & non-principal columns in row.
    pub(crate) fn degree(&self) -> i32 {
        self.shared1
    }

    pub(crate) fn set_degree(&mut self, degree: i32) {
        self.shared1 = degree
    }

    // Used as a row pointer in initRowsCols().
    pub(crate) fn p(&self) -> i32 {
        self.shared1
    }

    pub(crate) fn set_p(&mut self, p: i32) {
        self.shared1 = p
    }

    // For computing set differences and marking dead rows.
    pub(crate) fn mark(&self) -> i32 {
        self.shared2
    }

    pub(crate) fn set_mark(&mut self, mark: i32) {
        self.shared2 = mark
    }

    // First column in row (used in garbage collection).
    pub(crate) fn first_column(&self) -> i32 {
        self.shared2
    }

    pub(crate) fn set_first_column(&mut self, first_column: i32) {
        self.shared2 = first_column
    }
}

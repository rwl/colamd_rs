#[derive(Clone, Default)]
pub(crate) struct Row {
    // Index for A of first col in this row.
    pub(crate) start: i32,
    // Number of principal columns in this row.
    pub(crate) length: i32,

    pub(crate) shared1: RowShared1,
    pub(crate) shared2: RowShared2,
}

#[derive(Clone)]
pub(crate) enum RowShared1 {
    Degree(i32),
    P(i32),
}

impl RowShared1 {
    // Number of principal & non-principal columns in row.
    pub(crate) fn degree(&self) -> i32 {
        match self {
            RowShared1::Degree(val) => *val,
            RowShared1::P(val) => panic!("called `RowShared1::degree()` on a `P` value: {}", val),
        }
    }

    // Used as a row pointer in `init_rows_cols()`.
    pub(crate) fn p(&self) -> i32 {
        match self {
            RowShared1::P(val) => *val,
            RowShared1::Degree(val) => {
                panic!("called `RowShared1::p()` on a `Degree` value: {}", val)
            }
        }
    }
}

impl Default for RowShared1 {
    fn default() -> Self {
        RowShared1::Degree(0)
    }
}

#[derive(Clone)]
pub(crate) enum RowShared2 {
    Mark(i32),
    FirstColumn(i32),
}

impl RowShared2 {
    // For computing set differences and marking dead rows.
    pub(crate) fn mark(&self) -> i32 {
        match self {
            RowShared2::Mark(val) => *val,
            RowShared2::FirstColumn(val) => panic!(
                "called `RowShared2::mark()` on a `FirstColumn` value: {}",
                val
            ),
        }
    }

    // First column in row (used in garbage collection).
    pub(crate) fn first_column(&self) -> i32 {
        match self {
            RowShared2::FirstColumn(val) => *val,
            RowShared2::Mark(val) => {
                panic!(
                    "called `RowShared1::first_column()` on a `mark` value: {}",
                    val
                )
            }
        }
    }
}

impl Default for RowShared2 {
    fn default() -> Self {
        RowShared2::Mark(-1)
    }
}

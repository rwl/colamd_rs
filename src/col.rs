#[derive(Clone, Default)]
pub(crate) struct Col {
    /// Index for A of first row in this column, or DEAD if column is dead.
    pub(crate) start: i32,
    /// Number of rows in this column.
    pub(crate) length: i32,

    shared1: i32,
    shared2: i32,
    shared3: i32,
    shared4: i32,
}

pub(crate) enum ColShared1 {
    Thickness(i32),
    Parent(i32),
}

impl ColShared1 {
    // Number of original columns represented by this col, if the column is alive.
    pub(crate) fn thickness(self) -> i32 {
        match self {
            ColShared1::Thickness(val) => val,
            ColShared1::Parent(val) => panic!(
                "called `ColShared1::thickness()` on a `Parent` value: {}",
                val
            ),
        }
    }

    // Parent in parent tree super-column structure, if the column is dead.
    pub(crate) fn parent(self) -> i32 {
        match self {
            ColShared1::Parent(val) => val,
            ColShared1::Thickness(val) => {
                panic!(
                    "called `ColShared1::parent()` on a `Thickness` value: {}",
                    val
                )
            }
        }
    }
}

impl Col {
    // Number of original columns represented by this
    // col, if the column is alive.
    pub(crate) fn thickness(&self) -> i32 {
        self.shared1
    }

    pub(crate) fn set_thickness(&mut self, thickness: i32) {
        self.shared1 = thickness
    }

    // Parent in parent tree super-column structure, if the column is dead.
    pub(crate) fn parent(&self) -> i32 {
        self.shared1
    }

    pub(crate) fn set_parent(&mut self, parent: i32) {
        self.shared1 = parent
    }

    // The score used to maintain heap, if col is alive.
    pub(crate) fn score(&self) -> i32 {
        self.shared2
    }

    pub(crate) fn set_score(&mut self, score: i32) {
        self.shared2 = score
    }

    // Pivot ordering of this column, if col is dead.
    pub(crate) fn order(&self) -> i32 {
        self.shared2
    }

    pub(crate) fn set_order(&mut self, order: i32) {
        self.shared2 = order
    }

    // Head of a hash bucket, if col is at the head of a degree list.
    pub(crate) fn headhash(&self) -> i32 {
        self.shared3
    }

    pub(crate) fn set_headhash(&mut self, headhash: i32) {
        self.shared3 = headhash
    }

    // Hash value, if col is not in a degree list.
    pub(crate) fn hash(&self) -> i32 {
        self.shared3
    }

    pub(crate) fn set_hash(&mut self, hash: i32) {
        self.shared3 = hash
    }

    // Previous column in degree list, if col is in a
    // degree list (but not at the head of a degree list).
    pub(crate) fn prev(&self) -> i32 {
        self.shared3
    }

    pub(crate) fn set_prev(&mut self, prev: i32) {
        self.shared3 = prev
    }

    // Next column, if col is in a degree list.
    pub(crate) fn degree_next(&self) -> i32 {
        self.shared4
    }

    pub(crate) fn set_degree_next(&mut self, degreeNext: i32) {
        self.shared4 = degreeNext
    }

    // Next column, if col is in a hash list.
    pub(crate) fn hash_next(&self) -> i32 {
        self.shared4
    }

    pub(crate) fn set_hash_next(&mut self, hashNext: i32) {
        self.shared4 = hashNext
    }
}

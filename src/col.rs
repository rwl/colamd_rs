use crate::internal::{Int, EMPTY};

#[derive(Clone, Default)]
pub(crate) struct Col {
    /// Index for A of first row in this column, or DEAD if column is dead.
    pub(crate) start: Int,
    /// Number of rows in this column.
    pub(crate) length: Int,

    pub(crate) shared1: ColShared1,
    pub(crate) shared2: ColShared2,
    pub(crate) shared3: ColShared3,
    pub(crate) shared4: ColShared4,
}

#[derive(Clone)]
pub(crate) enum ColShared1 {
    Thickness(Int),
    Parent(Int),
}

impl ColShared1 {
    // Number of original columns represented by this col, if the column is alive.
    pub(crate) fn thickness(&self) -> Int {
        match self {
            ColShared1::Thickness(val) => *val,
            ColShared1::Parent(val) => panic!(
                "called `ColShared1::thickness()` on a `Parent` value: {}",
                val
            ),
        }
    }

    // Parent in parent tree super-column structure, if the column is dead.
    pub(crate) fn parent(&self) -> Int {
        match self {
            ColShared1::Parent(val) => *val,
            ColShared1::Thickness(val) => {
                panic!(
                    "called `ColShared1::parent()` on a `Thickness` value: {}",
                    val
                )
            }
        }
    }
}

impl Default for ColShared1 {
    fn default() -> Self {
        ColShared1::Thickness(1)
    }
}

#[derive(Clone)]
pub(crate) enum ColShared2 {
    Score(Int),
    Order(Int),
}

impl ColShared2 {
    // The score used to maintain heap, if col is alive.
    pub(crate) fn score(&self) -> Int {
        match self {
            ColShared2::Score(val) => *val,
            ColShared2::Order(val) => {
                panic!("called `ColShared2::score()` on a `Order` value: {}", val)
            }
        }
    }

    // Pivot ordering of this column, if col is dead.
    pub(crate) fn order(&self) -> Int {
        match self {
            ColShared2::Order(val) => *val,
            ColShared2::Score(val) => {
                panic!("called `ColShared2::order()` on a `Score` value: {}", val)
            }
        }
    }
}

impl Default for ColShared2 {
    fn default() -> Self {
        ColShared2::Score(0)
    }
}

#[derive(Clone)]
pub(crate) enum ColShared3 {
    HeadHash(Int),
    Hash(Int),
    Prev(Int),
}

impl ColShared3 {
    // Head of a hash bucket, if col is at the head of a degree list.
    pub(crate) fn head_hash(&self) -> Int {
        match self {
            ColShared3::HeadHash(val) => *val,
            ColShared3::Hash(val) => {
                panic!(
                    "called `ColShared3::head_hash()` on a `Hash` value: {}",
                    val
                )
            }
            ColShared3::Prev(val) => {
                panic!(
                    "called `ColShared3::head_hash()` on a `Prev` value: {}",
                    val
                )
            }
        }
    }

    // Hash value, if col is not in a degree list.
    pub(crate) fn hash(&self) -> Int {
        match self {
            ColShared3::Hash(val) => *val,
            ColShared3::HeadHash(val) => {
                panic!("called `ColShared3::hash()` on a `HeadHash` value: {}", val)
            }
            ColShared3::Prev(val) => {
                panic!("called `ColShared3::hash()` on a `Prev` value: {}", val)
            }
        }
    }

    // Previous column in degree list, if col is in a
    // degree list (but not at the head of a degree list).
    pub(crate) fn prev(&self) -> Int {
        match self {
            ColShared3::Prev(val) => *val,
            ColShared3::HeadHash(val) => {
                panic!("called `ColShared3::prev()` on a `HeadHash` value: {}", val)
            }
            ColShared3::Hash(val) => {
                panic!("called `ColShared3::prev()` on a `Hash` value: {}", val)
            }
        }
    }

    pub(crate) fn unwrap(&self) -> Int {
        match self {
            ColShared3::Prev(val) => *val,
            ColShared3::HeadHash(val) => *val,
            ColShared3::Hash(val) => *val,
        }
    }
}

impl Default for ColShared3 {
    fn default() -> Self {
        ColShared3::Prev(EMPTY)
    }
}

#[derive(Clone)]
pub(crate) enum ColShared4 {
    DegreeNext(Int),
    HashNext(Int),
}

impl ColShared4 {
    // Next column, if col is in a degree list.
    pub(crate) fn degree_next(&self) -> Int {
        match self {
            ColShared4::DegreeNext(val) => *val,
            ColShared4::HashNext(val) => {
                panic!(
                    "called `ColShared4::degree_next()` on a `HashNext` value: {}",
                    val
                )
            }
        }
    }

    // Next column, if col is in a hash list.
    pub(crate) fn hash_next(&self) -> Int {
        match self {
            ColShared4::HashNext(val) => *val,
            ColShared4::DegreeNext(val) => {
                panic!(
                    "called `ColShared4::hash_next()` on a `DegreeNext` value: {}",
                    val
                )
            }
        }
    }
}

impl Default for ColShared4 {
    fn default() -> Self {
        ColShared4::DegreeNext(EMPTY)
    }
}

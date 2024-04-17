#[derive(Debug, Clone)]
pub struct Plaintext {
    pub(crate) msg: Vec<u64>,
    pub(crate) slots: usize,
    pub(crate) l: usize,
}

#[derive(Debug, Clone)]
pub struct Ciphertext {
    /// (ax, bx) satisfies ax * s + bx \approx m
    pub(crate) ax: Vec<u64>,
    /// bx is of the form: ex + m - ax * s
    pub(crate) bx: Vec<u64>,
    /// dimension of ring
    pub(crate) n: usize,
    /// the length of plaintext vector
    pub(crate) slots: usize,
    /// the level of the ciphertext
    pub(crate) l: usize,
}

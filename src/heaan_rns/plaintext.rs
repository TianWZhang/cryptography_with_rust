#[derive(Debug, Clone)]
pub struct Plaintext {
    pub(crate) msg: Vec<u64>,
    pub(crate) slots: usize,
    pub(crate) l: usize,
}

#[derive(Debug, Clone)]
pub struct Ciphertext {
    pub(crate) ax: Vec<u64>,
    pub(crate) bx: Vec<u64>,
    pub(crate) n: usize,     // dimension of ring
    pub(crate) slots: usize, // the length of plaintext vector
    pub(crate) l: usize,     // the level of the ciphertext
}

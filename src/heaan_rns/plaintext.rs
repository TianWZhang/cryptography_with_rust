#[derive(Debug, Clone)]
pub struct Plaintext {
    msg: Vec<u64>,
    slots: usize,
    l: usize,
}

#[derive(Debug, Clone)]
pub struct Ciphertext {
    ax: Vec<u64>,
    bx: Vec<u64>,
    n: usize,     // dimension of ring
    slots: usize, // the length of plaintext vector
    l: usize,     // the level of the ciphertext
}

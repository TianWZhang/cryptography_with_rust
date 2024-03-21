use super::context::Context;

#[derive(Hash, PartialEq, Eq)]
pub enum KeyType {
    Encryption,
    Multiplication,
    Conjugation,
}

#[derive(Debug, Clone)]
pub struct Sk {
    pub(crate) sx: Vec<u64>,
}

impl Sk {
    pub fn new(context: &Context) -> Self {
        let mut sx = context.sample_hwt(context.max_level, context.num_special_modulus);
        context.ntt(&mut sx, context.max_level, context.num_special_modulus);
        Self { sx }
    }
}

pub struct Key {
    pub(crate) ax: Vec<u64>,
    pub(crate) bx: Vec<u64>,
}

impl Key {
    pub fn new(ax: Vec<u64>, bx: Vec<u64>) -> Self {
        Self { ax, bx }
    }
}

use rand::prelude::*;

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct ByteArray {
    data: Vec<u8>,
}

impl ByteArray {
    pub const fn new() -> Self {
        Self { data: Vec::new() }
    }

    pub fn from_bytes(data: &[u8]) -> Self {
        Self {
            data: data.to_vec(),
        }
    }

    pub fn random_bytes(len: usize) -> Self {
        let mut data = vec![0; len];
        let mut rng = rand::thread_rng();
        rng.fill_bytes(&mut data);

        Self { data }
    }

    pub fn append(&mut self, other: &Self) {
        self.data.extend_from_slice(&other.data);
    }

    pub fn concat(items: &[&Self]) -> Self {
        let len = items.iter().map(|s| s.data.len()).sum();
        let mut data = Vec::with_capacity(len);
        for item in items.iter() {
            data.extend_from_slice(&item.data)
        }

        Self { data }
    }

    pub fn get_bit(&self, pos: usize) -> bool {
        let (index, offset) = (pos / 8, pos % 8);
        let mask = 1 << offset;
        (self.data[index] & mask) != 0
    }

    pub fn skip(&self, num: usize) -> Self {
        let data = if num < self.data.len() {
            Vec::from(&self.data[..num])
        } else {
            Vec::new()
        };
        Self { data }
    }

    pub fn split_at(&self, pos: usize) -> (Self, Self) {
        let (d1, d2) = self.data.split_at(pos);
        (Self { data: d1.to_vec() }, Self { data: d2.to_vec() })
    }

    pub fn truncate(&self, len: usize) -> Self {
        let mut data = self.data.clone();
        data.truncate(len);
        Self { data }
    }
}

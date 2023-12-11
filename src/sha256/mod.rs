const H: [u32; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

pub const BLOCK_SIZE: usize = 64;

pub struct Sha256 {
    state: [u32; 8],
    raw_bytes: Vec<u8>,
    total_bits: usize,
    finalized: bool,
}

impl Sha256 {
    pub fn new() -> Self {
        Self {
            state: H,
            raw_bytes: Vec::new(),
            total_bits: 0,
            finalized: false,
        }
    }

    /// update with a block of 64 bytes
    fn update_block(&mut self, block: &[u8; BLOCK_SIZE]) {
        let mut w = [0u32; BLOCK_SIZE];
        for i in 0..16 {
            w[i] = bytes_to_u32(&block[4 * i..4 * (i + 1)]);
        }
        for i in 16..BLOCK_SIZE {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let mut a = self.state[0];
        let mut b = self.state[1];
        let mut c = self.state[2];
        let mut d = self.state[3];
        let mut e = self.state[4];
        let mut f = self.state[5];
        let mut g = self.state[6];
        let mut h = self.state[7];
        for i in 0..BLOCK_SIZE {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ (!e & g);
            let temp1 = h
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(maj);

            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }

        self.state[0] = self.state[0].wrapping_add(a);
        self.state[1] = self.state[1].wrapping_add(b);
        self.state[2] = self.state[2].wrapping_add(c);
        self.state[3] = self.state[3].wrapping_add(d);
        self.state[4] = self.state[4].wrapping_add(e);
        self.state[5] = self.state[5].wrapping_add(f);
        self.state[6] = self.state[6].wrapping_add(g);
        self.state[7] = self.state[7].wrapping_add(h);
    }

    /// consume blocks of raw bytes
    fn consume(&mut self) {
        assert_eq!(self.raw_bytes.len() % BLOCK_SIZE, 0);
        let n_blocks = self.raw_bytes.len() / BLOCK_SIZE;
        for i in 0..n_blocks {
            let mut block = [0u8; BLOCK_SIZE];
            block.copy_from_slice(&self.raw_bytes[i * BLOCK_SIZE..(i + 1) * BLOCK_SIZE]);
            self.update_block(&block);
        }
        self.raw_bytes.clear();
    }

    pub fn digest(&mut self) -> [u32; 8] {
        self.finalize();
        self.state
    }

    fn padding(&mut self) {
        self.raw_bytes.push(0x80);
        let n_padding_bytes = (512 - (self.total_bits + 8 + 64) % 512) / 8;
        self.raw_bytes.extend(vec![0; n_padding_bytes]);
        self.raw_bytes
            .extend((self.total_bits as u64).to_be_bytes());
    }

    /// update the hash state
    pub fn update(&mut self, input: &[u8]) {
        if self.finalized || input.len() == 0 {
            return;
        }

        self.total_bits += input.len() * 8;
        let len = input.len() + self.raw_bytes.len();
        if len >= BLOCK_SIZE {
            self.raw_bytes
                .extend(&input[..(input.len() - len % BLOCK_SIZE)]);
            self.consume();
        }
        if len % BLOCK_SIZE > 0 {
            self.raw_bytes
                .extend(&input[(input.len() - len % BLOCK_SIZE)..]);
        }
    }

    fn finalize(&mut self) {
        self.padding();
        self.consume();
        self.finalized = true;
    }

    pub fn hex_digest(&mut self) -> String {
        let digest = self.digest();
        format!(
            "{:08X?}{:08X?}{:08X?}{:08X?}{:08X?}{:08X?}{:08X?}{:08X?}",
            digest[0], digest[1], digest[2], digest[3], digest[4], digest[5], digest[6], digest[7]
        )
    }
}

/// bytes is in big endian order
fn bytes_to_u32(bytes: &[u8]) -> u32 {
    assert!(bytes.len() == 4);
    let mut res = 0;
    for i in 0..4 {
        res |= (bytes[i] as u32) << (8 * (3 - i));
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Run empty test input from FIPS 180-2
    #[test]
    fn sha256_nist_empty() {
        let mut hasher = Sha256::new();
        hasher.update(&[]);
        let digest = hasher.digest();
        assert_eq!(
            digest,
            [
                0xe3b0c442, 0x98fc1c14, 0x9afbf4c8, 0x996fb924, 0x27ae41e4, 0x649b934c, 0xa495991b,
                0x7852b855
            ]
        );
    }

    /// Run abc test from FIPS 180-2
    #[test]
    fn sha256_nist_abc() {
        let mut hasher = Sha256::new();
        hasher.update(b"abc");
        let digest = hasher.digest();
        assert_eq!(
            digest,
            [
                0xba7816bf, 0x8f01cfea, 0x414140de, 0x5dae2223, 0xb00361a3, 0x96177a9c, 0xb410ff61,
                0xf20015ad
            ]
        )
    }

    /// Run two-block test from FIPS 180-2
    #[test]
    fn sha256_nist_two_blocks() {
        let mut hasher = Sha256::new();
        hasher.update(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq");
        let digest = hasher.digest();
        assert_eq!(
            digest,
            [
                0x248d6a61, 0xd20638b8, 0xe5c02693, 0x0c3e6039, 0xa33ce459, 0x64ff2167, 0xf6ecedd4,
                0x19db06c1
            ]
        );
    }

    /// Run large input test (1,000,000 x a) from FIPS 180-2
    #[test]
    fn sha256_nist_large_input() {
        let input_str = std::iter::repeat("a").take(1_000_000).collect::<String>();
        let input = input_str.as_bytes();
        let mut hasher = Sha256::new();
        hasher.update(&input);
        let digest = hasher.digest();
        assert_eq!(
            digest,
            [
                0xcdc76e5c, 0x9914fb92, 0x81a1c7e2, 0x84d73e67, 0xf1809a48, 0xa497200e, 0x046d39cc,
                0xc7112cd0
            ]
        );
    }
}

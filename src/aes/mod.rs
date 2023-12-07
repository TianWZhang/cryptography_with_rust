mod aes_ecb;
mod aes_gcm;
mod constants;

use constants::*;
use std::vec;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum AESKeyType {
    AES128,
    AES192,
    AES256,
}

pub struct AESKey {
    key: Vec<u8>,
    key_type: AESKeyType,
}

impl AESKey {
    pub fn new(key: &[u8]) -> AESKey {
        match key.len() {
            16 => AESKey {
                key: key.to_vec(),
                key_type: AESKeyType::AES128,
            },
            24 => AESKey {
                key: key.to_vec(),
                key_type: AESKeyType::AES192,
            },
            32 => AESKey {
                key: key.to_vec(),
                key_type: AESKeyType::AES256,
            },
            _ => panic!("Invalid key length"),
        }
    }
}

pub struct AESScheme {
    key_type: AESKeyType,
    /// AES128: 11 subkeys, each of 128 bits
    /// AES192: 13 subkeys, each of 128 bits
    /// AES256: 15 subkeys, each of 128 bits
    round_keys: Vec<u32>,
}

fn get_u32(bytes: &[u8]) -> u32 {
    ((bytes[0] as u32) << 24)
        | ((bytes[1] as u32) << 16)
        | ((bytes[2] as u32) << 8)
        | (bytes[3] as u32)
}

impl AESScheme {
    pub fn new(key: &[u8]) -> Self {
        let key = AESKey::new(key);
        AESScheme {
            key_type: key.key_type,
            round_keys: key_schedule_aes(&key),
        }
    }

    pub fn encrypt_block(&self, plaintext: &[u8; 16]) -> [u8; 16] {
        let nb = 4;
        let nr = match self.key_type {
            AESKeyType::AES128 => 10,
            AESKeyType::AES192 => 12,
            AESKeyType::AES256 => 14,
        };
        let mut state = plaintext.clone();
        let w = &self.round_keys;
        add_round_key(&mut state, &w[0..nb]);
        for round in 1..nr {
            bytes_substitution(&mut state);
            shift_rows(&mut state);
            mix_column(&mut state);
            add_round_key(&mut state, &w[(round * nb)..((round + 1) * nb)]);
        }

        bytes_substitution(&mut state);
        shift_rows(&mut state);
        add_round_key(&mut state, &w[(nr * nb)..((nr + 1) * nb)]);

        state
    }

    pub fn encrypt_block_t_table(&self, plaintext: &[u8; 16]) -> [u8; 16] {
        let nr = match self.key_type {
            AESKeyType::AES128 => 10,
            AESKeyType::AES192 => 12,
            AESKeyType::AES256 => 14,
        };
        let w = &self.round_keys;
        let mut state = [
            get_u32(&plaintext[0..4]) ^ w[0],
            get_u32(&plaintext[4..8]) ^ w[1],
            get_u32(&plaintext[8..12]) ^ w[2],
            get_u32(&plaintext[12..16]) ^ w[3],
        ];

        for i in 1..nr {
            let temp = state.clone();
            state[0] = TE0[(temp[0] >> 24) as usize]
                ^ TE1[((temp[1] >> 16) as usize) & 0xff]
                ^ TE2[((temp[2] >> 8) as usize) & 0xff]
                ^ TE3[(temp[3] as usize) & 0xff]
                ^ w[4 * i];
            state[1] = TE0[(temp[1] >> 24) as usize]
                ^ TE1[((temp[2] >> 16) as usize) & 0xff]
                ^ TE2[((temp[3] >> 8) as usize) & 0xff]
                ^ TE3[(temp[0] as usize) & 0xff]
                ^ w[4 * i + 1];
            state[2] = TE0[(temp[2] >> 24) as usize]
                ^ TE1[((temp[3] >> 16) as usize) & 0xff]
                ^ TE2[((temp[0] >> 8) as usize) & 0xff]
                ^ TE3[(temp[1] as usize) & 0xff]
                ^ w[4 * i + 2];
            state[3] = TE0[(temp[3] >> 24) as usize]
                ^ TE1[((temp[0] >> 16) as usize) & 0xff]
                ^ TE2[((temp[1] >> 8) as usize) & 0xff]
                ^ TE3[(temp[2] as usize) & 0xff]
                ^ w[4 * i + 3];
        }

        let temp = state.clone();
        state[0] = ((AES_SBOX[(temp[0] >> 24) as usize] as u32) << 24
            | (AES_SBOX[((temp[1] >> 16) & 0xff) as usize] as u32) << 16
            | (AES_SBOX[((temp[2] >> 8) & 0xff) as usize] as u32) << 8
            | AES_SBOX[(temp[3] & 0xff) as usize] as u32)
            ^ w[4 * nr];
        state[1] = ((AES_SBOX[(temp[1] >> 24) as usize] as u32) << 24
            | (AES_SBOX[((temp[2] >> 16) & 0xff) as usize] as u32) << 16
            | (AES_SBOX[((temp[3] >> 8) & 0xff) as usize] as u32) << 8
            | AES_SBOX[(temp[0] & 0xff) as usize] as u32)
            ^ w[4 * nr + 1];
        state[2] = ((AES_SBOX[(temp[2] >> 24) as usize] as u32) << 24
            | (AES_SBOX[((temp[3] >> 16) & 0xff) as usize] as u32) << 16
            | (AES_SBOX[((temp[0] >> 8) & 0xff) as usize] as u32) << 8
            | AES_SBOX[(temp[1] & 0xff) as usize] as u32)
            ^ w[4 * nr + 2];
        state[3] = ((AES_SBOX[(temp[3] >> 24) as usize] as u32) << 24
            | (AES_SBOX[((temp[0] >> 16) & 0xff) as usize] as u32) << 16
            | (AES_SBOX[((temp[1] >> 8) & 0xff) as usize] as u32) << 8
            | AES_SBOX[(temp[2] & 0xff) as usize] as u32)
            ^ w[4 * nr + 3];

        let mut res = [0; 16];
        for i in 0..4 {
            res[4 * i] = (state[i] >> 24) as u8;
            res[4 * i + 1] = (state[i] >> 16) as u8;
            res[4 * i + 2] = (state[i] >> 8) as u8;
            res[4 * i + 3] = state[i] as u8;
        }
        res
    }

    pub fn decrypt_block(&self, ciphertext: &[u8; 16]) -> [u8; 16] {
        let nb = 4;
        let nr = match self.key_type {
            AESKeyType::AES128 => 10,
            AESKeyType::AES192 => 12,
            AESKeyType::AES256 => 14,
        };
        let mut state = ciphertext.clone();
        let w = &self.round_keys;

        add_round_key(&mut state, &w[(nr * nb)..((nr + 1) * nb)]);
        inv_shift_rows(&mut state);
        inv_byte_substitution(&mut state);

        for round in (1..nr).rev() {
            add_round_key(&mut state, &w[(round * nb)..((round + 1) * nb)]);
            inv_mix_column(&mut state);
            inv_shift_rows(&mut state);
            inv_byte_substitution(&mut state);
        }

        add_round_key(&mut state, &w[0..nb]);

        state
    }

    pub fn cbc_encrypt(&self, iv: &mut [u8], data: &[u8]) -> Vec<u8> {
        let mut padded = pad(data);
        for block in padded.chunks_mut(AES_BLOCK_SIZE) {
            let mut plaintext = [0u8; AES_BLOCK_SIZE];
            for i in 0..AES_BLOCK_SIZE {
                plaintext[i] = block[i] ^ iv[i];
            }
            let ciphertext = self.encrypt_block(&plaintext);
            for i in 0..AES_BLOCK_SIZE {
                block[i] = ciphertext[i];
                iv[i] = ciphertext[i];
            }
        }
        padded
    }

    pub fn cbc_decrypt(&self, iv: &mut [u8], data: &[u8]) -> Vec<u8> {
        assert!(data.len() % AES_BLOCK_SIZE == 0);
        let mut res = data.to_vec();

        for block in res.chunks_mut(AES_BLOCK_SIZE) {
            let mut ciphertext = [0u8; AES_BLOCK_SIZE];
            for i in 0..AES_BLOCK_SIZE {
                ciphertext[i] = block[i];
            }
            let plaintext = self.decrypt_block(&ciphertext);
            for i in 0..AES_BLOCK_SIZE {
                block[i] = plaintext[i] ^ iv[i];
                iv[i] = ciphertext[i];
            }
        }
        unpad(&mut res);
        res
    }
}

/// XOR the round key with the state and put the result in state
// State:
// +----+----+-----+-----+
// | b0 | b4 | b8  | b12 |
// +----+----+-----+-----+
// | b1 | b5 | b9  | b13 |
// +----+----+-----+-----+
// | b2 | b6 | b10 | b14 |
// +----+----+-----+-----+
// | b3 | b7 | b11 | b15 |
// +----+----+-----+-----+
// [b0, b1, b2, b3] = word[0]
fn add_round_key(state: &mut [u8; 16], round_key: &[u32]) {
    for i in 0..4 {
        state[4 * i] ^= (round_key[i] >> 24) as u8;
        state[4 * i + 1] ^= (round_key[i] >> 16) as u8;
        state[4 * i + 2] ^= (round_key[i] >> 8) as u8;
        state[4 * i + 3] ^= round_key[i] as u8;
    }
}

fn sbox(byte: u8) -> u8 {
    AES_SBOX[byte as usize]
}

/// replace each byte of the word with its substitution in the S-box
fn word_substitution(word: u32) -> u32 {
    let u0 = sbox((word >> 24) as u8) as u32;
    let u1 = sbox((word >> 16) as u8) as u32;
    let u2 = sbox((word >> 8) as u8) as u32;
    let u3 = sbox(word as u8) as u32;
    u0 << 24 | u1 << 16 | u2 << 8 | u3
}

fn bytes_substitution(state: &mut [u8; 16]) {
    for i in 0..16 {
        state[i] = sbox(state[i]);
    }
}

fn inv_sbox(byte: u8) -> u8 {
    INVERSE_AES_SBOX[byte as usize]
}

fn inv_byte_substitution(state: &mut [u8; 16]) {
    for i in 0..16 {
        state[i] = inv_sbox(state[i]);
    }
}

fn galois_field_mul(ap: u8, bp: u8) -> u8 {
    let mut res = 0;
    let mut a = ap;
    let mut b = bp;
    let mut high_bit;
    for _ in 0..8 {
        if (b & 1) == 1 {
            res ^= a;
        }
        high_bit = a & 0x80;
        a <<= 1;
        if high_bit == 0x80 {
            a ^= 0x1b;
        }
        b >>= 1;
    }
    res
}

fn shift_rows(state: &mut [u8; 16]) {
    let mut tem = vec![0; 16];
    for i in 0..16 {
        tem[i] = state[i];
    }
    for i in 1..4 {
        for j in 0..4 {
            // state[i][j] -> state[i][(i + j) % 4]
            state[i + 4 * j] = tem[(i + 4 * (j + i)) % 16];
        }
    }
}

fn inv_shift_rows(state: &mut [u8; 16]) {
    let mut tem = vec![0; 16];
    for i in 0..16 {
        tem[i] = state[i];
    }
    for i in 1..4 {
        for j in 0..4 {
            state[i + 4 * j] = tem[(i + 4 * (j + 4 - i)) % 16];
        }
    }
}

fn mix_column(state: &mut [u8; 16]) {
    let mut temp = vec![0; 16];
    for i in 0..16 {
        temp[i] = state[i];
    }
    for i in 0..4 {
        state[4 * i] = galois_field_mul(temp[4 * i], 2)
            ^ temp[2 + 4 * i]
            ^ temp[3 + 4 * i]
            ^ galois_field_mul(temp[1 + 4 * i], 3);
        state[1 + 4 * i] = galois_field_mul(temp[1 + 4 * i], 2)
            ^ temp[4 * i]
            ^ temp[3 + 4 * i]
            ^ galois_field_mul(temp[2 + 4 * i], 3);
        state[2 + 4 * i] = galois_field_mul(temp[2 + 4 * i], 2)
            ^ temp[1 + 4 * i]
            ^ temp[4 * i]
            ^ galois_field_mul(temp[3 + 4 * i], 3);
        state[3 + 4 * i] = galois_field_mul(temp[3 + 4 * i], 2)
            ^ temp[2 + 4 * i]
            ^ temp[1 + 4 * i]
            ^ galois_field_mul(temp[4 * i], 3);
    }
}

fn inv_mix_column(state: &mut [u8; 16]) {
    let mut temp = vec![0; 16];
    for i in 0..16 {
        temp[i] = state[i];
    }
    for i in 0..4 {
        state[4 * i] = galois_field_mul(temp[4 * i], 0x0e)
            ^ galois_field_mul(temp[1 + 4 * i], 0x0b)
            ^ galois_field_mul(temp[2 + 4 * i], 0x0d)
            ^ galois_field_mul(temp[3 + 4 * i], 9);
        state[1 + 4 * i] = galois_field_mul(temp[4 * i], 9)
            ^ galois_field_mul(temp[1 + 4 * i], 0x0e)
            ^ galois_field_mul(state[2 + 4 * i], 0x0b)
            ^ galois_field_mul(state[3 + 4 * i], 0x0d);
        state[2 + 4 * i] = galois_field_mul(temp[4 * i], 0x0d)
            ^ galois_field_mul(temp[1 + 4 * i], 9)
            ^ galois_field_mul(temp[2 + 4 * i], 0x0e)
            ^ galois_field_mul(temp[3 + 4 * i], 0x0b);
        state[3 + 4 * i] = galois_field_mul(temp[4 * i], 0x0b)
            ^ galois_field_mul(temp[1 + 4 * i], 0x0d)
            ^ galois_field_mul(temp[2 + 4 * i], 9)
            ^ galois_field_mul(temp[3 + 4 * i], 0x0e);
    }
}

fn g(word: u32, round: usize) -> u32 {
    // rotate: word = [u1 u2 u3 u4] => [u2 u3 u4 u1]
    let mut res = (word >> 24) | (word << 8);
    res = word_substitution(res);
    res ^= RC[round];
    res
}

fn key_schedule_aes(key: &AESKey) -> Vec<u32> {
    let nb = 4;
    // (the number of rounds, the number of words of key)
    let (nr, nk) = match key.key_type {
        AESKeyType::AES128 => (10, 4),
        AESKeyType::AES192 => (12, 6),
        AESKeyType::AES256 => (14, 8),
    };

    let mut w = vec![0u32; nb * (nr + 1)];
    for i in 0..nk {
        w[i] = (key.key[4 * i] as u32) << 24
            | (key.key[4 * i + 1] as u32) << 16
            | (key.key[4 * i + 2] as u32) << 8
            | (key.key[4 * i + 3] as u32);
    }

    for i in nk..nb * (nr + 1) {
        let mut temp = w[i - 1];
        if i % nk == 0 {
            temp = g(temp, i / nk);
        } else if (nk > 6) && (i % nk == 4) {
            temp = word_substitution(temp);
        }
        w[i] = w[i - nk] ^ temp;
    }
    w
}

/// PKCS7 padding
fn pad(input: &[u8]) -> Vec<u8> {
    let n = input.len();
    let padding_len = AES_BLOCK_SIZE - (n % AES_BLOCK_SIZE);
    let mut res = Vec::<u8>::with_capacity(n + padding_len);
    res.extend_from_slice(input);
    for _ in 0..padding_len {
        res.push(padding_len as u8);
    }
    res
}

fn unpad(input: &mut Vec<u8>) {
    let n = input.len();
    assert!(n >= AES_BLOCK_SIZE && n % AES_BLOCK_SIZE == 0);
    let pad_len = input[n - 1];
    assert!((pad_len as usize) < n);
    assert!(input
        .iter()
        .rev()
        .take(pad_len as usize)
        .all(|x| *x == pad_len));
    let unpad_len = n - (pad_len as usize);
    input.truncate(unpad_len);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn test_create_invalid_key() {
        let key = [0x7e, 0x2b, 0x15];
        let _ = AESKey::new(&key);
    }

    #[test]
    fn test_galois_field_mul() {
        assert_eq!(galois_field_mul(0xc2, 0x2f), 1);
        assert_eq!(galois_field_mul(0x53, 0xca), 1);
    }

    #[test]
    fn test_sbox() {
        assert_eq!(sbox(0x53), 0xed);
        assert_eq!(sbox(0xff), 0x16);
        assert_eq!(sbox(0xf0), 0x8c);
        assert_eq!(sbox(0x0f), 0x76);
        assert_eq!(sbox(0x00), 0x63);
    }

    #[test]
    fn test_sub_word() {
        assert_eq!(word_substitution(0x5355fcb0), 0xedfcb0e7);
        assert_eq!(word_substitution(0x00000000), 0x63636363);
        assert_eq!(word_substitution(0xffffffff), 0x16161616);
    }

    #[test]
    fn test_shift_rows() {
        let mut state = [
            0x00, 0x10, 0x20, 0x30, 0x01, 0x11, 0x21, 0x31, 0x02, 0x12, 0x22, 0x32, 0x03, 0x13,
            0x23, 0x33,
        ];
        let new_state = [
            0x00, 0x11, 0x22, 0x33, 0x01, 0x12, 0x23, 0x30, 0x02, 0x13, 0x20, 0x31, 0x03, 0x10,
            0x21, 0x32,
        ];
        shift_rows(&mut state);
        assert_eq!(state, new_state);
    }

    #[test]
    fn test_inv_shift_rows() {
        let new_state = [
            0x00, 0x10, 0x20, 0x30, 0x01, 0x11, 0x21, 0x31, 0x02, 0x12, 0x22, 0x32, 0x03, 0x13,
            0x23, 0x33,
        ];
        let mut state = [
            0x00, 0x11, 0x22, 0x33, 0x01, 0x12, 0x23, 0x30, 0x02, 0x13, 0x20, 0x31, 0x03, 0x10,
            0x21, 0x32,
        ];
        inv_shift_rows(&mut state);
        assert_eq!(state, new_state);
    }

    #[test]
    fn test_mix_column() {
        let mut state = [25; 16];
        let new_state = [25; 16];
        mix_column(&mut state);
        assert_eq!(state, new_state);
        inv_mix_column(&mut state);
        assert_eq!(state, new_state);
    }

    #[test]
    fn test_cipher_encrypt() {
        let aes128 = AESScheme::new(&[
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
            0x0e, 0x0f,
        ]);

        let aes192 = AESScheme::new(&[
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
            0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
        ]);

        let aes256 = AESScheme::new(&[
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
            0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b,
            0x1c, 0x1d, 0x1e, 0x1f,
        ]);

        let correct_output_128 = [
            0x69, 0xc4, 0xe0, 0xd8, 0x6a, 0x7b, 0x04, 0x30, 0xd8, 0xcd, 0xb7, 0x80, 0x70, 0xb4,
            0xc5, 0x5a,
        ];

        let correct_output_192 = [
            0xdd, 0xa9, 0x7c, 0xa4, 0x86, 0x4c, 0xdf, 0xe0, 0x6e, 0xaf, 0x70, 0xa0, 0xec, 0x0d,
            0x71, 0x91,
        ];

        let correct_output_256 = [
            0x8e, 0xa2, 0xb7, 0xca, 0x51, 0x67, 0x45, 0xbf, 0xea, 0xfc, 0x49, 0x90, 0x4b, 0x49,
            0x60, 0x89,
        ];

        let input = [
            0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd,
            0xee, 0xff,
        ];

        assert_eq!(aes128.encrypt_block(&input), correct_output_128);
        assert_eq!(aes192.encrypt_block(&input), correct_output_192);
        assert_eq!(aes256.encrypt_block(&input), correct_output_256);

        assert_eq!(aes128.encrypt_block_t_table(&input), correct_output_128);
        assert_eq!(aes192.encrypt_block_t_table(&input), correct_output_192);
        assert_eq!(aes256.encrypt_block_t_table(&input), correct_output_256);
    }

    #[test]
    fn test_cipher_decrypt() {
        let aes128 = AESScheme::new(&[
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
            0x0e, 0x0f,
        ]);

        let aes192 = AESScheme::new(&[
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
            0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
        ]);

        let aes256 = AESScheme::new(&[
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
            0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b,
            0x1c, 0x1d, 0x1e, 0x1f,
        ]);

        let input128 = [
            0x69, 0xc4, 0xe0, 0xd8, 0x6a, 0x7b, 0x04, 0x30, 0xd8, 0xcd, 0xb7, 0x80, 0x70, 0xb4,
            0xc5, 0x5a,
        ];
        let input192 = [
            0xdd, 0xa9, 0x7c, 0xa4, 0x86, 0x4c, 0xdf, 0xe0, 0x6e, 0xaf, 0x70, 0xa0, 0xec, 0x0d,
            0x71, 0x91,
        ];

        let input256 = [
            0x8e, 0xa2, 0xb7, 0xca, 0x51, 0x67, 0x45, 0xbf, 0xea, 0xfc, 0x49, 0x90, 0x4b, 0x49,
            0x60, 0x89,
        ];

        let output = [
            0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd,
            0xee, 0xff,
        ];

        assert_eq!(aes128.decrypt_block(&input128), output);
        assert_eq!(aes192.decrypt_block(&input192), output);
        assert_eq!(aes256.decrypt_block(&input256), output);
    }

    macro_rules! try_pad {
        ($n:expr, $len:expr) => {
            let arr = vec![$n; $len];
            let mut pm = pad(&arr);
            unpad(&mut pm);
            assert_eq!(pm, arr);
        };
    }

    #[test]
    fn test_pkcs7() {
        for i in 1..=254 {
            try_pad!(0, i);
            try_pad!(1, i);
            try_pad!(254, i);
            try_pad!(255, i);
        }
    }

    const NIST_AES_128_KEY: &[u8; 16] =
        b"\x2b\x7e\x15\x16\x28\xae\xd2\xa6\xab\xf7\x15\x88\x09\xcf\x4f\x3c";
    const NIST_AES_192_KEY: &[u8; 24] =
        b"\x8e\x73\xb0\xf7\xda\x0e\x64\x52\xc8\x10\xf3\x2b\x80\x90\x79\xe5\
      \x62\xf8\xea\xd2\x52\x2c\x6b\x7b";
    const NIST_AES_256_KEY: &[u8; 32] =
        b"\x60\x3d\xeb\x10\x15\xca\x71\xbe\x2b\x73\xae\xf0\x85\x7d\x77\x81\
      \x1f\x35\x2c\x07\x3b\x61\x08\xd7\x2d\x98\x10\xa3\x09\x14\xdf\xf4";
    const NIST_IV: &[u8; 16] = b"\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f";

    #[test]
    fn tes_aes_cbc() {
        let plaintext = b"\x6b\xc1\xbe\xe2\x2e\x40\x9f\x96\xe9\x3d\x7e\x11\x73\x93\x17\x2a\
                           \xae\x2d\x8a\x57\x1e\x03\xac\x9c\x9e\xb7\x6f\xac\x45\xaf\x8e\x51\
                           \x30\xc8\x1c\x46\xa3\x5c\xe4\x11\xe5\xfb\xc1\x19\x1a\x0a\x52\xef\
                           \xf6\x9f\x24\x45\xdf\x4f\x9b\x17\xad\x2b\x41\x7b\xe6\x6c\x37\x10";
        let ciphertext128 = b"\x76\x49\xab\xac\x81\x19\xb2\x46\xce\xe9\x8e\x9b\x12\xe9\x19\x7d\
                           \x50\x86\xcb\x9b\x50\x72\x19\xee\x95\xdb\x11\x3a\x91\x76\x78\xb2\
                           \x73\xbe\xd6\xb8\xe3\xc1\x74\x3b\x71\x16\xe6\x9e\x22\x22\x95\x16\
                           \x3f\xf1\xca\xa1\x68\x1f\xac\x09\x12\x0e\xca\x30\x75\x86\xe1\xa7";
        let ciphertext192 = b"\x4f\x02\x1d\xb2\x43\xbc\x63\x3d\x71\x78\x18\x3a\x9f\xa0\x71\xe8\
                           \xb4\xd9\xad\xa9\xad\x7d\xed\xf4\xe5\xe7\x38\x76\x3f\x69\x14\x5a\
                           \x57\x1b\x24\x20\x12\xfb\x7a\xe0\x7f\xa9\xba\xac\x3d\xf1\x02\xe0\
                           \x08\xb0\xe2\x79\x88\x59\x88\x81\xd9\x20\xa9\xe6\x4f\x56\x15\xcd";
        let ciphertext256 = b"\xf5\x8c\x4c\x04\xd6\xe5\xf1\xba\x77\x9e\xab\xfb\x5f\x7b\xfb\xd6\
                           \x9c\xfc\x4e\x96\x7e\xdb\x80\x8d\x67\x9f\x77\x7b\xc6\x70\x2c\x7d\
                           \x39\xf2\x33\x69\xa9\xd9\xba\xcf\xa5\x30\xe2\x63\x04\x23\x14\x61\
                           \xb2\xeb\x05\xe2\xc3\x9b\xe9\xfc\xda\x6c\x19\x07\x8c\x6a\x9d\x1b";
        let len_without_padding = 16 * 4;
        let padding_size = 16;

        let aes128 = AESScheme::new(NIST_AES_128_KEY);
        let mut iv = NIST_IV.clone();
        let encrypted = aes128.cbc_encrypt(&mut iv, plaintext);
        assert_eq!(encrypted.len(), len_without_padding + padding_size);
        assert_eq!(encrypted[..len_without_padding], ciphertext128[..]);
        iv = NIST_IV.clone();
        let decrypted = aes128.cbc_decrypt(&mut iv, &encrypted[..]);
        assert_eq!(decrypted[..], plaintext[..]);

        let aes192 = AESScheme::new(NIST_AES_192_KEY);
        iv = NIST_IV.clone();
        let encrypted = aes192.cbc_encrypt(&mut iv, plaintext);
        assert_eq!(encrypted.len(), len_without_padding + padding_size);
        assert_eq!(encrypted[..len_without_padding], ciphertext192[..]);
        iv = NIST_IV.clone();
        let decrypted = aes192.cbc_decrypt(&mut iv, &encrypted[..]);
        assert_eq!(decrypted[..], plaintext[..]);

        let aes256 = AESScheme::new(NIST_AES_256_KEY);
        iv = NIST_IV.clone();
        let encrypted = aes256.cbc_encrypt(&mut iv, plaintext);
        assert_eq!(encrypted.len(), len_without_padding + padding_size);
        assert_eq!(encrypted[..len_without_padding], ciphertext256[..]);
        iv = NIST_IV.clone();
        let decrypted = aes256.cbc_decrypt(&mut iv, &encrypted[..]);
        assert_eq!(decrypted[..], plaintext[..]);
    }
}

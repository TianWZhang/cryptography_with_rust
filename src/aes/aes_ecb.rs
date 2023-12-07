use crate::aes::constants::AES_BLOCK_SIZE;
use crate::aes::AESScheme;

/// assume the length of data is already a multiple of AES block size
/// no padding
pub fn ecb_encrypt(aes: &AESScheme, data: &[u8]) -> Vec<u8> {
    assert!(data.len() % AES_BLOCK_SIZE == 0);
    let mut res = data.to_vec();
    for block in res.chunks_mut(AES_BLOCK_SIZE) {
        let mut plaintext = [0u8; AES_BLOCK_SIZE];
        for i in 0..AES_BLOCK_SIZE {
            plaintext[i] = block[i];
        }
        let ciphertext = aes.encrypt_block(&plaintext);
        for i in 0..AES_BLOCK_SIZE {
            block[i] = ciphertext[i];
        }
    }
    res
}

pub fn ecb_decrypt(aes: &AESScheme, data: &[u8]) -> Vec<u8> {
    assert!(data.len() % AES_BLOCK_SIZE == 0);
    let mut res = data.to_vec();

    for block in res.chunks_mut(AES_BLOCK_SIZE) {
        let mut ciphertext = [0u8; AES_BLOCK_SIZE];
        for i in 0..AES_BLOCK_SIZE {
            ciphertext[i] = block[i];
        }
        let plaintext = aes.decrypt_block(&ciphertext);
        for i in 0..AES_BLOCK_SIZE {
            block[i] = plaintext[i];
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tes_aes_ecb() {
        // NIST ECB-AES-128 test vector F.1
        let plaintext = b"\x6b\xc1\xbe\xe2\x2e\x40\x9f\x96\xe9\x3d\x7e\x11\x73\x93\x17\x2a\
                                  \xae\x2d\x8a\x57\x1e\x03\xac\x9c\x9e\xb7\x6f\xac\x45\xaf\x8e\x51\
                                  \x30\xc8\x1c\x46\xa3\x5c\xe4\x11\xe5\xfb\xc1\x19\x1a\x0a\x52\xef\
                                  \xf6\x9f\x24\x45\xdf\x4f\x9b\x17\xad\x2b\x41\x7b\xe6\x6c\x37\x10";

        let key: [u8; 16] = [
            0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf,
            0x4f, 0x3c,
        ];

        let expected_ciphertext = b"\x3a\xd7\x7b\xb4\x0d\x7a\x36\x60\xa8\x9e\xca\xf3\x24\x66\xef\x97\
                                          \xf5\xd3\xd5\x85\x03\xb9\x69\x9d\xe7\x85\x89\x5a\x96\xfd\xba\xaf\
                                          \x43\xb1\xcd\x7f\x59\x8e\xce\x23\x88\x1b\x00\xe3\xed\x03\x06\x88\
                                          \x7b\x0c\x78\x5e\x27\xe8\xad\x3f\x82\x23\x20\x71\x04\x72\x5d\xd4";

        let aes128 = AESScheme::new(&key);
        let encrypted = ecb_encrypt(&aes128, plaintext);
        assert_eq!(encrypted, expected_ciphertext);
        let decrypted = ecb_decrypt(&aes128, &encrypted[..]);
        assert_eq!(decrypted, plaintext);
    }
}
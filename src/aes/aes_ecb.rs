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
    use hex::FromHex;

    use super::*;

    #[test]
    fn tes_aes_ecb() {
        // NIST ECB-AES-128 test vector F.1
        let plaintext = Vec::from_hex(
            "6bc1bee22e409f96e93d7e117393172a\
            ae2d8a571e03ac9c9eb76fac45af8e51\
            30c81c46a35ce411e5fbc1191a0a52ef\
            f69f2445df4f9b17ad2b417be66c3710",
        )
        .expect("Couldn't parse plaintext");

        let key = Vec::from_hex("2b7e151628aed2a6abf7158809cf4f3c").expect("Couldn't parse key");

        let expected_ciphertext = Vec::from_hex(
            "3ad77bb40d7a3660a89ecaf32466ef97\
             f5d3d58503b9699de785895a96fdbaaf\
             43b1cd7f598ece23881b00e3ed030688\
             7b0c785e27e8ad3f8223207104725dd4",
        )
        .expect("Could not parse ciphertext");

        let aes128 = AESScheme::new(&key);
        let encrypted = ecb_encrypt(&aes128, &plaintext);
        assert_eq!(encrypted, expected_ciphertext);
        let decrypted = ecb_decrypt(&aes128, &encrypted[..]);
        assert_eq!(decrypted, plaintext);
    }
}

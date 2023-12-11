use super::{constants::AES_BLOCK_SIZE, pad, unpad, AESScheme};

pub fn cbc_encrypt(aes: &AESScheme, iv: &mut [u8], data: &[u8]) -> Vec<u8> {
    let mut padded = pad(data);
    for block in padded.chunks_mut(AES_BLOCK_SIZE) {
        let mut plaintext = [0u8; AES_BLOCK_SIZE];
        for i in 0..AES_BLOCK_SIZE {
            plaintext[i] = block[i] ^ iv[i];
        }
        let ciphertext = aes.encrypt_block(&plaintext);
        for i in 0..AES_BLOCK_SIZE {
            block[i] = ciphertext[i];
            iv[i] = ciphertext[i];
        }
    }
    padded
}

pub fn cbc_decrypt(aes: &AESScheme, iv: &mut [u8], data: &[u8]) -> Vec<u8> {
    assert!(data.len() % AES_BLOCK_SIZE == 0);
    let mut res = data.to_vec();

    for block in res.chunks_mut(AES_BLOCK_SIZE) {
        let mut ciphertext = [0u8; AES_BLOCK_SIZE];
        for i in 0..AES_BLOCK_SIZE {
            ciphertext[i] = block[i];
        }
        let plaintext = aes.decrypt_block(&ciphertext);
        for i in 0..AES_BLOCK_SIZE {
            block[i] = plaintext[i] ^ iv[i];
            iv[i] = ciphertext[i];
        }
    }
    unpad(&mut res);
    res
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let encrypted = cbc_encrypt(&aes128, &mut iv, plaintext);
        assert_eq!(encrypted.len(), len_without_padding + padding_size);
        assert_eq!(encrypted[..len_without_padding], ciphertext128[..]);
        iv = NIST_IV.clone();
        let decrypted = cbc_decrypt(&aes128, &mut iv, &encrypted[..]);
        assert_eq!(decrypted[..], plaintext[..]);

        let aes192 = AESScheme::new(NIST_AES_192_KEY);
        iv = NIST_IV.clone();
        let encrypted = cbc_encrypt(&aes192, &mut iv, plaintext);
        assert_eq!(encrypted.len(), len_without_padding + padding_size);
        assert_eq!(encrypted[..len_without_padding], ciphertext192[..]);
        iv = NIST_IV.clone();
        let decrypted = cbc_decrypt(&aes192, &mut iv, &encrypted[..]);
        assert_eq!(decrypted[..], plaintext[..]);

        let aes256 = AESScheme::new(NIST_AES_256_KEY);
        iv = NIST_IV.clone();
        let encrypted = cbc_encrypt(&aes256, &mut iv, plaintext);
        assert_eq!(encrypted.len(), len_without_padding + padding_size);
        assert_eq!(encrypted[..len_without_padding], ciphertext256[..]);
        iv = NIST_IV.clone();
        let decrypted = cbc_decrypt(&aes256, &mut iv, &encrypted[..]);
        assert_eq!(decrypted[..], plaintext[..]);
    }
}

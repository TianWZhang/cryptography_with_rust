use super::AESScheme;

const R: u128 = 0xe1 << 120;

fn g128_mul(x: u128, y: u128) -> u128 {
    let mut z = 0;
    let mut v = y;

    for i in 0..128 {
        let xi = (x >> (127 - i)) & 1;
        if xi == 1 {
            z ^= v;
        }
        let lsb_v = v & 1;
        if lsb_v == 0 {
            v >>= 1;
        } else {
            v = (v >> 1) ^ R;
        }
    }
    z
}

fn ghash(hash_key: u128, x: &[u128]) -> u128 {
    let mut y = 0;
    let n = x.len();
    for i in 0..n {
        y = g128_mul(y ^ x[i], hash_key);
    }
    y
}

/// get the left most `s` bits of bytes
fn msb_s(s: usize, bytes: &[u8]) -> Vec<u8> {
    assert!(8 * bytes.len() >= s);
    let mut res = vec![];
    let num_bytes = s / 8;
    let bits_rem = s % 8;
    for i in 0..num_bytes {
        res.push(bytes[i]);
    }

    if bits_rem > 0 {
        res.push(bytes[num_bytes] >> (8 - bits_rem));
    }
    res
}

/// increment the right-most 32 bits of x, regarding as
/// the binary representation of an integer modulo 2^32
/// the remaining left-most 96 bits remain unchanged.
fn inc_32(x: u128) -> u128 {
    let msb = x >> 32;
    let mut lsb = x as u32;
    lsb = lsb.wrapping_add(1);
    msb << 32 | (lsb as u128)
}

/// bytes is in big endian order
fn bytes_to_128(bytes: &[u8]) -> u128 {
    let n = bytes.len();
    assert!(n <= 16);
    let mut res = 0;
    for i in 0..n {
        res |= (bytes[i] as u128) << (8 * (n - 1 - i));
    }
    res
}

fn gctr(aes: &AESScheme, icb: u128, x: &[u8]) -> Vec<u8> {
    let mut cb = icb; // counter block
    let mut res = vec![];
    for block in x.chunks(16) {
        if block.len() < 16 {
            let msb = msb_s(block.len() * 8, &aes.encrypt_block(&cb.to_be_bytes()));
            let y = bytes_to_128(block) ^ bytes_to_128(&msb);
            res.extend_from_slice(&y.to_be_bytes()[16 - block.len()..16]);
        } else {
            let y = bytes_to_128(block) ^ bytes_to_128(&aes.encrypt_block(&cb.to_be_bytes()));
            res.extend_from_slice(&y.to_be_bytes());
        }
        cb = inc_32(cb);
    }
    res
}

/// convert bytes to blocks, block size is 16 bytes
/// pad if needed
fn bytes_to_blocks(bytes: &[u8]) -> Vec<u128> {
    let blocks = bytes.len() / 16;
    let mut res = vec![];
    for i in 0..blocks {
        res.push(bytes_to_128(&bytes[i * 16..(i + 1) * 16]));
    }
    if bytes.len() % 16 > 0 {
        let mut rem_block = vec![];
        rem_block.extend_from_slice(&bytes[blocks * 16..bytes.len()]);
        for _ in 0..(16 - bytes.len() % 16) {
            rem_block.push(0);
        }
        res.push(bytes_to_128(&rem_block));
    }
    res
}

pub fn gcm_authenticated_encrypt(
    aes: &AESScheme,
    iv: &[u8],
    plaintext: &[u8],
    additional_data: &[u8],
    tag_len: usize,
) -> (Vec<u8>, Vec<u8>) {
    let hash_key = bytes_to_128(&aes.encrypt_block(&[0u8; 16]));
    let j0 = if iv.len() == 12 {
        (bytes_to_128(iv) << 32) | 0x00000001
    } else {
        let mut pad_iv = bytes_to_blocks(iv);
        pad_iv.push((iv.len() * 8) as u128);
        ghash(hash_key, &pad_iv)
    };
    let ciphertext = gctr(aes, inc_32(j0), plaintext);

    let mut pad_data = bytes_to_blocks(&additional_data);
    pad_data.extend(bytes_to_blocks(&ciphertext));
    pad_data.push((((additional_data.len() * 8) as u128) << 64) | (ciphertext.len() * 8) as u128);

    let s = ghash(hash_key, &pad_data);
    let tag = msb_s(tag_len, &gctr(aes, j0, &s.to_be_bytes()));
    (ciphertext, tag)
}

pub fn gcm_authenticated_decrypt(
    aes: &AESScheme,
    iv: &[u8],
    ciphertext: &[u8],
    additional_data: &[u8],
    tag: &[u8],
    tag_len: usize,
) -> Result<Vec<u8>, String> {
    if tag.len() * 8 != tag_len {
        return Err("incorrect tag length".to_string());
    }
    let hash_key = bytes_to_128(&aes.encrypt_block(&[0u8; 16]));
    let j0 = if iv.len() == 12 {
        (bytes_to_128(iv) << 32) | 0x00000001
    } else {
        let mut pad_iv = bytes_to_blocks(iv);
        pad_iv.push((iv.len() * 8) as u128);
        ghash(hash_key, &pad_iv)
    };
    let plaintext = gctr(aes, inc_32(j0), ciphertext);

    let mut pad_data = bytes_to_blocks(&additional_data);
    pad_data.extend(bytes_to_blocks(&ciphertext));
    pad_data.push((((additional_data.len() * 8) as u128) << 64) | (ciphertext.len() * 8) as u128);
    let s = ghash(hash_key, &pad_data);
    let computed_tag = msb_s(tag_len, &gctr(aes, j0, &s.to_be_bytes()));

    if tag == computed_tag {
        Ok(plaintext)
    } else {
        Err("incorrect unthenticated tag".to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bytes_to_128() {
        let mut bytes = [
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x01,
        ];
        assert_eq!(bytes_to_128(&bytes), 1);
        bytes = [
            0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00,
        ];
        assert_eq!(bytes_to_128(&bytes), 2u128.pow(120));
    }

    #[test]
    fn test_inc_32() {
        let mut x = 0xf;
        assert_eq!(inc_32(x), 0x10);

        x = 0;
        assert_eq!(inc_32(x), 1);

        x = 0xffffffff;
        assert_eq!(inc_32(x), 0);

        x = 0xeffffffff;
        assert_eq!(inc_32(x), 0xe00000000);
    }

    #[test]
    fn test_bytes_to_blocks() {
        let bytes = [
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x01, 0x01,
        ];
        assert_eq!(bytes_to_blocks(&bytes), [1, 2u128.pow(120)]);

        let bytes = [
            0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x01,
        ];
        assert_eq!(bytes_to_blocks(&bytes), [2u128.pow(120), 1, 256]);
    }

    #[test]
    fn test_gcm_ae() {
        use hex::FromHex;

        // NIST test vector 2
        let key = Vec::from_hex("fe47fcce5fc32665d2ae399e4eec72ba").expect("Couldn't parse key");
        let iv = Vec::from_hex("5adb9609dbaeb58cbd6e7275").expect("Couldn't parse IV");
        let expected_ct = Vec::from_hex("98f4826f05a265e6dd2be82db241c0fbbbf9ffb1c173aa83964b7cf5393043736365253ddbc5db8778371495da76d269e5db3e").expect("Couldn't parse ciphertext");
        let aad =
            Vec::from_hex("88319d6e1d3ffa5f987199166c8a9b56c2aeba5a").expect("Couldn't parse AAD");
        let expected_tag =
            Vec::from_hex("291ef1982e4defedaa2249f898556b47").expect("Couldn't parse tag");
        let pt = Vec::from_hex("7c0e88c88899a779228465074797cd4c2e1498d259b54390b85e3eef1c02df60e743f1b840382c4bccaf3bafb4ca8429bea063").expect("Couldn't parse plaintext");
        let tag_length = 128;

        let aes128 = AESScheme::new(&key);
        let (test_ct, test_tag) = gcm_authenticated_encrypt(&aes128, &iv, &pt, &aad, tag_length);

        assert_eq!(test_ct, expected_ct);
        assert_eq!(test_tag, expected_tag);
    }

    #[test]
    fn test_gcm_ad() {
        use hex::FromHex;
        // authenticated decryption that fails
        let mut key =
            Vec::from_hex("f5a0b1639c67c7760109056a3a329804").expect("Couldn't parse key");
        let mut iv = Vec::from_hex("e1b75506d66509a52f0960f7").expect("Couldn't parse IV");
        let mut ct = Vec::from_hex("4d8738341660f7e49ca1ddf7db1255c1eca46b947fa80134340d364e611255194f3261413a82e763720ef81dedc8b10bed3b30").expect("Couldn't parse ciphertext");
        let mut aad =
            Vec::from_hex("8421f67419d3d37cc9e97b712b8b0924").expect("Couldn't parse AAD");
        let mut tag = Vec::from_hex("d7c586892b2e6ad60c2106a8").expect("Couldn't parse tag");
        let mut tag_length = 96;

        let mut aes128 = AESScheme::new(&key);
        let result = gcm_authenticated_decrypt(&aes128, &iv, &ct, &aad, &tag, tag_length);
        assert!(result.is_err());

        // authenticated decryption test that passes
        // NIST test vector 2
        key = Vec::from_hex("fe47fcce5fc32665d2ae399e4eec72ba").expect("Couldn't parse key");
        iv = Vec::from_hex("5adb9609dbaeb58cbd6e7275").expect("Couldn't parse IV");
        ct = Vec::from_hex("98f4826f05a265e6dd2be82db241c0fbbbf9ffb1c173aa83964b7cf5393043736365253ddbc5db8778371495da76d269e5db3e").expect("Couldn't parse ciphertext");
        aad =
            Vec::from_hex("88319d6e1d3ffa5f987199166c8a9b56c2aeba5a").expect("Couldn't parse AAD");
        tag = Vec::from_hex("291ef1982e4defedaa2249f898556b47").expect("Couldn't parse tag");
        let expected_pt = Vec::from_hex("7c0e88c88899a779228465074797cd4c2e1498d259b54390b85e3eef1c02df60e743f1b840382c4bccaf3bafb4ca8429bea063").expect("Couldn't parse plaintext");
        tag_length = 128;

        aes128 = AESScheme::new(&key);
        let result = gcm_authenticated_decrypt(&aes128, &iv, &ct, &aad, &tag, tag_length);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), expected_pt);
    }
}

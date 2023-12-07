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

fn ghash(H: u128, x: &[u128]) -> u128 {
    let mut y = 0;
    let n = x.len();
    for i in 0..n {
        y = g128_mul(y ^ x[i], H);
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
fn u8s_to_u128(bytes: &[u8]) -> u128 {
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
            let y = u8s_to_u128(block) ^ u8s_to_u128(&msb);
            res.extend_from_slice(&y.to_be_bytes()[16 - block.len()..16]);
        } else {
            let y = u8s_to_u128(block) ^ u8s_to_u128(&aes.encrypt_block(&cb.to_be_bytes()));
            res.extend_from_slice(&y.to_be_bytes());
        }
        cb = inc_32(cb);
    }
    res
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_u8s_to_u128() {
        let mut bytes = [
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01
        ];
        assert_eq!(u8s_to_u128(&bytes), 1);
        bytes = [
            0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
        ];
        assert_eq!(u8s_to_u128(&bytes), 2u128.pow(120));
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
}
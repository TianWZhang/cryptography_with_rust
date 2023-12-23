use sha3::{Digest, Sha3_256, Sha3_512, Shake128, Shake256};

/// shake-128 wrapper
/// SHAKE128 is an extendable-output function (XOF) in the SHA-3 family. As a XOF, SHAKE128
/// is a generalization of a cryptographic hash function. Instead of creating a fixed-length digest,
/// it can produce outputs of any desired length.
// The 128 in its name indicates its maximum security level (in bits).
pub fn shake_128(data: &[u8], len: usize) -> Vec<u8> {
    use sha3::digest::{ExtendableOutput, Update, XofReader};
    let mut hasher = Shake128::default();
    hasher.update(data);
    let mut reader = hasher.finalize_xof();

    let mut buffer = vec![0; len];
    reader.read(&mut buffer);
    buffer
}

/// shake-256 wrapper
pub fn shake_256(data: &[u8], len: usize) -> Vec<u8> {
    use sha3::digest::{ExtendableOutput, Update, XofReader};
    let mut hasher = Shake256::default();
    hasher.update(data);
    let mut reader = hasher.finalize_xof();

    let mut buffer = vec![0; len];
    reader.read(&mut buffer);
    buffer
}

/// sha3-256 wrapper
pub fn sha3_256(data: &[u8]) -> Vec<u8> {
    let mut hasher = Sha3_256::new();

    // write input message
    hasher.update(data);

    // read hash digest
    hasher.finalize().to_vec()
}

/// sha3-512 wrapper
pub fn sha3_512(data: &[u8]) -> Vec<u8> {
    let mut hasher = Sha3_512::new();
    hasher.update(data);
    hasher.finalize().to_vec()
}

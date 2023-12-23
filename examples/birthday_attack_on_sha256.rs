use cryptography_with_rust::hash::Sha256;
use std::collections::HashMap;

const BYTES: usize = 2;

/// Truncate the output of SHA-256(message) to n bytes
fn sha_256_n(msg: &[u8], n: usize) -> Vec<u8> {
    assert!((1..=8).contains(&n));
    let mut hasher = Sha256::new();
    hasher.update(msg);
    let digest = hasher.digest();
    let mut res = vec![];
    if n >= 4 {
        res.extend(digest[0].to_be_bytes());
        let second_bytes = digest[1].to_be_bytes();
        for i in 4..n {
            res.push(second_bytes[i - 4]);
        }
    } else {
        let first_bytes = digest[0].to_be_bytes();
        for i in 0..n {
            res.push(first_bytes[i]);
        }
    }
    res
}

fn birthday_attack_on_sha256n(n_bytes: usize) -> (usize, usize) {
    assert!((1..=8).contains(&n_bytes));
    let mut map = HashMap::new();
    // return a collision pair (i1, i2)
    (1usize..100000)
        .find_map(|i| {
            //Applies function to the elements of iterator and returns
            // the first non-none result.
            let hash = sha_256_n(&i.to_be_bytes(), n_bytes);
            //If the map did have this key present, the value is updated, and the old value is returned.
            map.insert(hash.clone(), i).map(|i2| (i, i2))
        })
        .unwrap()
}

fn preimage_sha256_2(target: &[u8], start: usize) -> Option<usize> {
    assert!(target.len() == BYTES);
    (start..(start + 1000000)).find(|i| {
        let hash = sha_256_n(&i.to_be_bytes(), BYTES);
        hash == target
    })
}

fn main() {
    let (i1, i2) = birthday_attack_on_sha256n(BYTES);
    let i1_hash = sha_256_n(&i1.to_be_bytes(), BYTES);
    println!(
        "i1 = {}, the first 2 bytes of the sha256 hash value of i1: {}",
        i1,
        format!("{:08X?}{:08X?}", i1_hash[0], i1_hash[1])
    );
    let i2_hash = sha_256_n(&i2.to_be_bytes(), BYTES);
    println!(
        "i2 = {}, the first 2 bytes of the sha256 hash value of i2: {}",
        i2,
        format!("{:08X?}{:08X?}", i2_hash[0], i2_hash[1])
    );

    let s = [0x3D, 0x4B];
    if let Some(n) = preimage_sha256_2(&s, 0) {
        println!("The algorithm need {} tries", n);
    } else {
        println!("Could not find the preimage")
    }
}

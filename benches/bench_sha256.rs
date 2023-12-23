use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cryptography_with_rust::hash::sha256;
use ring::digest::SHA256;

const TEXT: &[u8] = b"some text to test hash algorithms";

// 304.13 ns on my machine
pub fn ben_my_sha2(c: &mut Criterion) {
    c.bench_function("my sha256", |b| {
        b.iter(|| {
            let mut hasher = sha256::Sha256::new();
            hasher.update(black_box(TEXT));
            hasher.digest();
        })
    });
}

// 68.468 ns on my machine
pub fn ben_ring_sha256(c: &mut Criterion) {
    c.bench_function("ring sha256", |b| {
        b.iter(|| {
            ring::digest::digest(&SHA256, black_box(TEXT));
        })
    });
}

criterion_group!(benches, ben_my_sha2, ben_ring_sha256);
criterion_main!(benches);

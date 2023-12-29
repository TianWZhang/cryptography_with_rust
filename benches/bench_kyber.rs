use criterion::{criterion_group, criterion_main, Criterion};
use cryptography_with_rust::kyber::{
    KYBER1024KEM, KYBER1024PKE, KYBER512KEM, KYBER512PKE, KYBER768KEM, KYBER768PKE,
};
use rand::RngCore;

fn random_bytes(len: usize) -> Vec<u8> {
    let mut data = vec![0; len];
    let mut rng = rand::thread_rng();
    rng.fill_bytes(&mut data);
    data
}

//  ns on my machine
pub fn bench_kyber512_pke(c: &mut Criterion) {
    let pke = KYBER512PKE;
    let msg = random_bytes(32);
    let random_coins = random_bytes(32);

    let mut group = c.benchmark_group("KYBER 512 PKE");
    let (sk, pk) = pke.key_gen();
    let ciphertext = pke.encrypt(&pk, &msg, &random_coins);
    let _msg_dec = pke.decrypt(&sk, &ciphertext);

    group.bench_function("KeyGen", |b| b.iter(|| pke.key_gen()));
    group.bench_function("Encrypt", |b| {
        b.iter(|| pke.encrypt(&pk, &msg, &random_coins))
    });
    group.bench_function("Decrypt", |b| b.iter(|| pke.decrypt(&sk, &ciphertext)));
    group.finish();
}

pub fn bench_kyber768_pke(c: &mut Criterion) {
    let pke = KYBER768PKE;
    let msg = random_bytes(32);
    let random_coins = random_bytes(32);

    let mut group = c.benchmark_group("KYBER 768 PKE");
    let (sk, pk) = pke.key_gen();
    let ciphertext = pke.encrypt(&pk, &msg, &random_coins);
    let _msg_dec = pke.decrypt(&sk, &ciphertext);

    group.bench_function("KeyGen", |b| b.iter(|| pke.key_gen()));
    group.bench_function("Encrypt", |b| {
        b.iter(|| pke.encrypt(&pk, &msg, &random_coins))
    });
    group.bench_function("Decrypt", |b| b.iter(|| pke.decrypt(&sk, &ciphertext)));
    group.finish();
}

pub fn bench_kyber1024_pke(c: &mut Criterion) {
    let pke = KYBER1024PKE;
    let msg = random_bytes(32);
    let random_coins = random_bytes(32);

    let mut group = c.benchmark_group("KYBER 1024 PKE");
    let (sk, pk) = pke.key_gen();
    let ciphertext = pke.encrypt(&pk, &msg, &random_coins);
    let _msg_dec = pke.decrypt(&sk, &ciphertext);

    group.bench_function("KeyGen", |b| b.iter(|| pke.key_gen()));
    group.bench_function("Encrypt", |b| {
        b.iter(|| pke.encrypt(&pk, &msg, &random_coins))
    });
    group.bench_function("Decrypt", |b| b.iter(|| pke.decrypt(&sk, &ciphertext)));
    group.finish();
}

//  ns on my machine
pub fn bench_kyber512_kem(c: &mut Criterion) {
    let kem = KYBER512KEM;
    let mut group = c.benchmark_group("KYBER 512 KEM");
    let (sk, pk) = kem.key_gen();
    let (ciphertext, _shared_key) = kem.encapsulate(&pk);
    let _shared_key = kem.decapsulate(&ciphertext, &sk);

    group.bench_function("KeyGen", |b| b.iter(|| kem.key_gen()));
    group.bench_function("Encapsulate", |b| b.iter(|| kem.encapsulate(&pk)));
    group.bench_function("Decapsulate", |b| {
        b.iter(|| kem.decapsulate(&ciphertext, &sk))
    });
    group.finish();
}

pub fn bench_kyber768_kem(c: &mut Criterion) {
    let kem = KYBER768KEM;
    let mut group = c.benchmark_group("KYBER 768 KEM");
    let (sk, pk) = kem.key_gen();
    let (ciphertext, _shared_key) = kem.encapsulate(&pk);
    let _shared_key = kem.decapsulate(&ciphertext, &sk);

    group.bench_function("KeyGen", |b| b.iter(|| kem.key_gen()));
    group.bench_function("Encapsulate", |b| b.iter(|| kem.encapsulate(&pk)));
    group.bench_function("Decapsulate", |b| {
        b.iter(|| kem.decapsulate(&ciphertext, &sk))
    });
    group.finish();
}

pub fn bench_kyber1024_kem(c: &mut Criterion) {
    let kem = KYBER1024KEM;
    let mut group = c.benchmark_group("KYBER 1024 KEM");
    let (sk, pk) = kem.key_gen();
    let (ciphertext, _shared_key) = kem.encapsulate(&pk);
    let _shared_key = kem.decapsulate(&ciphertext, &sk);

    group.bench_function("KeyGen", |b| b.iter(|| kem.key_gen()));
    group.bench_function("Encapsulate", |b| b.iter(|| kem.encapsulate(&pk)));
    group.bench_function("Decapsulate", |b| {
        b.iter(|| kem.decapsulate(&ciphertext, &sk))
    });
    group.finish();
}

criterion_group! {
    name = kyber512;
    config = Criterion::default().sample_size(100);
    targets = bench_kyber512_pke, bench_kyber512_kem
}

criterion_group! {
    name = kyber768;
    config = Criterion::default().sample_size(100);
    targets = bench_kyber768_pke, bench_kyber768_kem
}
criterion_group! {
    name = kyber1024;
    config = Criterion::default().sample_size(100);
    targets = bench_kyber1024_pke, bench_kyber1024_kem
}

criterion_main!(kyber512, kyber768, kyber1024);

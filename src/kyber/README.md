This is an educational implementation of the CRYSTALS-KYBER lattice-based key encapsulation suite (KYBER [1]) in Rust Programming Language, which won the NIST competition for the first post-quantum cryptography (PQ) standard [2].

This implementation prioritizes accuracy but sacrifices security and efficiency for simplicity. It is NOT suitable for production use. In particular, it is not guaranteed to be constant time.

### Supported features and algorithms

* Key encapsulation mechanism (`KEM`)
* Public-key encryption (`PKE`)
* All the parameters described in the NIST submission: `kyber-512`, `kyber-768` and `kyber-1024`.

The 3rd round updated specification (4 August 2021) is used as a basis for implementation.

### Dev options
#### Test
```cargo test```


#### Benchmark
Benchmarks with criterion:

```cargo bench --bench bench_kyber```


### References and documentation

* https://pq-crystals.org/kyber/
* https://csrc.nist.gov/Projects/Post-Quantum-Cryptography
* https://github.com/rust-crypto-labs/kybe-rs/tree/master

[1]: https://pq-crystals.org/kyber/
[2]: https://csrc.nist.gov/Projects/Post-Quantum-Cryptography
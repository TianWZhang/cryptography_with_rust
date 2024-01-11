use num_bigint::BigUint;
use rand::Rng;

fn prime_factors(mut n: u64) -> Vec<u64> {
    let mut res = vec![];
    let mut i = 2;
    while i * i <= n {
        if n % i == 0 {
            res.push(i);
            n /= i;
            while n % i == 0 {
                n /= i;
            }
        }
        i += 1;
    }
    if n > 1 {
        res.push(n);
    }
    res
}

pub(crate) fn gcd(mut n: u64, mut m: u64) -> u64 {
    assert!(n != 0 && m != 0);
    while m != 0 {
        if m < n {
            std::mem::swap(&mut m, &mut n);
        }
        m %= n;
    }
    n
}

/// miller rabin prime test
pub(crate) fn is_prime(n: u64) -> bool {
    if n == 2 || n == 3 {
        return true;
    }
    if n % 2 == 0 || n == 1 {
        return false;
    }

    // n - 1 = 2^k * q, with q odd
    let (k, q) = {
        let mut k = 0;
        let mut q = n - 1;
        while q % 2 == 0 {
            k += 1;
            q /= 2;
        }
        (k, q)
    };

    // 3. Looping 200 times, choose a random value for the witness a, with 2= < a < n-1.
    let mut rng = rand::thread_rng();
    for _ in 0..200 {
        let a = rng.gen_range(2..n - 1);

        if gcd(a, n) != 1 {
            return false;
        }

        // a = a^q mod n. If a == 1, no information.
        let mut a_big = BigUint::from(a);
        a_big = a_big.modpow(&BigUint::from(q), &BigUint::from(n));
        if a_big == BigUint::from(1u64) {
            continue;
        }

        // for i in 0..=k-1, loop:
        let mut unbroken = true;
        for _ in 0..k {
            if a_big == BigUint::from(n - 1) {
                unbroken = false;
                break;
            }
            a_big = a_big.modpow(&BigUint::from(2u64), &BigUint::from(n));
        }
        // a \in 2..n-1 such that a^q != 1 mod n, and for all i \in 0..k, a^{2^i * q} != n - 1 mod n,
        // hence a is a Miller-Rabin witness for n => n is a composite number
        if unbroken {
            return false;
        }
    }

    true
}

fn is_primitive_root(val: u64, p: u64) -> bool {
    let val = BigUint::from(val);
    if val.modpow(&BigUint::from(p - 1), &BigUint::from(p)) != BigUint::from(1u64) {
        return false;
    }
    for i in prime_factors(p - 1) {
        if val.modpow(&BigUint::from((p - 1) / i), &BigUint::from(p)) == BigUint::from(1u64) {
            return false;
        }
    }
    true
}

fn get_generator(p: u64) -> u64 {
    if !is_prime(p) {
        panic!("P is not a prime number");
    }

    if p == 3329 {
        return 3;
    }
    for i in 1..p {
        if is_primitive_root(i, p) {
            return i;
        }
    }
    unreachable!("P is not a prime number")
}

/// p mod n == 1
/// get n-th primitive root of unity
pub(crate) fn get_primitive_root_of_unity(n: u64, p: u64) -> u64 {
    if p % n != 1 {
        panic!("p mod n == 1 not satisfied");
    }

    if p == 3329 && n == 256 {
        return 17;
    }

    let g = get_generator(p);
    let power = (p - 1) / n as u64;
    BigUint::from(g)
        .modpow(&BigUint::from(power), &BigUint::from(p))
        .try_into()
        .unwrap()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_prime_factors() {
        let mut n = 128;
        assert_eq!(prime_factors(n), vec![2]);

        n = 182;
        assert_eq!(prime_factors(n), vec![2, 7, 13]);

        n = 1;
        assert_eq!(prime_factors(n), vec![]);
    }

    #[test]
    fn test_is_prime() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(!is_prime(4));
        assert!(is_prime(5));
        assert!(!is_prime(6));
        assert!(is_prime(7));
        assert!(!is_prime(8));
        assert!(!is_prime(9));
        assert!(!is_prime(10));
        assert!(is_prime(11));
        assert!(!is_prime(12));
        assert!(is_prime(13));
        assert!(!is_prime(14));
        assert!(!is_prime(15));
        assert!(!is_prime(16));
        assert!(is_prime(17));
        assert!(!is_prime(18));
        assert!(is_prime(19));
        assert!(is_prime(10001231));
        assert!(is_prime(100001029));
        assert!(is_prime(2305843009214414849));
    }

    #[test]
    fn test_generator() {
        assert_eq!(get_generator(17), 3);
        assert_eq!(get_generator(113), 3);
        assert_eq!(get_generator(22273), 5);
        assert_eq!(get_generator(31489), 7);
        assert_eq!(get_generator(26881), 11);
        assert_eq!(get_generator(3329), 3);
    }
}

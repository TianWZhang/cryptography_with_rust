use rand::Rng;

pub(crate) fn prime_factors(mut n: u64) -> Vec<u64> {
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

pub(crate) fn modpow(mut a: u64, mut n: u64, p: u64) -> u64 {
    let mut res = 1;
    while n > 0 {
        if n % 2 == 1 {
            res = (res * a) % p;
        }
        a = (a * a) % p;
        n /= 2;
    }
    res
}

/// miller rabin prime test
pub(crate) fn is_prime(n: u64) -> bool {
    if n == 2 || n == 3 {
        return true;
    }
    if n % 2 == 0 || n == 1 {
        return false;
    }

    // n-1 = 2^k * q, with q odd
    let (k, q) = {
        let mut k = 0;
        let mut q = n - 1;
        while q % 2 == 0 {
            k += 1;
            q /= 2;
        }
        (k, q)
    };

    // 3. Looping 200 times, choose a random value for the witness a, with 2 < a < n-2.
    let mut rng = rand::thread_rng();
    for _ in 0..200 {
        let mut a = rng.gen_range(2..n - 1);

        if gcd(a, n) != 1 {
            return false;
        }

        // a = a^q mod n. If a == 1, no information.
        a = modpow(a, q, n);
        if a == 1 {
            continue;
        }

        // for i in 0..=k-1, loop:
        let mut unbroken = true;
        for _ in 0..=(k - 1) {
            if a == n - 1 {
                unbroken = false;
                break;
            }
            a = modpow(a, 2, n);
        }
        if unbroken {
            return false;
        }
    }

    true
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
    }
}

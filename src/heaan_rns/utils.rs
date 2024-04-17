use num_complex::Complex64;
use rand::Rng;

#[allow(dead_code)]
pub(crate) fn equal_up_to_epsilon(nums1: &[Complex64], nums2: &[Complex64], epsilon: f64) -> bool {
    for (x1, x2) in nums1.iter().zip(nums2.iter()) {
        if (x1.re - x2.re).abs() > epsilon || (x1.im - x2.im).abs() > epsilon {
            return false;
        }
    }
    true
}

#[allow(dead_code)]
pub(crate) fn gen_random_complex_vector(length: usize) -> Vec<Complex64> {
    let mut rng = rand::thread_rng();
    (0..length)
        .map(|_| {
            let re = rng.gen_range(-1.0..1.0);
            let im = rng.gen_range(-1.0..1.0);
            Complex64::new(re, im)
        })
        .collect()
}

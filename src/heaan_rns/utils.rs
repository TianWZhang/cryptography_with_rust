use num_complex::Complex64;

pub fn equal_up_to_epsilon(nums1: &[Complex64], nums2: &[Complex64], epsilon: f64) -> bool {
    for (x1, x2) in nums1.iter().zip(nums2.iter()) {
        if (x1.re - x2.re).abs() > epsilon || (x1.im - x2.im).abs() > epsilon {
            return false;
        }
    }
    true
}

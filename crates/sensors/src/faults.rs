use rand::Rng;
use rand_distr::{Distribution, Normal};

pub fn apply_gaussian(value: f64, mean: f64, std_dev: f64, rng: &mut impl Rng) -> f64 {
    let normal = Normal::new(mean, std_dev).expect("Invalid mean or std_dev");
    value + normal.sample(rng)
}

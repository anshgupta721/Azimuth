use rand::SeedableRng;
use rand_chacha::ChaCha12Rng;
use rand_distr::{Distribution, Normal};

pub fn apply_gaussian(value: f64, mean: f64, std_dev: f64, seed: u64) -> f64 {
    let mut rng = ChaCha12Rng::seed_from_u64(seed);
    let normal = Normal::new(mean, std_dev).expect("Invalid mean or std_dev");
    value + normal.sample(&mut rng)
}

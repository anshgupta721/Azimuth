use nalgebra::{SMatrix, SVector};

pub fn squared_mahalanobis_distance<const N: usize>(
    error: SVector<f64, N>,
    covariance: SMatrix<f64, N, N>,
) -> Option<f64> {
    let chol = covariance.cholesky()?;
    let y = chol.l().solve_lower_triangular(&error)?;
    Some(y.norm_squared())
}

use nalgebra::{SMatrix, SVector};

use crate::gravity::spherical_harmonics::{HarmonicCoeffs, LegendreCache};
pub trait GravityField {
    fn acceleration(&self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3>;
    /// Optional: f
    fn potential(&self, pos_body_fixed: SVector<f64, 3>) -> Option<f64>;

    fn gravity_gradient(&self, pos_body_fixed: SVector<f64, 3>) -> SMatrix<f64, 3, 3> {
        let mut gradient = SMatrix::<f64, 3, 3>::zeros();
        let h = SVector::<f64, 3>::new(1.0, 1.0, 1.0); // Step Size

        let accel_plus = self.acceleration(pos_body_fixed + h);
        let accel_minus = self.acceleration(pos_body_fixed - h);

        for i in 0..3 {
            gradient.set_column(i, &(accel_plus - accel_minus).component_div(&(2.0 * h)));
        }

        gradient
    }
}

pub struct PointMass {
    pub mu: f64,
}

impl GravityField for PointMass {
    fn acceleration(&self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3> {
        let r = pos_body_fixed.norm();
        let common = -self.mu / r.powf(3.0);
        pos_body_fixed * common
    }

    fn potential(&self, pos_body_fixed: SVector<f64, 3>) -> Option<f64> {
        let r = pos_body_fixed.norm();
        Some((pos_body_fixed * (self.mu / r)).norm())
    }
}

pub struct SphericalHarmonics {
    mu: f64,
    r_ref: f64,
    pub coeffs: HarmonicCoeffs,
}

impl SphericalHarmonics {
    pub fn new(mu: f64, r_ref: f64, coeffs: HarmonicCoeffs) -> Self {
        SphericalHarmonics { mu, r_ref, coeffs }
    }

    fn spherical(pos: SVector<f64, 3>) -> (f64, f64, f64) {
        let r = pos.norm();
        let phi = (pos[2] / r).asin();
        let lambda = pos[1].atan2(pos[0]);
        (r, phi, lambda)
    }
}

impl GravityField for SphericalHarmonics {
    fn acceleration(&self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3> {
        let (r, phi, lambda) = Self::spherical(pos_body_fixed);
        let leg = LegendreCache::new(self.coeffs.degree(), phi);
        let (sin_phi, cos_phi) = phi.sin_cos();
        let (sin_lam, cos_lam) = lambda.sin_cos();

        let mut dv_dr = 0.0;
        let mut dv_dphi = 0.0;
        let mut dv_dlambda = 0.0;

        for n in 0..=self.coeffs.degree() {
            let rn = (self.r_ref / r).powi(n as i32);
            for m in 0..=n {
                let (sin_ml, cos_ml) = (m as f64 * lambda).sin_cos();
                let cs = self.coeffs.get_c(n, m) * cos_ml + self.coeffs.get_s(n, m) * sin_ml;
                let cs_dl = -self.coeffs.get_c(n, m) * sin_ml + self.coeffs.get_s(n, m) * cos_ml;

                dv_dr += (n as f64 + 1.0) * rn * leg.p[n][m] * cs;
                dv_dphi += rn * leg.dp[n][m] * cs;
                dv_dlambda += rn * leg.p[n][m] * (m as f64) * cs_dl;
            }
        }

        let dv_dr = -self.mu / (r * r) * dv_dr;
        let dv_dphi = self.mu / r * dv_dphi;
        let dv_dlambda = self.mu / r * dv_dlambda;

        // Spherical -> local components
        let a_r = dv_dr;
        let a_phi = dv_dphi / r;

        // Note: singular at the poles (cos_phi -> 0)
        let a_lambda = dv_dlambda / (r * cos_phi);

        let ax = (a_r * cos_phi - a_phi * sin_phi) * cos_lam - a_lambda * sin_lam;
        let ay = (a_r * cos_phi - a_phi * sin_phi) * sin_lam + a_lambda * cos_lam;
        let az = a_r * sin_phi + a_phi * cos_phi;

        [ax, ay, az].into()
    }

    fn potential(&self, pos_body_fixed: SVector<f64, 3>) -> Option<f64> {
        let (r, phi, lambda) = Self::spherical(pos_body_fixed);
        let leg = LegendreCache::new(self.coeffs.degree(), phi);

        let mut v = 0.0;
        for n in 0..=self.coeffs.degree() {
            let mut term = 0.0;
            for m in 0..=n {
                let (sin_l, cos_l) = (m as f64 * lambda).sin_cos();
                term += leg.p[n][m]
                    * (self.coeffs.get_c(n, m) * cos_l + self.coeffs.get_s(n, m) * sin_l);
            }
            v += (self.r_ref / r).powi(n as i32) * term;
        }
        Some(self.mu / r * v)
    }
    fn gravity_gradient(&self, pos_body_fixed: SVector<f64, 3>) -> nalgebra::Matrix3<f64> {
        let (r, phi, lambda) = Self::spherical(pos_body_fixed);
        let leg = LegendreCache::new(self.coeffs.degree(), phi);
        let (sin_phi, cos_phi) = phi.sin_cos();
        let (sin_lam, cos_lam) = lambda.sin_cos();
        let tan_phi = sin_phi / cos_phi;
        let sec2_phi = 1.0 / (cos_phi * cos_phi);

        let mut dv_dr = 0.0;
        let mut dv_dphi = 0.0;
        let mut dv_dlambda = 0.0;

        let mut d2v_drdr = 0.0;
        let mut d2v_drdphi = 0.0;
        let mut d2v_drdlambda = 0.0;
        let mut d2v_dphidphi = 0.0;
        let mut d2v_dphidlambda = 0.0;
        let mut d2v_dlambdadlambda = 0.0;

        for n in 0..=self.coeffs.degree() {
            let rn = (self.r_ref / r).powi(n as i32);
            let n_f = n as f64;
            for m in 0..=n {
                let m_f = m as f64;
                let (sin_ml, cos_ml) = (m_f * lambda).sin_cos();
                let cs = self.coeffs.get_c(n, m) * cos_ml + self.coeffs.get_s(n, m) * sin_ml;
                let cs_dl = -self.coeffs.get_c(n, m) * sin_ml + self.coeffs.get_s(n, m) * cos_ml;

                let p = leg.p[n][m];
                let dp = leg.dp[n][m];
                let d2p = tan_phi * dp + m_f * m_f * sec2_phi * p - n_f * (n_f + 1.0) * p;

                dv_dr += (n_f + 1.0) * rn * p * cs;
                dv_dphi += rn * dp * cs;
                dv_dlambda += rn * p * m_f * cs_dl;

                d2v_drdr += (n_f + 1.0) * (n_f + 2.0) * rn * p * cs;
                d2v_drdphi += (n_f + 1.0) * rn * dp * cs;
                d2v_drdlambda += (n_f + 1.0) * rn * p * m_f * cs_dl;
                d2v_dphidphi += rn * d2p * cs;
                d2v_dphidlambda += rn * dp * m_f * cs_dl;
                d2v_dlambdadlambda += rn * p * m_f * m_f * cs;
            }
        }

        let dv_dr = -self.mu / (r * r) * dv_dr;
        let dv_dphi = self.mu / r * dv_dphi;
        let dv_dlambda = self.mu / r * dv_dlambda;

        let d2v_drdr = self.mu / (r * r * r) * d2v_drdr;
        let d2v_drdphi = -self.mu / (r * r) * d2v_drdphi;
        let d2v_drdlambda = -self.mu / (r * r) * d2v_drdlambda;
        let d2v_dphidphi = self.mu / r * d2v_dphidphi;
        let d2v_dphidlambda = self.mu / r * d2v_dphidlambda;
        let d2v_dlambdadlambda = -self.mu / r * d2v_dlambdadlambda;

        // Hessian in the local (r, phi, lambda) orthonormal frame
        let g_rr = d2v_drdr;
        let g_pp = dv_dr / r + d2v_dphidphi / (r * r);
        let g_ll = dv_dr / r - dv_dphi * tan_phi / (r * r)
            + d2v_dlambdadlambda / (r * r * cos_phi * cos_phi);
        let g_rp = d2v_drdphi / r - dv_dphi / (r * r);
        let g_rl = d2v_drdlambda / (r * cos_phi) - dv_dlambda / (r * r * cos_phi);
        let g_pl = (d2v_dphidlambda + dv_dlambda * tan_phi) / (r * r * cos_phi);

        let gamma_local =
            nalgebra::Matrix3::new(g_rr, g_rp, g_rl, g_rp, g_pp, g_pl, g_rl, g_pl, g_ll);

        // exact same basis vectors your ax/ay/az lines build implicitly
        let r_hat: SVector<f64, 3> = [cos_phi * cos_lam, cos_phi * sin_lam, sin_phi].into();
        let phi_hat: SVector<f64, 3> = [-sin_phi * cos_lam, -sin_phi * sin_lam, cos_phi].into();
        let lambda_hat: SVector<f64, 3> = [-sin_lam, cos_lam, 0.0].into();

        let rot = nalgebra::Matrix3::from_columns(&[r_hat, phi_hat, lambda_hat]);

        rot * gamma_local * rot.transpose()
    }
}

#[cfg(test)]
mod gravity_gradient_tests {
    use super::*;
    use nalgebra::{SMatrix, SVector};

    /// Builds a `HarmonicCoeffs` with only C(0,0) = 1, S(0,0) = 0, degree 0.
    /// This makes the spherical-harmonics potential collapse to mu/r, i.e.
    /// a plain point mass.
    fn point_mass_coeffs() -> HarmonicCoeffs {
        let mut coeffs = HarmonicCoeffs::zeros(0);
        coeffs.set(0, 0, 1.0, 0.0);
        coeffs
    }

    fn analytic_point_mass_gradient(mu: f64, pos: SVector<f64, 3>) -> SMatrix<f64, 3, 3> {
        let r = pos.norm();
        let mut gamma = SMatrix::<f64, 3, 3>::identity() * (-mu / r.powi(3));
        gamma += (pos * pos.transpose()) * (3.0 * mu / r.powi(5));
        gamma
    }

    fn assert_matrix_close(a: &SMatrix<f64, 3, 3>, b: &SMatrix<f64, 3, 3>, tol: f64) {
        for i in 0..3 {
            for j in 0..3 {
                let diff = (a[(i, j)] - b[(i, j)]).abs();
                assert!(
                    diff < tol,
                    "mismatch at ({i},{j}): got {}, expected {}, diff {diff}",
                    a[(i, j)],
                    b[(i, j)]
                );
            }
        }
    }

    #[test]
    fn spherical_harmonics_point_mass_matches_closed_form() {
        let mu = 3.986e14; // m^3/s^2, Earth-like
        let r_ref = 6.378e6; // m
        let field = SphericalHarmonics::new(mu, r_ref, point_mass_coeffs());
        let pos = SVector::<f64, 3>::new(6_000_000.0, 2_000_000.0, 1_500_000.0);

        let expected = analytic_point_mass_gradient(mu, pos);
        let actual = field.gravity_gradient(pos);

        assert_matrix_close(&actual, &expected, 1e-12);
    }

    #[test]
    fn spherical_harmonics_gradient_is_symmetric_and_traceless() {
        let mu = 3.986e14;
        let r_ref = 6.378e6;
        let field = SphericalHarmonics::new(mu, r_ref, point_mass_coeffs());
        let pos = SVector::<f64, 3>::new(7_000_000.0, -500_000.0, 3_200_000.0);

        let gamma = field.gravity_gradient(pos);

        assert_matrix_close(&gamma, &gamma.transpose(), 1e-12);
        let trace = gamma[(0, 0)] + gamma[(1, 1)] + gamma[(2, 2)];
        assert!(trace.abs() < 1e-12, "trace should be ~0, got {trace}");
    }
}

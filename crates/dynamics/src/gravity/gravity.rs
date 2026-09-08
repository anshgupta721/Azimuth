use nalgebra::SVector;

use crate::gravity::spherical_harmonics::{HarmonicCoeffs, LegendreCache};
pub trait GravityField {
    fn acceleration(&mut self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3>;
    /// Optional: f
    fn potential(&mut self, pos_body_fixed: SVector<f64, 3>) -> Option<f64>;
}

pub struct PointMass {
    pub mu: f64,
}

impl GravityField for PointMass {
    fn acceleration(&mut self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3> {
        let r = pos_body_fixed.norm();
        let common = -self.mu / r.powf(3.0);
        pos_body_fixed * common
    }

    fn potential(&mut self, pos_body_fixed: SVector<f64, 3>) -> Option<f64> {
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
    fn spherical(pos: SVector<f64, 3>) -> (f64, f64, f64) {
        let r = pos.norm();
        let phi = (pos[2] / r).asin();
        let lambda = pos[1].atan2(pos[0]);
        (r, phi, lambda)
    }
}

impl GravityField for SphericalHarmonics {
    fn acceleration(&mut self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3> {
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

    fn potential(&mut self, pos_body_fixed: SVector<f64, 3>) -> Option<f64> {
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
}

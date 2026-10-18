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
    degree: usize,
}

impl SphericalHarmonics {
    pub fn new(mu: f64, r_ref: f64, degree: usize) -> Self {
        SphericalHarmonics { mu, r_ref, degree }
    }

    fn spherical(pos: SVector<f64, 3>) -> (f64, f64, f64) {
        let r = pos.norm();
        let phi = (pos[2] / r).asin();
        let lambda = pos[1].atan2(pos[0]);
        (r, phi, lambda)
    }

    pub fn acceleration_and_gradient(
        &self,
        pos_body_fixed: SVector<f64, 3>,
        coeffs: HarmonicCoeffs
    ) -> (SVector<f64, 3>, SMatrix<f64, 3, 3>) {
        let (r, phi, lambda) = Self::spherical(pos_body_fixed);
        let leg = LegendreCache::new(self.degree, phi);
        let (sin_phi, cos_phi) = phi.sin_cos();
        let (sin_lam, cos_lam) = lambda.sin_cos();

        let mut dv_dr = 0.0;
        let mut dv_dphi = 0.0;
        let mut dv_dlambda = 0.0;

        let mut d2v_drr = 0.0;
        let mut d2v_drdphi = 0.0;
        let mut d2v_drdlambda = 0.0;
        let mut d2v_dphiphi = 0.0;
        let mut d2v_dphidlambda = 0.0;
        let mut d2v_dlambdalambda = 0.0;

        for n in 0..=coeffs.degree() {
            let rn = (self.r_ref / r).powi(n as i32);
            for m in 0..=n {
                let (sin_ml, cos_ml) = (m as f64 * lambda).sin_cos();
                let cs = coeffs.get_c(n, m) * cos_ml + coeffs.get_s(n, m) * sin_ml;
                let cs_dl = -coeffs.get_c(n, m) * sin_ml + coeffs.get_s(n, m) * cos_ml;

                dv_dr += (n as f64 + 1.0) * rn * leg.p[n][m] * cs;
                dv_dphi += rn * leg.dp[n][m] * cs;
                dv_dlambda += rn * leg.p[n][m] * (m as f64) * cs_dl;

                let (nf, mf) = (n as f64, m as f64);

                d2v_drr += rn * (nf + 2.0) * (nf + 1.0) * cs * leg.p[n][m];
                d2v_drdphi += rn * (nf + 1.0) * cs * leg.dp[n][m];
                d2v_drdlambda += rn * (nf + 1.0) * mf * cs_dl * leg.p[n][m];
                d2v_dphiphi += rn * cs * leg.d2p[n][m];
                d2v_dphidlambda += rn * mf * cs_dl * leg.dp[n][m];
                d2v_dlambdalambda += rn * mf * mf * cs * leg.p[n][m];
            }
        }

        let dv_dr = -self.mu / (r * r) * dv_dr;
        let dv_dphi = self.mu / r * dv_dphi;
        let dv_dlambda = self.mu / r * dv_dlambda;

        let d2v_drr = self.mu / (r * r * r) * d2v_drr;
        let d2v_drdphi = -self.mu / (r * r) * d2v_drdphi;
        let d2v_drdlambda = -self.mu / (r * r) * d2v_drdlambda;
        let d2v_dphiphi = self.mu / r * d2v_dphiphi;
        let d2v_dphidlambda = self.mu / r * d2v_dphidlambda;
        let d2v_dlambdalambda = -self.mu / r * d2v_dlambdalambda;

        // Spherical -> local components
        let a_r = dv_dr;
        let a_phi = dv_dphi / r;

        // Note: singular at the poles (cos_phi -> 0)
        let a_lambda = dv_dlambda / (r * cos_phi);

        let tan_phi = sin_phi / cos_phi;
        let g_rr = d2v_drr;
        let g_rp = d2v_drdphi / r - dv_dphi / (r * r);
        let g_rl = d2v_drdlambda / (r * cos_phi) - dv_dlambda / (r * r * cos_phi);
        let g_pp = d2v_dphiphi / (r * r) + dv_dr / r;
        let g_pl = d2v_dphidlambda / (r * r * cos_phi) + dv_dlambda * tan_phi / (r * r * cos_phi); //Potentially a single cos instead of squared
        let g_ll = d2v_dlambdalambda / (r * r * cos_phi * cos_phi) - tan_phi * dv_dphi / (r * r)
            + dv_dr / r;

        let gamma = SMatrix::<f64, 3, 3>::new(
            g_rr, g_rp, g_rl, //
            g_rp, g_pp, g_pl, //
            g_rl, g_pl, g_ll,
        );

        let e_r = SVector::<f64, 3>::new(cos_phi * cos_lam, cos_phi * sin_lam, sin_phi);
        let e_phi = SVector::<f64, 3>::new(-sin_phi * cos_lam, -sin_phi * sin_lam, cos_phi);
        let e_lam = SVector::<f64, 3>::new(-sin_lam, cos_lam, 0.0);

        let ax = (a_r * cos_phi - a_phi * sin_phi) * cos_lam - a_lambda * sin_lam;
        let ay = (a_r * cos_phi - a_phi * sin_phi) * sin_lam + a_lambda * cos_lam;
        let az = a_r * sin_phi + a_phi * cos_phi;

        let rot = SMatrix::<f64, 3, 3>::from_columns(&[e_r, e_phi, e_lam]);

        let gradient = rot * gamma * rot.transpose();

        ([ax, ay, az].into(), gradient)
    }
}


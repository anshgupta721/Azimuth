use nalgebra::SVector;

pub trait GravityField {
    fn acceleration(&self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3>;
    /// Optional: f
    fn potential(&self, _pos_body_fixed: SVector<f64, 3>) -> Option<f64> {
        None
    }
}

pub struct PointMass {
    mu: f64,
}

impl GravityField for PointMass {
    fn acceleration(&self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3> {
        let r = pos_body_fixed.norm();
        let common = -self.mu / r.powf(3.0);
        pos_body_fixed * common
    }
}

// pub struct SphericalHarmonics{
//     mu: f64,
//     r_ref: f64,
//     pub coeffs: HarmonicCoeffs,
// }

// impl GravityField for SphericalHarmonics{
//     fn acceleration(&self, pos_body_fixed: SVector<f64, 3>) -> SVector<f64, 3> {

//         pos_body_fixed
//     }

//     // pub fn potential(&self, _pos_body_fixed: SVector<f64, 3>) -> Option<f64> {
//     //     None
//     // }
// }

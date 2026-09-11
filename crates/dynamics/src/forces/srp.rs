use nalgebra::SVector;
use std::f64::consts::PI;

pub fn srp_acceleration(
    r_sc: SVector<f64, 3>,
    r_sun: SVector<f64, 3>,
    area_mass_ratio: f64,
    reflectivity: f64,
) -> SVector<f64, 3> {
    let r_sun_to_sc = r_sc - r_sun;
    let d_sun_to_sc = r_sun_to_sc.norm();
    let u_sun_to_sc = r_sun_to_sc / d_sun_to_sc;

    let i_s = 3.8e26 / (4.0 * PI * d_sun_to_sc * d_sun_to_sc);

    let p = i_s / 299792458.0;

    p * area_mass_ratio * reflectivity * u_sun_to_sc
}

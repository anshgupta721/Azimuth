use nalgebra::SVector;

pub fn rk4<const NX: usize, const NU: usize>(
    f: &dyn Fn(f64, SVector<f64, NX>, SVector<f64, NU>) -> SVector<f64, NX>,
    t: f64,
    dt: f64,
    x: SVector<f64, NX>,
    u: SVector<f64, NU>,
) -> SVector<f64, NX> {
    let k1 = dt * f(t, x, u);
    let k2 = dt * f(t + dt / 2.0, x + k1 / 2.0, u);
    let k3 = dt * f(t + dt / 2.0, x + k2 / 2.0, u);
    let k4 = dt * f(t + dt, x + k3, u);
    x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
}

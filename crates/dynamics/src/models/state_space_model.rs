use nalgebra::{SMatrix, SVector};

pub struct StateSpace<const NX: usize, const NU: usize, const NY: usize> {
    a: SMatrix<f64, NX, NX>,
    b: SMatrix<f64, NX, NU>,
    c: SMatrix<f64, NY, NX>,
    d: SMatrix<f64, NY, NU>,
}

impl<const NX: usize, const NU: usize, const NY: usize> StateSpace<NX, NU, NY> {
    pub fn new(
        a: SMatrix<f64, NX, NX>,
        b: SMatrix<f64, NX, NU>,
        c: SMatrix<f64, NY, NX>,
        d: SMatrix<f64, NY, NU>,
    ) -> StateSpace<NX, NU, NY> {
        StateSpace { a, b, c, d }
    }

    pub fn derivative(&self, x: SVector<f64, NX>, u: SVector<f64, NU>) -> SVector<f64, NX> {
        return self.a * x + self.b * u;
    }

    pub fn output(&self, x: SVector<f64, NX>, u: SVector<f64, NU>) -> SVector<f64, NY> {
        return self.c * x + self.d * u;
    }
}

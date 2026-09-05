use nalgebra::{SMatrix, SVector};

pub trait LTVSystem<const NX: usize, const NU: usize, const NY: usize> {
    fn a(&self, t: f64) -> SMatrix<f64, NX, NX>;
    fn b(&self, t: f64) -> SMatrix<f64, NX, NU>;
    fn c(&self, t: f64) -> SMatrix<f64, NY, NX>;
    fn d(&self, t: f64) -> SMatrix<f64, NY, NU>;
}

// pub trait LTISystem<const NX: usize, const NU: usize, const NY: usize

pub struct StateSpace<S, const NX: usize, const NU: usize, const NY: usize> {
    system: S,
}

impl<S, const NX: usize, const NU: usize, const NY: usize> StateSpace<S, NX, NU, NY>
where
    S: LTVSystem<NX, NU, NY>,
{
    pub fn new(system: S) -> StateSpace<S, NX, NU, NY> {
        StateSpace { system }
    }

    pub fn derivative(&self, t: f64, x: SVector<f64, NX>, u: SVector<f64, NU>) -> SVector<f64, NX> {
        self.system.a(t) * x + self.system.b(t) * u
    }

    pub fn output(&self, t: f64, x: SVector<f64, NX>, u: SVector<f64, NU>) -> SVector<f64, NY> {
        return self.system.c(t) * x + self.system.d(t) * u;
    }
}

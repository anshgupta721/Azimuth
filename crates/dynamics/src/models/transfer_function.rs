use super::state_space_model::StateSpace;
use nalgebra::SMatrix;

pub struct TransferFunction<const N: usize> {
    num: [f64; N],
    den: [f64; N],
}

impl<const N: usize> TransferFunction<N> {
    pub fn new(num: [f64; N], den: [f64; N]) -> TransferFunction<N> {
        TransferFunction { num, den }
    }

    pub fn to_ss<const NX: usize, const NU: usize, const NY: usize>() -> StateSpace<NX, NU, NY> {
        // Turn this jawn into a state space model
        // StateSpace::<NX,NU,NY>::new(_,_,_,_)
        todo!()
    }
}

pub struct MimoTransferFunction<const INPUTS: usize, const OUTPUTS: usize, const N: usize> {
    mimo_tf: SMatrix<TransferFunction<N>, OUTPUTS, INPUTS>,
}

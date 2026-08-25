use nalgebra::SMatrix;
use super::state_space_model::StateSpace;

pub struct TransferFunction<const MAX_ORDER: usize> {
    num: [f64; MAX_ORDER + 1],
    den: [f64; MAX_ORDER + 1],
}

impl<const MAX_ORDER: usize> TransferFunction<MAX_ORDER> {
    pub fn new(num: [f64; MAX_ORDER + 1], den: [f64; MAX_ORDER + 1]) -> TransferFunction<MAX_ORDER> {
        TransferFunction { num, den }
    }

    pub fn to_ss() -> StateSpace {
        // Turn this jawn into a state space model
    }
}

pub struct MimoTransferFunction<const INPUTS: usize, const OUTPUTS: usize, const MAX_ORDER: usize> {
    mimo_tf: SMatrix<TransferFunction<MAX_ORDER>, OUTPUTS, INPUTS>,
}


pub trait EquationOfState {
    fn primitive_to_conserved(&mut self);
    fn conserved_to_primitive(&mut self);
}

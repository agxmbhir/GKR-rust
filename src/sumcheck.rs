use ark_poly::{Evaluations};
use ark_ff::FftField;

// How do you represent a 
struct Prover {
    // The polynomial to be verified
    f: Evaluations<FftField>, // This needs to be multilinear // TODO 
    // The domain of the polynomial
    domain: EvaluationDomain<FftField>,
    // The current degree
    degree: usize,
}

impl Prover { 
    // Create a new prover
    pub fn new(poly: Evaluations<FftField>, domain: EvaluationDomain<FftField>) -> Self {
        Self {
            f: poly,
            domain,
            degree: poly.len() - 1,
        }
    }  

    // At each step, prover computes g_i where g(x) = Summation of f(x, a_i) 
    // for all a_i in domain D which is a boolean hypercube
    // boolean hypercube is a set of all possible binary strings of length n
    // For example, g_0 = f(x, 0, 0) + f(x, 0, 1) + f(x, 1, 0) + f(x, 1, 1)
    pub fn step(&mut self, step: usize) -> Evaluations<FftField> { 
        // Cool beans
        
    }

    
}
pub fn verify(sum: Evaluations<FftField>, domain: EvaluationDomain<FftField>) -> bool {
    // TODO
    // 1. Create a prover
    // 2. Create a verifier
    // 3. Run the protocol
    // 4. Return the result
    true
}
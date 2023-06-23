use std::ops::Add;

use ark_poly::{ multivariate::{SparsePolynomial, Term}, univariate::SparsePolynomial as UniSparsePolynomial};
use ark_ff::{Field, Zero};

// 
struct Prover<F: Field, T: Term> {
    // The polynomial to be verified
    f: SparsePolynomial<F, T>, // This needs to be multilinear // TODO 
    // The domain of the polynomial
   // domain: EvaluationDomain<FftField>, Do we need this? 
}

impl<F: Field, T: Term>  Prover<F, T> { 
    // Create a new prover
    pub fn new(poly: SparsePolynomial<F, T>) -> Self {
        Self {
            f: poly.clone(),
        }
    }  

    // At each step, prover computes g_i where g(x) = Summation of f(x, a_i) 
    // for all a_i in domain D which is a boolean hypercube
    // boolean hypercube is a set of all possible binary strings of length n
    // For example, g_0 = f(x, 0, 0) + f(x, 0, 1) + f(x, 1, 0) + f(x, 1, 1)
    pub fn get_univar_poly(&mut self, step: usize) -> UniSparsePolynomial<F>{ 
        // Cool beans.
        let params = get_permutations(self.f.num_vars);
        // Constructing a sparse polynomial
        // evaluating each sparse monomial at point (X, ... a_i) where a_i is a permutation of 0s and 1s
        // Example: f(x, y, z) = 2x2y + 3xy2z^2 + 4z^3
        // We construct a vector of monomials: [2x2y, 3xy2z^2, 4z^3]
        // Evaluate each monomial at point (X, 0, 0), (X, 0, 1), (X, 1, 0), (X, 1, 1)
        // We get: [2X2, 3X, 4X]
        // We then interpolate the above vector to get g_0
        // We repeat the above process for g_1, g_2, ... g_n
        // We then interpolate the above vector to get g 
        let mut g = UniSparsePolynomial::<F>::zero();
        let f = self.f.clone();
        // Vars = (degree, index)
        for param in params { 
            // Evaluation of f at point (X, ... a_i)
            let mut f_i: Vec<(usize, F)> = Vec::<(usize, F)>::new();
            // Going through each term in f, and evaluating it at point (..r_i, X, ... a_i)
            for term in f.terms.iter() {
                let (coeff, vars) = term.clone();
                let mut term_eval = coeff.clone();
                let mut prod = F::one();
                let mut contains_fixed_var = false;
                let mut fixed_var_degree = 0;
                // for each var in term, we multiply it by a_i except for the one that is X
                if vars.is_empty() {
                    f_i.push((0, coeff));
                    break;
                }
                println!("vars: {:?}", vars);
                for var in vars.iter()  {
                    let (var_index, var_degree) = var.clone();
                    if var_index == step {
                        contains_fixed_var = true;
                        fixed_var_degree = var_degree;
                    } else {
                        prod *= if param[var_index] {F::one()} else {F::zero()};
                    }
                }
                term_eval *= prod;
                if term_eval != F::zero() {
                    if contains_fixed_var {
                        f_i.push((fixed_var_degree, term_eval));
                    } else {
                        f_i.push((0, term_eval));
                    }
                }
                
            } 
            g = g.add(UniSparsePolynomial::<F>::from_coefficients_vec(f_i));
        }
        return g;
    }
} 


// Returns a vector of all possible permutations of a vector of length n
fn get_permutations(n: usize) -> Vec<Vec<bool>> {
    let mut permutations = Vec::new();
    for i in 0..2usize.pow(n as u32) {
        let mut permutation = Vec::new();
        for j in 0..n {
            permutation.push((i >> j) & 1 == 1);
        }
        permutations.push(permutation);
    }
    permutations
}

// Returns the variables of a polynomial
fn get_vars<F: Field,T: Term> (f: SparsePolynomial<F, T>) -> Vec<(usize, usize)>{
    let mut vars = vec![];
    for term in f.terms {
        let (_, term) = term.clone();
        for var in term.iter() {
            if vars.contains(var) {
                continue;
            } else {
                vars.push(var.clone());
            }
        }
    }
    vars
}
 
pub fn verify() -> bool {
    // TODO
    // 1. Create a prover
    // 2. Create a verifier
    // 3. Run the protocol
    // 4. Return the result
    true
} 

#[cfg(test)]
mod tests {
    use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term},
     DenseMVPolynomial,
    univariate::SparsePolynomial as UniSparsePolynomial
    };
    use ark_test_curves::fp128::Fq;
    use super::Prover;

    #[test]
    fn test_get_univar_poly() {
        let poly = SparsePolynomial::from_coefficients_vec(
            3,
            // 2*x_0^3 + x_0*x_2 + x_1*x_2 + 5
            vec![
                (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
                (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
                (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
                (Fq::from(5), SparseTerm::new(vec![])),
            ]
        );
        let mut p = Prover::new(
            poly
        );

        let uni_1 = p.get_univar_poly(0);
        let exp = 
        UniSparsePolynomial::from_coefficients_vec(
            // 8x^3 + 2x + 21
            vec![
                (3, Fq::from(8)),
                (1, Fq::from(2)),
                (0, Fq::from(21)),
               ]
        );
        assert_eq!(exp, uni_1);
    
    }
}
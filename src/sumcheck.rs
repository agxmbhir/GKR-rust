use std::ops::Add;
use ark_poly::{ multivariate::{SparsePolynomial, Term}, univariate::SparsePolynomial as UniSparsePolynomial, Polynomial};
use ark_ff::{Field, Zero};
use rand::Rng;

struct Prover<F: Field, T: Term> {
    f: SparsePolynomial<F, T>,  
}

impl<F: Field, T: Term>  Prover<F, T> { 
    // Create a new prover
    pub fn new(poly: SparsePolynomial<F, T>) -> Self {
        Self {
            f: poly.clone(),
        }
    }  

    // Returns the univariate polynomial g_i(x_i) for a given index i
    pub fn get_univar_poly(&mut self, step: usize) -> UniSparsePolynomial<F>{ 

        let params = get_permutations(self.f.num_vars - 1, step); 
        let mut g = Vec::<(usize, F)>::new();
        let f = self.f.clone();

        for param in params { 
            let mut f_i = Vec::new();
            // Each term is a monomial in the polynomial
            for term in f.terms.iter() {
                // Getting the vars in the term
                let (coeff, vars) = term.clone();
                let mut term_eval = coeff.clone();

                let mut prod = F::one();
                let mut contains_fixed_var = false;
                let mut fixed_var_degree = 0;

                // If the term is a constant, we don't need to do anything
                if vars.is_empty() {
                    f_i.push((0, coeff));
                } 

                else { 
                    // we iterate over the variables 
                  for var in vars.iter()  {
                    let (var_index, var_degree) = var.clone();
                    // We skip the one we want to keep and create the uni-variate polynomial
                    if var_index == step {
                        contains_fixed_var = true;
                        fixed_var_degree = var_degree;
                    } else {
                        // evaluting the term for the current permutation
                        prod *= if param[var_index] {F::one()} else {F::zero()};
                    }
                }
                // completing the evaluation of the term by multiplying it with the coeff
                term_eval *= prod;
                if term_eval != F::zero() {
                    if contains_fixed_var {
                        f_i.push((fixed_var_degree, term_eval));
                    } else {
                        f_i.push((0, term_eval));
                    }
                } 
            }
        } 
            // We construct the uni-variate polynomial by adding the terms with the same degree
            f_i.iter().for_each(|(degree, coeff)| {
                let mut found = false;
                for (i, (d, c)) in g.iter_mut().enumerate() {
                    if *d == *degree {
                        *c += coeff.clone();
                        found = true;
                        break;
                    }
                }
                if !found {
                    g.push((*degree, coeff.clone()));
                }
            });
        }

        return UniSparsePolynomial::<F>::from_coefficients_vec(g);
    }
   

} 



 // Returns a vector of all possible permutations of a vector of length n + a placeholder for the fixed index
fn get_permutations(n: usize, fixed_index: usize) -> Vec<Vec<bool>> {
    let mut permutations = Vec::new();
    for i in 0..2usize.pow(n as u32) {
        let mut permutation = Vec::new();
        for j in 0..n  {
              if j == fixed_index {
                permutation.push(true);
              }
              permutation.push((i >> j) & 1 == 1);
        }
        permutations.push(permutation);
    }
    permutations
}

fn get_random() ->  Option<Field>{
	let mut rng = rand::thread_rng();
	let r: Field = rng.gen();
	Some(r)
}
 
pub fn verify<F: Field, T: Term>(poly: SparsePolynomial<F, T>, claim: F) -> bool {
    	// 1st round
	let mut p = Prover::new(poly);
	let mut gi = p.get_univar_poly(0);
	let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
    assert_eq!(expected_c, claim);


    true
} 

#[cfg(test)]
mod tests {
    use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term},
     DenseMVPolynomial,
    univariate::SparsePolynomial as UniSparsePolynomial,
    Polynomial
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
        let uni_2 = p.get_univar_poly(1);
        let uni_3 = p.get_univar_poly(2);

        let exp1 = 
        UniSparsePolynomial::from_coefficients_vec(
            // 8x^3 + 2x + 21
            vec![
                (3, Fq::from(8)),
                (1, Fq::from(2)),
                (0, Fq::from(21)),
               ]
        );

        let exp2 = 
        UniSparsePolynomial::from_coefficients_vec(
             // 2x + 25
            vec![
                (1, Fq::from(2)),
                (0, Fq::from(25)),
               ]
        );
        let exp3 =
        UniSparsePolynomial::from_coefficients_vec(
            // 4x + 24
           vec![
               (1, Fq::from(4)),
               (0, Fq::from(24)),
              ]
        );  

        assert_eq!(exp1, uni_1);
        assert_eq!(exp2, uni_2);
        assert_eq!(exp3, uni_3);
    
    }

    #[test]
    fn test_verifier_first_step() {
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
        let claimed_sum = Fq::from(52);
        
        let mut p = Prover::new(
            poly
        );
        let uni_1 = p.get_univar_poly(0);
        let eval = uni_1.evaluate(&0u32.into()) + uni_1.evaluate(&1u32.into());
        assert_eq!(eval, claimed_sum);
    }
}
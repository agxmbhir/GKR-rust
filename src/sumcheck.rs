
use ark_ff::{Field, Zero};
use ark_poly::{
    multivariate::{SparsePolynomial, Term, SparseTerm},
    univariate::{SparsePolynomial as UniSparsePolynomial},
    Polynomial,
};
use ark_bls12_381::{Fr, fr};
use rand::distributions::{Distribution, Standard};
use rand::Rng;

struct Prover<T: Term> {
    f: SparsePolynomial<Fr, T>,
    steps: Vec<Fr>,
}

impl<T: Term> Prover<T> {
    // Create a new prover
    pub fn new(poly: SparsePolynomial<Fr, T>) -> Self {
        Self {
            f: poly.clone(),
            steps: Vec::new(),
        }
    }

    pub fn gen_univar_poly(&mut self, random: Option<Fr>) -> UniSparsePolynomial<Fr> {

        if random.is_some() {
            self.steps.push(random.unwrap());
         }

        let params = get_permutations(self.f.num_vars - self.steps.len() - 1);

        let mut poly:UniSparsePolynomial<Fr> = UniSparsePolynomial::from_coefficients_vec(vec![(0, Fr::from(0))]);

        for param in params {
           poly = poly + self.evaluate_at_point(&param);
        }
        poly
    }

    // loops through all the terms in the polynomial and evaluates them except the "fixed" variable
    fn evaluate_at_point(&self, param: &Vec<Fr>) -> UniSparsePolynomial<Fr> {

        let mut uni_terms: Vec<(usize, Fr)> = Vec::new();
        for term in self.f.terms.clone().into_iter() {

            let (coeff, vars) = term;
            let prod = coeff;
            let mut term_eval = (0, prod);

            for (var, power) in vars.iter() {
                match *var {
                    i if i > self.steps.len() =>  { 
                        let p = param[*var - self.steps.len() - 1];
                        if p == Fr::from(0) {
                            term_eval = (0, Fr::from(0));
                        }
                    }
                    i if i < self.steps.len() => {
                        term_eval.1 = term_eval.1 * self.steps[i].pow([*power as u64]);
                    },
                   _ => {  
                        term_eval.0 = *power;
                    },
                  
                }
            }
            uni_terms.push(
                term_eval
            );
            // Runs in O(n^2) time :( // TODO: optimize // Sums up all the terms in uni_terms with the same power
            // for i in 0..uni_terms.len() {
            //     for j in i+1..uni_terms.len() {
            //         if uni_terms[i].0 == uni_terms[j].0 {
            //             uni_terms[i].1 = uni_terms[i].1 + uni_terms[j].1;
            //             uni_terms[j].1 = Fr::from(0);
            //         }
            //     }
            // }
        }
        let poly = UniSparsePolynomial::from_coefficients_vec(uni_terms);
        poly
    }

    fn evaluate(&self, steps: &Vec<Fr>) -> Fr {
        let mut poly = self.f.clone();
        let mut eval = Fr::zero();
        for term in poly.terms.iter_mut() {
            let (coeff, vars) = term;
            let mut prod = coeff.clone();
            for (var, power) in vars.iter() {
                prod *= steps[*var - 1].pow([*power as u64]);
            }
            eval += prod;
        }
        return eval;
    }

    }


// Gives the permuations of a boolean vector
fn get_permutations<>(n: usize) -> Vec<Vec<Fr>> {
    let mut permutations = Vec::new();
    for i in 0..2usize.pow(n as u32) {
        let mut permutation = Vec::new();
        for j in 0..n {
            if (i >> j) & 1 == 1 {
                permutation.push(Fr::from(1));
            } else {
                permutation.push(Fr::from(0));
            }
        }
        permutations.push(permutation);
    }
    permutations
}

pub fn get_random() -> Option<Fr> {
	let mut rng = rand::thread_rng();
	let r: Fr = rng.gen();
	Some(r)
}

// pub fn verify<F: Field, T: Term>(poly: SparsePolynomial<F, T>, claim: F) -> bool
// where
//     Standard: Distribution<F>,
// {
//     // // 1st round
//     // let mut p = Prover::new(poly);
//     // let mut gi = p.gen_univar_poly(0, get_random());
//     // let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
//     // assert_eq!(expected_c, claim);
//     // // middle rounds
//     // for j in 1..p.f.num_vars {
//     //     let r = get_random();
//     //     expected_c = gi.evaluate(&r);
//     //     gi = p.get_univar_poly(j, r);
//     //     let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
//     //     assert_eq!(expected_c, new_c);
//     // }
//     // // final round
//     // let r = get_random();
//     // expected_c = gi.evaluate(&r);

//     // assert_eq!(expected_c, new_c);

//     true
// }

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_poly::{
        multivariate::{SparsePolynomial, SparseTerm, Term},
        DenseMVPolynomial, Polynomial,
    };

    use crate::sumcheck::get_random;
    use super::Prover;


    #[test]
    fn test_verifier() {
        // 1st round
        let poly = SparsePolynomial::from_coefficients_vec(
            3,
            // 2*x_0^3 + x_0*x_2 + x_1*x_2 + 5
            vec![
                (Fr::from(2), SparseTerm::new(vec![(0, 3)])),
                (Fr::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
                (Fr::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
                (Fr::from(5), SparseTerm::new(vec![])),
            ],
        );

        let mut p = Prover::new(poly);
        let mut gi = p.gen_univar_poly(None);
        let claim = Fr::from(52);
        let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
        assert_eq!(expected_c, claim);

        // middle rounds
        for j in 1..p.f.num_vars {
            let r = get_random();
            expected_c = gi.evaluate(&r.unwrap());
            gi = p.gen_univar_poly(r);
            let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
            assert_eq!(expected_c, new_c);
        }


        // final round
        let r = get_random();
        expected_c = gi.evaluate(&r.unwrap());
        println!("len {}", &p.steps.len());
        p.steps.push(r.unwrap());
        let new_c = p.f.evaluate(&p.steps);
        println!("expected_c: {:?}, new_c: {:?}", expected_c, new_c);
        assert_eq!(expected_c, new_c);
    }
}

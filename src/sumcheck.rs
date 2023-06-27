
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
           println!("params {:?}", param);
           poly = poly + self.evaluatePoly(&param);
        }
        println!("poly: {:?}", poly);
        poly
    }

    fn evaluatePoly(&self, param: &Vec<Fr>) -> UniSparsePolynomial<Fr> {
        let mut uni_terms: Vec<(usize, Fr)> = Vec::new();
        for term in self.f.terms.clone().into_iter() {
            let (mut coeff, vars) = term;
            let mut prod = coeff;
            for (var, power) in vars.iter() {
                match *var {
                    i if i == self.steps.len() => {
                        uni_terms.push(
                            (*power, prod),
                        )
                    },
                    i if i < self.steps.len() =>
                        uni_terms.push(
                            (0, prod * self.steps[i].pow([*power as u64])),
                        ),
                    _ => uni_terms.push(
                        (0, prod * param[*var - self.steps.len() - 1].pow([*power as u64])),
                    ),
                }
            }
        }
        let poly = UniSparsePolynomial::from_coefficients_vec(uni_terms);
        println!("poly: {:?}", poly);
        poly
    }

    }


// Gives the permuations of a boolean vector
fn get_permutations<>(n: usize) -> Vec<Vec<Fr>> {
    let mut permutations = Vec::new();
    let F0 = Fr::from(0);
    let F1 = Fr::from(1);
    for i in 0..2usize.pow(n as u32) {
        let mut permutation = Vec::new();
        for j in 0..n {
            if (i >> j) & 1 == 1 {
                permutation.push(F1);
            } else {
                permutation.push(F0);
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

pub fn verify<F: Field, T: Term>(poly: SparsePolynomial<F, T>, claim: F) -> bool
where
    Standard: Distribution<F>,
{
    // // 1st round
    // let mut p = Prover::new(poly);
    // let mut gi = p.gen_univar_poly(0, get_random());
    // let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
    // assert_eq!(expected_c, claim);
    // // middle rounds
    // for j in 1..p.f.num_vars {
    //     let r = get_random();
    //     expected_c = gi.evaluate(&r);
    //     gi = p.get_univar_poly(j, r);
    //     let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
    //     assert_eq!(expected_c, new_c);
    // }
    // // final round
    // let r = get_random();
    // expected_c = gi.evaluate(&r);

    // assert_eq!(expected_c, new_c);

    true
}

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
            println!("new_c: {}", new_c);
            println!("expected_c: {}", expected_c);
        }
        // final round
        let r = get_random();
        expected_c = gi.evaluate(&r.unwrap());
        let new_c = p.f.evaluate(&p.steps);
        // assert_eq!(expected_c, new_c);
    }
}

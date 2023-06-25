use ark_ff::Field;
use ark_poly::{
    multivariate::{SparsePolynomial, Term},
    univariate::SparsePolynomial as UniSparsePolynomial,
    Polynomial,
};
use rand::distributions::{Distribution, Standard};
use rand::Rng;

struct Prover<F: Field, T: Term> {
    f: SparsePolynomial<F, T>,
    steps: Vec<F>,
}

impl<F: Field, T: Term> Prover<F, T> {
    // Create a new prover
    pub fn new(poly: SparsePolynomial<F, T>) -> Self {
        Self {
            f: poly.clone(),
            steps: Vec::new(),
        }
    }

    // Returns the univariate polynomial g_i(x_i) for a given index i
    // random is the verifier's random value
    pub fn get_univar_poly(&mut self, step: usize, random: F) -> UniSparsePolynomial<F> {
        self.steps.push(random);
        let params = self.get_params(self.f.num_vars - self.steps.len(), random);

        let mut g = Vec::<(usize, F)>::new();
        let f = self.f.clone();
        println!("params len: {}", params.len());
        println!("paras {:?} ", params);
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
                } else {
                    // we iterate over the variables
                    for var in vars.iter() {
                        let (var_index, var_degree) = var.clone();
                        // We skip the one we want to keep and create the uni-variate polynomial
                        if var_index == step {
                            contains_fixed_var = true;
                            fixed_var_degree = var_degree;
                        } else {
                            // evaluting the term for the current permutation
                            prod *= param[var_index].pow(&[var_degree as u64]);
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
            println!("g: {:?}", g);
        }
        return UniSparsePolynomial::<F>::from_coefficients_vec(g);
    }

    fn get_params(&self, n: usize, r: F) -> Vec<Vec<F>> {
        let mut params = Vec::new();
        let mut perms = get_permutations(n);
        perms.iter().for_each(|perm| {
            let mut param = self.steps.clone();
            param.extend_from_slice(perm);
            params.push(param);
        });
        return params;
    }
}

// Gives the permuations of a boolean vector
fn get_permutations<F: Field>(n: usize) -> Vec<Vec<F>> {
    let mut permutations = Vec::new();
    let F0 = F::zero();
    let F1 = F::one();
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

// A function to get a random value from the field
fn get_random<F: Field>() -> F
where
    Standard: Distribution<F>,
{
    let mut rng = rand::thread_rng();
    let distribution = Standard;
    let random = rng.sample(distribution);
    random
}

pub fn verify<F: Field, T: Term>(poly: SparsePolynomial<F, T>, claim: F) -> bool
where
    Standard: Distribution<F>,
{
    // 1st round
    let mut p = Prover::new(poly);
    let mut gi = p.get_univar_poly(0, get_random());
    let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
    assert_eq!(expected_c, claim);
    // middle rounds
    for j in 1..p.f.num_vars {
        let r = get_random();
        expected_c = gi.evaluate(&r);
        gi = p.get_univar_poly(j, r);
        let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
        assert_eq!(expected_c, new_c);
    }
    // final round
    let r = get_random();
    expected_c = gi.evaluate(&r);
    let new_c = p.f.evaluate(&p.steps);
    assert_eq!(expected_c, new_c);

    true
}

#[cfg(test)]
mod tests {
    use ark_poly::{
        multivariate::{SparsePolynomial, SparseTerm, Term},
        univariate::SparsePolynomial as UniSparsePolynomial,
        DenseMVPolynomial, Polynomial,
    };

    use crate::sumcheck::get_random;
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
            ],
        );
        let mut p = Prover::new(poly);
        let uni_0 = p.get_univar_poly(0, get_random());
        let uni_1 = p.get_univar_poly(1, get_random());
        let uni_2 = p.get_univar_poly(2, get_random());

        let exp1 = UniSparsePolynomial::from_coefficients_vec(
            // 8x^3 + 2x + 21
            vec![(3, Fq::from(8)), (1, Fq::from(2)), (0, Fq::from(21))],
        );

        assert_eq!(exp1, uni_1);
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
            ],
        );
        let claimed_sum = Fq::from(52);

        let mut p = Prover::new(poly);
        let uni_1 = p.get_univar_poly(0, get_random());
        println!("uni_1: {:?}", uni_1);
        let eval = uni_1.evaluate(&0u32.into()) + uni_1.evaluate(&1u32.into());
        assert_eq!(eval, claimed_sum);
    }

    #[test]
    fn test_verifier() {
        // 1st round
        let mut p = Prover::new(poly);
        let mut gi = p.get_univar_poly(0, get_random());
        let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
        assert_eq!(expected_c, claim);
        // middle rounds
        for j in 1..p.f.num_vars {
            let r = get_random();
            expected_c = gi.evaluate(&r);
            gi = p.get_univar_poly(j, r);
            let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
            assert_eq!(expected_c, new_c);
        }
        // final round
        let r = get_random();
        expected_c = gi.evaluate(&r);
        let new_c = p.f.evaluate(&p.steps);
        assert_eq!(expected_c, new_c);
    }
}

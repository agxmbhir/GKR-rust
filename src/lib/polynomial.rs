use crate::lib::field::Field;

struct Polynomial<F: Field> { 
    coefficients: Vec<F>, // degree is coefficients.len() - 1
}

impl<F: Field> Polynomial<F> {
    fn new(coefficients: Vec<F>) -> Self {
        Self { coefficients }
    }
    fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }

    fn add(&self, other: &Self) -> Self {
        let mut coefficients = Vec::new();
        let mut i = 0;
        while i < self.coefficients.len() && i < other.coefficients.len() {
            coefficients.push(self.coefficients[i].add(&other.coefficients[i]));
            i += 1;
        }
        while i < self.coefficients.len() {
            coefficients.push(self.coefficients[i].clone());
            i += 1;
        }
        while i < other.coefficients.len() {
            coefficients.push(other.coefficients[i].clone());
            i += 1;
        }
        Self::new(coefficients)
    }

    fn sub(&self, other: &Self) -> Self {
        let mut coefficients = Vec::new();
        let mut i = 0;
        while i < self.coefficients.len() && i < other.coefficients.len() {
            coefficients.push(self.coefficients[i].sub(&other.coefficients[i]));
            i += 1;
        }
        while i < self.coefficients.len() {
            coefficients.push(self.coefficients[i].clone());
            i += 1;
        }
        while i < other.coefficients.len() {
            coefficients.push(other.coefficients[i].clone());
            i += 1;
        }
        Self::new(coefficients)
    }

    fn mul(&self, other: &Self) -> Self {
        let mut coefficients = Vec::new();
        for _ in 0..self.coefficients.len() + other.coefficients.len() - 1 {
            coefficients.push(F::zero());
        }
        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                coefficients[i + j] = coefficients[i + j].add(&self.coefficients[i].mul(&other.coefficients[j]));
            }
        }
        Self::new(coefficients)
    }
    // Create a multilinear polynomial extension of a field. 
    
    
}
mod field { 
    pub trait Field {
        fn add(&self, other: &Self) -> Self;
        fn sub(&self, other: &Self) -> Self;
        fn mul(&self, other: &Self) -> Self;
        fn div(&self, other: &Self) -> Self;
        fn zero() -> Self;
        fn one() -> Self;
    }
}
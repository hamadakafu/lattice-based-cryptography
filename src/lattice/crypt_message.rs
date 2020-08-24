use crate::cyclotomic::Cyclotomic;
use std::ops::{Add, Mul};

pub struct CryptMessage {
    pub q: i64,
    pub c0: Cyclotomic,
    pub c1: Cyclotomic,
    /// モジュラブースト用
    pub p: i64,
    pub switch_key: (Cyclotomic, Cyclotomic),
}

/// ビット同士の和に相当する
impl Add for CryptMessage {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        CryptMessage {
            c0: (self.c0 + rhs.c0) % self.q,
            c1: (self.c1 + rhs.c1) % self.q,
            ..self
        }
    }
}

impl Mul for CryptMessage {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let d0: Cyclotomic = (self.c0.clone() * rhs.c0.clone()) % self.q;
        let d1: Cyclotomic =
            (self.c1.clone() * rhs.c0.clone() + self.c0.clone() * rhs.c1.clone()) % self.q;
        let d2: Cyclotomic = -(self.c1.clone() * rhs.c1.clone()) % self.q;

        let mut d0_: Cyclotomic =
            (d0 * self.p + self.switch_key.1.clone() * d2.clone()) % (self.p * self.q);
        let mut d1_: Cyclotomic =
            (d1 * self.p + self.switch_key.0.clone() * d2) % (self.p * self.q);
        // dbg!(&d0_, &d1_, self.p);
        d0_.coefficients = d0_
            .coefficients
            .iter()
            .map(|c| {
                let mut rem = c % self.p;
                if rem % 2 != 0 && rem < 0 {
                    rem += self.p;
                } else if rem % 2 != 0 && rem > 0 {
                    rem -= self.p;
                }
                assert!((c - rem) % self.p == 0);
                (c - rem) / self.p
            })
            .collect();
        d1_.coefficients = d1_
            .coefficients
            .iter()
            .map(|c| {
                let mut rem = c % self.p;
                if rem % 2 != 0 && rem < 0 {
                    rem += self.p;
                } else if rem % 2 != 0 && rem > 0 {
                    rem -= self.p;
                }
                assert!((c - rem) % self.p == 0);
                (c - rem) / self.p
            })
            .collect();
        d0_ = d0_ % self.q;
        d1_ = d1_ % self.q;
        // dbg!(&d0_, &d1_, self.p);
        CryptMessage {
            c0: d0_,
            c1: d1_,
            ..self
        }
    }
}

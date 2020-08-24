use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, Sub, SubAssign},
};

use num::Zero;
use rand_distr::{Distribution, Normal};

#[derive(Debug, Clone)]
pub struct Cyclotomic {
    /// TODO: 今mは素数のみしか対応していない
    pub m: usize,
    /// 係数は格子上の座標そのもの
    /// FIXME: $\phi(m)$次元であるべきだが $m-1$次元になっている
    pub coefficients: Vec<i64>,
    /// モジュラス
    pub q: Option<i64>,
}

impl Cyclotomic {
    /// m分体整数
    pub fn with_coef(m: usize, coefficients: Vec<i64>) -> Self {
        assert_eq!(m - 1, coefficients.len(), "coef: {:?}", coefficients);
        Cyclotomic {
            m,
            coefficients,
            q: None,
        }
    }

    /// 一様分布 乱数生成
    pub fn sample_with_module(m: usize, q: i64) -> Self {
        Cyclotomic {
            m,
            coefficients: (0..m - 1)
                .into_iter()
                .map(|_| {
                    let output = rand::random::<i64>() % q;
                    if output as f64 <= -q as f64 / 2. {
                        output + q
                    } else if output as f64 > q as f64 / 2. {
                        output - q
                    } else {
                        output
                    }
                })
                .collect(),
            q: Some(q),
        }
    }
    /// ガウス分布で乱数生成し、結果をround
    /// delta: 標準偏差
    pub fn gaussian_with_delta(m: usize, delta: f64) -> Self {
        let normal = Normal::new(0., delta).unwrap();
        Cyclotomic {
            m,
            coefficients: (0..m - 1)
                .into_iter()
                .map(|_| {
                    let v = normal.sample(&mut rand::thread_rng());
                    v.round() as i64
                })
                .collect(),
            q: None,
        }
    }
}

impl Add for Cyclotomic {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.m, rhs.m);

        Cyclotomic {
            m: self.m,
            coefficients: self
                .coefficients
                .iter()
                .zip(rhs.coefficients.iter())
                .map(|(a, b)| *a + *b)
                .collect(),
            q: self.q,
        }
    }
}

impl AddAssign for Cyclotomic {
    fn add_assign(&mut self, rhs: Self) {
        assert_eq!(self.m, rhs.m);

        self.coefficients = self
            .coefficients
            .iter()
            .zip(rhs.coefficients.iter())
            .map(|(a, b)| *a + *b)
            .collect();
    }
}

impl Sub for Cyclotomic {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.m, rhs.m);

        Cyclotomic {
            m: self.m,
            coefficients: self
                .coefficients
                .iter()
                .zip(rhs.coefficients.iter())
                .map(|(a, b)| *a - *b)
                .collect(),
            q: self.q,
        }
    }
}

impl SubAssign for Cyclotomic {
    fn sub_assign(&mut self, rhs: Self) {
        assert_eq!(self.m, rhs.m);

        self.coefficients = self
            .coefficients
            .iter()
            .zip(rhs.coefficients.iter())
            .map(|(a, b)| *a - *b)
            .collect();
    }
}

/// TODO: 多項式の掛け算は通常$O(m^2)$であるが、fftを用いることにより$O(mlogm)$で実行可能
impl Mul for Cyclotomic {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        assert_eq!(self.m, rhs.m);

        let mut ans_coef = vec![i64::zero(); self.m];
        for (ia, a) in self.coefficients.iter().enumerate() {
            for (ib, b) in rhs.coefficients.iter().enumerate() {
                ans_coef[(ia + ib) % self.m] += *a * *b;
            }
        }
        // FIXME: m-1次は消せる
        for i in 0..self.m - 1 {
            ans_coef[i] -= ans_coef[self.m - 1];
        }
        ans_coef.pop();

        Cyclotomic {
            m: self.m,
            coefficients: ans_coef,
            q: self.q,
        }
    }
}

/// TODO: 多項式の掛け算は通常$O(m^2)$であるが、fftを用いることにより$O(mlogm)$で実行可能
impl MulAssign for Cyclotomic {
    fn mul_assign(&mut self, rhs: Self) {
        assert_eq!(self.m, rhs.m);

        let mut ans_coef = vec![i64::zero(); self.m];
        for (ia, a) in self.coefficients.iter().enumerate() {
            for (ib, b) in rhs.coefficients.iter().enumerate() {
                ans_coef[(ia + ib) % self.m] += *a * *b;
            }
        }
        // m-1次は消せる
        // TODO: 確認、円分整数を表現するのに、m-1次の係数は本当にいらないのか
        for i in 0..self.m - 1 {
            ans_coef[i] -= ans_coef[self.m - 1];
        }
        ans_coef.pop();
        self.coefficients = ans_coef;
    }
}

impl Mul<i64> for Cyclotomic {
    type Output = Self;
    fn mul(self, rhs: i64) -> Self::Output {
        Cyclotomic {
            coefficients: self.coefficients.iter().map(|c| c * rhs).collect(),
            ..self
        }
    }
}

impl MulAssign<i64> for Cyclotomic {
    fn mul_assign(&mut self, rhs: i64) {
        self.coefficients = self.coefficients.iter().map(|c| c * rhs).collect();
    }
}

/// mod p -> (- p / 2, p/2]
/// p = 5 -> {-2, -1, 0, 1, 2}
impl Rem<i64> for Cyclotomic {
    type Output = Self;
    fn rem(self, rhs: i64) -> Self::Output {
        let mut ans_coef = Vec::with_capacity(self.m - 1);
        for c in self.coefficients {
            let output = {
                let out = c % rhs;
                if out as f64 <= -rhs as f64 / 2. {
                    out + rhs
                } else if out as f64 > rhs as f64 / 2. {
                    out - rhs
                } else {
                    out
                }
            };
            ans_coef.push(output);
        }
        Cyclotomic {
            coefficients: ans_coef,
            q: Some(rhs),
            ..self
        }
    }
}

impl Neg for Cyclotomic {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        self.coefficients = self.coefficients.iter().map(|a| -a).collect();
        self
    }
}

impl PartialEq for Cyclotomic {
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(self.m, other.m);
        self.coefficients == other.coefficients
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_add_m3_q65() {
        let m = 3;
        let q = 65;
        let s = Cyclotomic::with_coef(m, vec![1, 1]);
        let a = Cyclotomic::with_coef(m, vec![-19, -8]);
        let e = Cyclotomic::with_coef(m, vec![1, -1]);
        let b = (s.clone() * a.clone() + e.clone()) % q;
        assert_eq!(
            b,
            Cyclotomic {
                m: m,
                coefficients: vec![-10, -20],
                q: None
            }
        );
    }

    #[test]
    fn test_mul_m3_p65() {
        for t in vec![
            (vec![11, -6], vec![21, 15], vec![-4, -1]),
            (vec![-11, -21], vec![21, 15], vec![19, -31]),
            (vec![11, 21], vec![12, -11], vec![-27, -28]),
        ] {
            let c0 = Cyclotomic {
                coefficients: t.0,
                m: 3,
                q: None,
            };

            let c1 = Cyclotomic {
                coefficients: t.1,
                m: 3,
                q: None,
            };

            assert_eq!(Cyclotomic::with_coef(3, t.2), (c0 * c1) % 65);
        }
    }
    #[test]
    fn test_neg_mul_m3_p65() {
        for t in vec![(vec![-11, -21], vec![12, -11], vec![-27, -28]), (vec![-23, -6], vec![-15, -23], vec![-12, -26])] {
            let c0 = Cyclotomic {
                coefficients: t.0,
                m: 3,
                q: None,
            };

            let c1 = Cyclotomic {
                coefficients: t.1,
                m: 3,
                q: None,
            };

            assert_eq!(Cyclotomic::with_coef(3, t.2), -(c0 * c1) % 65);
        }
    }
}

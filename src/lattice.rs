use crate::cyclotomic::Cyclotomic;

mod crypt_message;
use crypt_message::CryptMessage;

struct Lattice {
    m: usize,
    q: i64,
    /// 正規分布の標準偏差
    delta: f64,
    /// モジュラブースト用
    /// TODO: 奇数でないといけない、多分$\delta$で引く際に$\delta$が2とpと合同である必要があるため
    p: i64,
}

impl Lattice {
    pub fn new(m: usize, q: i64, delta: f64, p: i64) -> Self {
        assert!(
            p % 2 == 1,
            "モジュラブースト用のpは奇数でないといけない. p: {}",
            p
        );
        Lattice { m, q, delta, p: p }
    }
    /// $s$
    pub fn gen_sk(&self) -> Cyclotomic {
        Cyclotomic::sample_with_module(self.m, 2)
    }

    /// $(a, b)$
    pub fn gen_pk(&self, sk: &Cyclotomic) -> (Cyclotomic, Cyclotomic) {
        let a = Cyclotomic::sample_with_module(self.m, self.q);
        let e = Cyclotomic::gaussian_with_delta(self.m, self.delta);
        (a.clone(), (a * sk.clone() + e * 2) % self.q)
    }

    /// $(A, B)$、モジュラブースト用
    pub fn gen_switch_key(&self, sk: &Cyclotomic) -> (Cyclotomic, Cyclotomic) {
        let a = Cyclotomic::sample_with_module(self.m, self.p * self.q);
        let e = Cyclotomic::gaussian_with_delta(self.m, self.delta);
        (
            a.clone(),
            (a * sk.clone() - sk.clone() * sk.clone() * self.p + e * 2) % (self.p * self.q),
        )
    }

    /// output: $(c_0, c_1)$
    pub fn encode(
        &self,
        mes: Cyclotomic,
        pk: (Cyclotomic, Cyclotomic),
        switch_key: &(Cyclotomic, Cyclotomic),
    ) -> CryptMessage {
        let v = Cyclotomic::sample_with_module(self.m, 2);
        let e0 = Cyclotomic::gaussian_with_delta(self.m, self.delta);
        let e1 = Cyclotomic::gaussian_with_delta(self.m, self.delta);
        let c0 = (pk.1 * v.clone() + e0 * 2 + mes) % self.q;
        let c1 = (pk.0 * v + e1 * 2) % self.q;
        CryptMessage {
            c0,
            c1,
            q: self.q,
            p: self.p,
            switch_key: switch_key.clone(),
        }
    }

    /// output: $\hat m$
    pub fn decode(&self, sk: Cyclotomic, c: CryptMessage) -> Cyclotomic {
        (c.c0 - sk * c.c1) % self.q % 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_m3_q65_delta2_p67() {
        let m = 3;
        let q = 65;
        let delta = 2.;
        let p = 67;
        let lattice = Lattice::new(m, q, delta, p);
        let mes = Cyclotomic::with_coef(m, vec![1, 0]);
        let sk = Cyclotomic::with_coef(m, vec![1, 1]);
        let pk = (
            Cyclotomic::with_coef(m, vec![-19, -8]),
            Cyclotomic::with_coef(m, vec![-9, -21]),
        );
        let switch_key = lattice.gen_switch_key(&sk);
        let c = lattice.encode(mes.clone(), pk, &switch_key);
        let hat_mes = lattice.decode(sk, c);
        dbg!(&hat_mes);
        assert_eq!(mes, hat_mes);
    }

    #[test]
    fn test_add_m3_q257_delta2_p259() {
        let m = 3;
        let q = 257;
        let delta = 2.;
        let p = 259;
        let lattice = Lattice::new(m, q, delta, p);

        for t in vec![
            (vec![1, 0], vec![1, 0], vec![0, 0]),
            (vec![1, 1], vec![1, 0], vec![0, 1]),
            (vec![1, 1], vec![1, 0], vec![0, 1]),
        ] {
            let mes_a = Cyclotomic::with_coef(m, t.0);
            let mes_b = Cyclotomic::with_coef(m, t.1);

            let sk = Cyclotomic::with_coef(m, vec![1, 1]);
            let pk = (
                Cyclotomic::with_coef(m, vec![-19, -8]),
                Cyclotomic::with_coef(m, vec![-9, -21]),
            );
            let switch_key = lattice.gen_switch_key(&sk);
            let c_a = lattice.encode(mes_a.clone(), pk.clone(), &switch_key);
            let c_b = lattice.encode(mes_b.clone(), pk, &switch_key);
            let hat_mes = lattice.decode(sk, c_a + c_b);
            dbg!(&hat_mes);
            assert_eq!(Cyclotomic::with_coef(m, t.2), hat_mes);
        }
    }

    #[test]
    fn test_mul_m3_q257_delta2_p259() {
        let m = 3;
        let q = 257;
        let delta = 2.;
        let p = 259;
        let lattice = Lattice::new(m, q, delta, p);
        for t in vec![(vec![0, 1], vec![1, 1], vec![1, 0])] {
            let mes_a = Cyclotomic::with_coef(m, t.0);
            let mes_b = Cyclotomic::with_coef(m, t.1);

            let sk = Cyclotomic::with_coef(m, vec![1, 1]);
            let pk = (
                Cyclotomic::with_coef(m, vec![-19, -8]),
                Cyclotomic::with_coef(m, vec![-9, -21]),
            );
            let switch_key = lattice.gen_switch_key(&sk);
            let c_a = lattice.encode(mes_a.clone(), pk.clone(), &switch_key);
            let c_b = lattice.encode(mes_b.clone(), pk, &switch_key);
            let hat_mes = lattice.decode(sk, c_a * c_b);
            dbg!(&hat_mes);
            assert_eq!(Cyclotomic::with_coef(m, t.2), hat_mes);
        }
    }

    extern crate test;
    use test::Bencher;
    #[bench]
    fn bench_m3_q101_delta2_p101(b: &mut Bencher) {
        let m = 3;
        let q = 101;
        let delta = 2.;
        let p = 101;
        let lattice = Lattice::new(m, q, delta, p);
        let mes = Cyclotomic::with_coef(m, vec![1, 0]);
        dbg!(mes.clone());
        let sk = Cyclotomic::with_coef(m, vec![1, 1]);
        let pk = (
            Cyclotomic::with_coef(m, vec![-19, -8]),
            Cyclotomic::with_coef(m, vec![-9, -21]),
        );
        b.iter(|| {
            let switch_key = lattice.gen_switch_key(&sk);
            let c = lattice.encode(mes.clone(), pk.clone(), &switch_key);
            let hat_mes = lattice.decode(sk.clone(), c);
            assert_eq!(mes, hat_mes);
        });
    }

    #[bench]
    fn bench_add_m3_q257_delta2_p503(b: &mut Bencher) {
        let m = 3;
        let q = 257;
        let delta = 2.;
        let p = 503;
        let lattice = Lattice::new(m, q, delta, p);

        b.iter(|| {
            for t in vec![
                (vec![1, 0], vec![1, 0], vec![0, 0]),
                (vec![1, 1], vec![1, 0], vec![0, 1]),
                (vec![1, 1], vec![1, 0], vec![0, 1]),
                (vec![1, 1], vec![1, 1], vec![0, 0]),
            ] {
                let mes_a = Cyclotomic::with_coef(m, t.0);
                let mes_b = Cyclotomic::with_coef(m, t.1);

                let sk = Cyclotomic::with_coef(m, vec![1, 1]);
                let pk = (
                    Cyclotomic::with_coef(m, vec![-19, -8]),
                    Cyclotomic::with_coef(m, vec![-9, -21]),
                );
                let switch_key = lattice.gen_switch_key(&sk);
                let c_a = lattice.encode(mes_a.clone(), pk.clone(), &switch_key);
                let c_b = lattice.encode(mes_b.clone(), pk, &switch_key);
                let hat_mes = lattice.decode(sk, c_a + c_b);
                assert_eq!(Cyclotomic::with_coef(m, t.2), hat_mes);
            }
        });
    }

    /// TODO: パラメータが割とシビアでどういうときに計算結果が合わないかがわからない
    /// deltaを1.5とかにするとわりとこける。1.2とかもたまにこける
    #[bench]
    fn bench_mul_m3_q257_delta1_p259(b: &mut Bencher) {
        let m = 3;
        let q = 507;
        let delta = 1.1;
        let p = 509;
        let lattice = Lattice::new(m, q, delta, p);

        let mut count = 0;
        let sk = Cyclotomic::with_coef(m, vec![1, 1]);
        let pk = (
            Cyclotomic::with_coef(m, vec![-19, -8]),
            Cyclotomic::with_coef(m, vec![-9, -21]),
        );
        let switch_key = lattice.gen_switch_key(&sk);
        b.iter(|| {
            count += 1;
            for t in vec![(vec![0, 1], vec![1, 1], vec![1, 0])] {
                let mes_a = Cyclotomic::with_coef(m, t.0);
                let mes_b = Cyclotomic::with_coef(m, t.1);

                let c_a = lattice.encode(mes_a.clone(), pk.clone(), &switch_key);
                let c_b = lattice.encode(mes_b.clone(), pk.clone(), &switch_key);
                let hat_mes = lattice.decode(sk.clone(), c_a * c_b);
                assert_eq!(Cyclotomic::with_coef(m, t.2), hat_mes, "count: {:?}", count);
            }
        });
    }
}

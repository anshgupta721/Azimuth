#[derive(Debug, Clone)]
pub struct HarmonicCoeffs {
    degree: usize,
    c: Vec<Vec<f64>>,
    s: Vec<Vec<f64>>,
}

impl HarmonicCoeffs {
    pub fn zeros(degree: usize) -> Self {
        let c = (0..=degree).map(|n| vec![0.0; n + 1]).collect();
        let s = (0..=degree).map(|n| vec![0.0; n + 1]).collect();
        Self { degree, c, s }
    }

    pub fn set(&mut self, n: usize, m: usize, cnm: f64, snm: f64) {
        self.c[n][m] = cnm;
        self.s[n][m] = snm;
    }
}

pub struct LegendreCache {
    p: Vec<Vec<f64>>,
    dp: Vec<Vec<f64>>,
}

impl LegendreCache {
    fn new(n_max: usize, phi: f64) -> Self {
        let u = phi.sin();
        let t = phi.cos();
        let dim = n_max + 1;
        let mut p = vec![vec![0.0; dim + 1]; dim + 1];
        let mut dp = vec![vec![0.0; dim + 1]; dim + 1];

        p[0][0] = 1.0;

        // Sectorial (diagonal) recursion: P_nn from P_(n-1)(n-1)
        for n in 1..=n_max {
            let prev = p[n - 1][n - 1];
            let factor = ((2.0 * n as f64 + 1.0) / (2.0 * n as f64)).sqrt();
            p[n][n] = factor * t * prev;
        }

        for n in 1..n_max {
            for m in 0..n {
                let n_f = n as f64;
                let m_f = m as f64;

                let a_nm =
                    (((2.0 * n_f - 1.0) * (2.0 * n_f + 1.0)) / ((n_f - m_f) * (n_f + m_f))).sqrt();

                let p_prev1 = p[n - 1][m];
                let p_prev2 = if n >= 2 && n - 2 >= m {
                    p[n - 2][m]
                } else {
                    0.0
                };

                let b_nm = if n >= 2 && n - 2 >= m {
                    (((2.0 * n_f + 1.0) * (n_f + m_f - 1.0) * (n_f - m_f - 1.0))
                        / ((n_f - m_f) * (n_f + m_f) * (2.0 * n_f - 3.0)))
                        .sqrt()
                } else {
                    0.0
                };
                p[n][m] = a_nm * u * p_prev1 - b_nm * p_prev2;
            }
        }

        for n in 0..=n_max {
            for m in 0..=n{
                let n_f = n as f64;
                let m_f = m as f64;

                let n_nm = if m==n {
                    0.0
                } else if m == 0{
                    (n_f * (n_f + 1.0) / 2.0).sqrt()
                } else {
                    ((n_f - m_f) * (n_f + m_f + 1.0)).sqrt()
                };

                let p_next = if m < n { p[n][m+1] } else { 0.0 };

                dp[n][m] = n_nm * p_next - m_f * t.tan().recip() * 0.0;
                dp[n][m] = n_nm * p_next - m_f * (u / t) * p[n][m];
            }
        }

        LegendreCache { p, dp }
    }
}

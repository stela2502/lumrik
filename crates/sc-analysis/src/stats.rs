pub(crate) fn pearson(a: &[f64], b: &[f64]) -> f64 {
    if a.len() != b.len() || a.is_empty() {
        return 0.0;
    }
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let mut xy = 0.0;
    let mut aa = 0.0;
    let mut bb = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        let dx = x - ma;
        let dy = y - mb;
        xy += dx * dy;
        aa += dx * dx;
        bb += dy * dy;
    }
    let d = (aa * bb).sqrt();
    if d == 0.0 { 0.0 } else { xy / d }
}

pub(crate) fn pearson_f32(a: &[f32], b: &[f32]) -> f64 {
    let aa = a.iter().map(|&x| x as f64).collect::<Vec<_>>();
    let bb = b.iter().map(|&x| x as f64).collect::<Vec<_>>();
    pearson(&aa, &bb)
}

pub(crate) fn mann_whitney(a: &[f64], b: &[f64]) -> (f64, f64) {
    if a.is_empty() || b.is_empty() {
        return (0.0, 1.0);
    }
    let mut v = Vec::with_capacity(a.len() + b.len());
    for &x in a {
        v.push((x, 0u8));
    }
    for &x in b {
        v.push((x, 1u8));
    }
    v.sort_by(|x, y| x.0.total_cmp(&y.0));
    let mut rank_sum = 0.0;
    let mut i = 0;
    while i < v.len() {
        let mut j = i + 1;
        while j < v.len() && v[j].0 == v[i].0 {
            j += 1;
        }
        let rank = (i + j + 1) as f64 / 2.0;
        for q in i..j {
            if v[q].1 == 0 {
                rank_sum += rank;
            }
        }
        i = j;
    }
    let n1 = a.len() as f64;
    let n2 = b.len() as f64;
    let u = rank_sum - n1 * (n1 + 1.0) / 2.0;
    let mu = n1 * n2 / 2.0;
    let sd = (n1 * n2 * (n1 + n2 + 1.0) / 12.0).sqrt();
    let z = if sd > 0.0 { (u - mu).abs() / sd } else { 0.0 };
    let p = erfc(z / std::f64::consts::SQRT_2).clamp(0.0, 1.0);
    (u, p)
}
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398
                                    + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    if x >= 0.0 { r } else { 2.0 - r }
}
pub(crate) fn bh_adjust(p: &[f64]) -> Vec<f64> {
    let m = p.len();
    let mut o = (0..m).collect::<Vec<_>>();
    o.sort_by(|&a, &b| p[a].total_cmp(&p[b]));
    let mut out = vec![1.0; m];
    let mut prev = 1.0f64;
    for (rank, &i) in o.iter().enumerate().rev() {
        let q = (p[i] * m as f64 / (rank + 1) as f64).min(prev).min(1.0);
        out[i] = q;
        prev = q;
    }
    out
}

pub fn quantile(
    values: &[f64],
    q: f64,
) -> f64 {
    assert!(
        !values.is_empty(),
        "cannot calculate quantile of empty data"
    );

    assert!(
        (0.0..=1.0).contains(&q),
        "quantile must be between 0 and 1"
    );

    if values.len() == 1 {
        return values[0];
    }

    let position =
        q * (values.len() - 1) as f64;

    let lower =
        position.floor() as usize;

    let upper =
        position.ceil() as usize;

    if lower == upper {
        return values[lower];
    }

    let fraction =
        position - lower as f64;

    values[lower]
        + fraction
            * (values[upper] - values[lower])
}



pub fn percent(
    value: usize,
    total: usize,
) -> f64 {
    if total == 0 {
        0.0
    } else {
        value as f64
            / total as f64
            * 100.0
    }
}
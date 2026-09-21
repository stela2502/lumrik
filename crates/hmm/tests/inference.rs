use hmm::{CategoricalEmission, GaussianEmission, Hmm, StateId};

#[test]
fn categorical_viterbi_and_posterior_find_two_regimes() {
    let emissions = vec![
        CategoricalEmission::new(vec![0.95, 0.05]).unwrap(),
        CategoricalEmission::new(vec![0.05, 0.95]).unwrap(),
    ];
    let hmm = Hmm::new(vec![0.9, 0.1], vec![0.97, 0.03, 0.03, 0.97], emissions).unwrap();
    let observations = vec![0usize, 0, 0, 0, 1, 1, 1, 1];
    let result = hmm.infer(&observations).unwrap();
    assert_eq!(&result.viterbi()[0..4], &[StateId(0); 4]);
    assert_eq!(&result.viterbi()[4..8], &[StateId(1); 4]);
    assert!(result.posterior(1, StateId(0)) > 0.9);
    assert!(result.posterior(6, StateId(1)) > 0.9);
    for i in 0..observations.len() {
        assert!((result.posterior_row(i).iter().sum::<f64>() - 1.0).abs() < 1e-9);
    }
}

#[test]
fn gaussian_same_engine_finds_low_and_high_signal() {
    let emissions = vec![
        GaussianEmission::new(0.0, 0.25).unwrap(),
        GaussianEmission::new(5.0, 0.25).unwrap(),
    ];
    let hmm = Hmm::new(vec![0.99, 0.01], vec![0.98, 0.02, 0.02, 0.98], emissions).unwrap();
    let observations = vec![0.0, 0.1, -0.2, 0.2, 4.8, 5.1, 5.0, 5.2];
    let result = hmm.infer(&observations).unwrap();
    assert!(result.viterbi()[..4].iter().all(|x| *x == StateId(0)));
    assert!(result.viterbi()[4..].iter().all(|x| *x == StateId(1)));
}

#[test]
fn baum_welch_moves_gaussian_emissions_toward_data() {
    let emissions = vec![
        GaussianEmission::new(-1.0, 4.0).unwrap(),
        GaussianEmission::new(6.0, 4.0).unwrap(),
    ];
    let mut hmm = Hmm::new(vec![0.9, 0.1], vec![0.9, 0.1, 0.1, 0.9], emissions).unwrap();
    let observations = vec![0.0, 0.2, -0.1, 0.1, 5.0, 5.2, 4.9, 5.1];
    hmm.baum_welch(&observations, 8, 1e-9).unwrap();
    assert!(hmm.emissions()[0].mean().abs() < 0.5);
    assert!((hmm.emissions()[1].mean() - 5.05).abs() < 0.5);
}

# hmm AI Contract

## Role

This is a generic log-space hidden Markov model engine. It deliberately contains no genome/protein/assay-specific biology.

## Generic boundary

Biological meaning enters through observation, emission, and result-processing traits. Do not bake a specific biological modality into the core HMM implementation.

## Numerics

Probability calculations remain in log space for numerical stability. Avoid conversions to ordinary probability space inside core recurrences unless mathematically required and tested.

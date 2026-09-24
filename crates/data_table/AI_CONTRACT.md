# data_table AI Contract

## Role

A lightweight tabular data structure for numeric columns and categorical factors.

## Generic boundary

Keep the core table representation domain-agnostic. Biological interpretation belongs in caller crates.

## Factors

Categorical/factor level identity and row alignment are data semantics. Filtering/reordering rows must keep every column/factor synchronized.

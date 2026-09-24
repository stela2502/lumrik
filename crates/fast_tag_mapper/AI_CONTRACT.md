# fast_tag_mapper AI Contract

## Role

Fast supplemental-feature mapping uses exact short seeds to nominate candidates and a fuller verification step to decide matches.

## Seeds are candidates, not calls

A seed hit MUST NOT by itself become a feature assignment. Full candidate verification remains authoritative.

## Orientation and quality

Forward/reverse handling and quality-weighted verification are part of matching semantics. Performance optimizations must preserve them.

## Thresholds

`min_hits` and similar thresholds affect sensitivity/specificity and must remain explicit/configurable. Do not hard-code dataset-specific values into the core mapper.

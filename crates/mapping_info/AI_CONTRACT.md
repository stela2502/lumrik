# mapping_info AI Contract

## Role

`mapping_info` is shared run/accounting state used to carry mapping/analysis measurements across Lumrik stages.

## Metric semantics

A metric name must retain one meaning across producers and consumers. Do not reuse a counter for a superficially similar population.

Counters used by reports/health servers are part of the observability contract. If the underlying population changes, update the metric name/consumer rather than silently changing its meaning.

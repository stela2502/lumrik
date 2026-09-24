# lumrik-status AI Contract

## Role

`lumrik-status` provides reusable live status/dashboard infrastructure. It renders state supplied by tools; it should not invent biological interpretations.

## Truthful UI

Do not display unavailable measurements as if measured zeros. Prefer hiding/not-applicable state when a mode cannot produce a metric.

Status rendering must remain lightweight enough not to perturb the analysis it observes.

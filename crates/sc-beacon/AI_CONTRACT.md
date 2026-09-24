# sc-beacon AI Contract

## Role

`sc-beacon` distinguishes genuine cellular feature signal from experimentally observed ambient/background signal. It is a statistical caller for feature observations, not a GEX cell-definition engine unless explicitly invoked for that purpose by a caller.

## Background matters

Ambient/background observations are part of the model, not optional decoration. Do not silently substitute filtered cells for measured background or invent background when none exists.

## Empty populations

Empty filtered/background populations are meaningful failure/edge states. Handle them explicitly; do not manufacture calls merely to keep the pipeline moving.

## Outputs

Posterior/model statistics, guide calls, and cell assignments must remain traceable to the observations/model used. Threshold changes must not silently alter the meaning of existing output columns.

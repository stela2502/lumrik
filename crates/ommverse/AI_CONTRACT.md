# ommverse AI Contract

## Role

Ommverse is Lumrik's integrated genome-to-protein biological reference model. It separates stable biological identity from source/import formats.

## Source versus model

GTF, twoBit, UniProt/bigBed and similar resources are import sources. Downstream callers should consume normalized genes/transcripts/proteins/features rather than depend on source-file quirks.

## Provenance

Reference source/version and mapping provenance must remain traceable. Missing or conflicting source information must not be silently converted into certainty.

## Identity

Do not conflate gene, transcript, protein, isoform, and feature identifiers merely because a source offers convenient cross-references.

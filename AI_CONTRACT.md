# Lumrik AI Contract

MANDATORY FOR AI/CODING AGENTS:

1. Read this file BEFORE modifying this repository.
2. Before modifying any crate, read `crates/<crate>/AI_CONTRACT.md` if it exists.
3. If a change crosses crate boundaries, read ALL affected crate contracts.
4. These contracts describe intentional architecture and biological semantics.
   Do not override them based on assumptions inferred from the implementation.
5. When debugging reveals a new non-obvious invariant, update the appropriate
   crate contract as part of the fix.

## Maintaining contracts

When debugging reveals a non-obvious invariant, failure mode, or design
decision that would be useful to a future coding agent, update the relevant
`AI_CONTRACT.md` as part of the fix.

Do not turn contracts into changelogs.

Record the lesson/invariant, not the debugging history.

Bad:
"On September 23 we discovered that cell calling returned 275 cells."

Good:
"Canonical GEX cell calling MUST use only GEX-origin (`G`) molecules.
VDJ-origin (`V`) molecules will be quantified by the specific sc-vdj crate."

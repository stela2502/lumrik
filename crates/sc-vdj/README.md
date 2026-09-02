sc-vdj
sc-vdj is Lumrik’s posterior single-cell antigen-receptor reconstruction tool.

It is designed to work after normal gene-expression mapping. It does not require a dedicated V(D)J library. Instead, it revisits the retained Nelrune BAM together with the called-cell expression matrix and reconstructs immunoglobulin and T-cell receptor rearrangements from receptor-locus evidence that is already present in the mapped RNA reads.

The seven receptor chains handled explicitly are:

IGH
IGK
IGL
TRA
TRB
TRG
TRD
The production entry point is nelrune-vdj.

What sc-vdj needs

A normal analysis requires:

the Nelrune exonic MEX directory;
the retained Nelrune mapper BAM containing cell and UMI identity; and
a V(D)J reference index (*.vdjidx).
The V(D)J index is generated from the same genome FASTA and GTF annotation used for mapping.

A reference can also be rebuilt from GTF + genome at run time, but the persisted .vdjidx is the preferred production path.

Building a V(D)J reference index

Build the index once:

vdj-index \
    --gtf annotation.gtf \
    --genome genome.fa \
    --out mouse.vdjidx
The index stores the receptor germline reference and the segment ordering used by the compact HC/LC identifiers described below.

Because the packed identifiers refer to segment indices in this reference, keep the .vdjidx used for an analysis together with the results. An HC/LC identifier is intentionally decoded through the reference that created it.

Running sc-vdj

Typical production run:

nelrune-vdj \
    --exonic /path/to/nelrune/exonic \
    --bam /path/to/nelrune.mapper.bam \
    --index /path/to/mouse.vdjidx \
    --out /path/to/vdj
To also write observed and inferred-naive rearrangement FASTA files:

nelrune-vdj \
    --exonic /path/to/nelrune/exonic \
    --bam /path/to/nelrune.mapper.bam \
    --index /path/to/mouse.vdjidx \
    --out /path/to/vdj \
    --write-sequences
For legacy BD matrices whose barcodes were padded after conversion, the original BAM cell-barcode length can be specified, for example:

--cell-barcode-len 27
nelrune-vdj starts the Lumrik live status server by default on port 8787. In batch environments where listening sockets are undesirable:

--no-health-server
How the BAM is processed

The current implementation is deliberately single-pass and does not create temporary BAM shards.

For each alignment, the fast path is approximately:

BAM record
  |
  +-- secondary alignment? ----------------------> reject
  |
  +-- chromosome has indexed receptor segments? -> reject if no
  |
  +-- aligned block overlaps an actual V/D/J/C
  |   segment interval? --------------------------> reject if no
  |
  +-- resolve cell + UMI
  |
  +-- cell is present in the called-cell matrix? -> reject if no
  |
  `-- route evidence to the matching receptor chain(s)
Evidence is kept separately for IGH, IGK, IGL, TRA, TRB, TRG and TRD. It is not first mixed into one general per-cell receptor pool.

Chromosome identity is only the first inexpensive rejection step. Final routing uses the actual indexed segment intervals. This matters because receptor loci can share chromosomes and, in particular, TRA and TRD occupy overlapping/nested genomic territory.

After the BAM pass, cells are reconstructed in parallel.

What constitutes a receptor call?

A receptor call is built from UMI-aware evidence and local alignment to the indexed germline segments.

For VJ chains:

IGK, IGL, TRA, TRG
V ---------------- J
For D-bearing chains:

IGH, TRB, TRD
V -------- D -------- J
The D segment is treated as a bounded hypothesis inside the established V-to-J junction. All compatible D candidates in that bounded region can be compared, and the best-supported D is selected deterministically. D evidence refines a heavy-chain identity; it is not allowed to invent a V/J rearrangement by itself.

A cell may therefore have convincing segment evidence but still lack a compact HC/LC identifier if the junction geometry could not be resolved sufficiently. A blank HC/LC value does not necessarily mean “no receptor evidence”.

A useful interpretation is:

V/J segment evidence present
        !=
junction-resolved recombination identity
The compact HC/LC identifiers are only emitted for junction-resolved rearrangements.

Structural clonotyping: clones are not defined by observed sequence
The central design choice in sc-vdj is that a recombination identity is structural and nucleotide-independent.

The observed RNA sequence is retained as evidence and can contain:

sequencing errors;
unresolved bases;
somatic mutations;
transcript-specific sequence variation.
Those bases are valuable for inspection, but they are deliberately not used as the canonical clone identifier.

Instead, the stable identity describes the underlying recombination event:

which V segment was selected;
which D segment was selected, when applicable;
which J segment was selected;
how much germline sequence was deleted at the junction boundaries;
the lengths of inferred P and N additions;
the amount of D retained for D-bearing chains; and
whether an alternative P/N interpretation was selected.
This is the information represented by the HC: and LC: identifiers.

Therefore two observations do not need byte-for-byte identical receptor sequence to represent the same structural recombination.

Conversely, two cells that use the same V and J genes are not automatically the same clone. Their recombination geometry can differ.

For clone-oriented downstream analysis, the packed recombination IDs are the canonical structural keys:

same HC ID  -> same resolved heavy-role recombination structure
same LC ID  -> same resolved light-role recombination structure
A strict paired receptor clone can be defined downstream by the pair:

(HC ID, LC ID)
when both are available.

This separation is intentional: sc-vdj reconstructs and reports receptor identities; downstream clone/tree analysis can choose whether to require both receptor roles, use one role, incorporate biological lineage information, or tolerate unresolved secondary chains.

The current cell-level output reports the strongest heavy-role and strongest light-role rearrangement rather than silently collapsing all receptor evidence into a sequence-derived clone hash.

What are HC: and LC: HEX values?
Examples:

HC:10A70F6105200156006111
LC:11A91DB011031
These strings are not hashes.

They are compact, reversible encodings.

HC and LC describe recombination geometry:

HC = V-D-J geometry: IGH, TRB or TRD
LC = V-J geometry: IGK, IGL, TRA or TRG
The prefix therefore means “heavy-role / D-bearing geometry” or “light-role / VJ geometry”. The actual receptor chain is recovered from the referenced germline segments when the identifier is decoded.

Constant-region identity is not part of the HC/LC ID. Isotype/constant information can change independently of the original V(D)J recombination and is reported separately.

Packed segment identities

Each V/D/J segment is represented by a 12-bit index into the persisted V(D)J reference:

000 .. FFF
The decoder verifies that:

referenced indices exist;
V really resolves to a V segment;
D really resolves to a matching D segment for HC;
J really resolves to a J segment; and
all segments belong to a compatible receptor chain.
This is why the matching .vdjidx is required to decode an ID safely.

Packed measurements

Junction measurements are stored compactly in hexadecimal.

Values 0..14 are stored as one hex digit. Larger values are represented as:

Fxxxx
where xxxx is a four-digit hexadecimal u16.

The first nibble is the encoding version. The final 0/1 field records the P/N alternative flag.

Conceptually an LC contains:

version
V segment
J segment
V 3' deletion
P-V length
N length
P-J length
J 5' deletion
P/N alternative flag
An HC additionally contains:

D segment
P-D5 length
D 5' deletion
D retained length
D 3' deletion
P-D3 length
second N length
The exact human-readable values should normally be obtained with vdj-decode, rather than by manually counting hex digits.

Decoding HC/LC identifiers
Use the same V(D)J index that was used for the analysis:

vdj-decode \
    --index mouse.vdjidx \
    HC:10A70F6105200156006111
or:

vdj-decode \
    --index mouse.vdjidx \
    LC:11A91DB011031
vdj-decode prints fields such as:

type
id
chain
v
v_del_3
p_v3_len
n1_len
d
d_del_5
d_retained_len
d_del_3
n2_len
p_j5_len
j_del_5
j
pn_alternative
complete
Multiple IDs can also be sent through standard input:

printf '%s\n' \
    'HC:10A70F6105200156006111' \
    'LC:11A91DB011031' \
    | vdj-decode --index mouse.vdjidx
The decoder is the authoritative way to inspect the packed identifiers.

Main output tables
vdj_receptors.tsv

This is the most convenient cell-level table for biological inspection.

For each cell it contains the strongest reconstructed heavy-role and light-role receptor, including:

heavy HC ID;
light LC ID;
receptor chain;
V/D/J/C assignments;
D-hypothesis diagnostics for heavy chains;
observed V/D/J pieces;
inferred-naive recombination;
RAG1/RAG2/DNTT expression and recombination activity; and
light-chain status.
This is usually the first table to show to a domain expert.

vdj_calls.tsv

This is the detailed rearrangement table.

It contains one reconstructed biological rearrangement per row, including:

the packed recombination ID when the junction is resolved;
chain and recombination stage;
V/D/J/C support;
read and UMI support;
local-alignment scores;
locus position;
D-hypothesis diagnostics;
all measured deletion/P/N fields;
germline segment sequence;
observed segment sequence;
inferred-naive segment sequence;
complete observed rearrangement; and
complete inferred-naive recombination.
This is the table to inspect when asking why a particular HC/LC ID was assigned.

vdj-mapping-info.txt

Contains performance counters and timing information for the run.

The live dashboard additionally exposes the BAM routing funnel, including how many alignments were:

scanned;
on receptor-bearing chromosomes;
overlapping actual indexed receptor segments;
assigned to called cells; and
routed to each of the seven receptor chains.
Optional FASTA

With --write-sequences:

vdj_observed.fasta
vdj_naive.fasta
These are audit/inspection artifacts. They are not the canonical clone identifiers.

Inspecting a suspicious or interesting call
For a biological inspection workflow:

Find the cell in vdj_receptors.tsv.
Note its HC: and/or LC: identifier.
Decode the identifier using vdj-decode --index.
Find the corresponding row(s) in vdj_calls.tsv.
Inspect UMI support, V/D/J assignments, junction geometry and observed/inferred-naive sequence.
If needed, rerun with --write-sequences for FASTA-level inspection.
Do not interpret an empty HC/LC field as a hard negative without checking the call table. It may mean that segment evidence exists but the complete junction geometry remained unresolved.

Useful auxiliary tools
vdj-index

Builds the persisted *.vdjidx reference.

nelrune-vdj

Production single-cell V(D)J reconstruction from Nelrune exonic MEX + retained BAM.

vdj-decode

Decodes reversible HC: / LC: structural recombination identifiers. Requires the matching .vdjidx.

vdj-summary

Developer/detailed report path using the same BAM + MEX + reference inputs.

vdj-rich-cell

Creates a reproducible one-cell integration fixture from a real Nelrune dataset by selecting a receptor-rich cell (or a requested barcode), extracting its BAM records, and writing a one-cell exonic MEX. This is mainly intended for regression/debugging.

Important interpretation boundaries
HC/LC ID is not a sequence hash

Never interpret the hex payload as a hash of CDR3, observed receptor sequence, or amino-acid sequence. It is a reversible encoding of reference segment identity plus inferred recombination geometry.

Constant region is separate

C/isotype is reported, but is not part of the structural recombination identifier.

Mutation is separate from clone identity

Observed mutations remain visible in sequence outputs. They do not redefine the canonical structural recombination event.

Missing ID can still have evidence

A V/J or V/D/J segment call can exist without an HC/LC ID when no sufficiently resolved continuous junction measurement was obtained.

Reference index matters

Packed segment numbers are meaningful only with the persisted V(D)J reference used to create them. Archive the .vdjidx with results that contain HC/LC IDs.

Minimal reproducible archive
For someone who needs to interpret a completed run, the most useful bundle is:

README.md
vdj_receptors.tsv
vdj_calls.tsv
vdj-mapping-info.txt
<reference>.vdjidx
vdj-decode
Optionally add:

vdj_observed.fasta
vdj_naive.fasta
The decoder binary must match the platform on which it will be used. If sharing only tables for review, the binary is unnecessary; include it when the recipient should be able to decode HC/LC IDs themselves.

License

sc-vdj is part of Lumrik and is licensed under AGPL-3.0-or-later unless otherwise specified by the Lumrik project.

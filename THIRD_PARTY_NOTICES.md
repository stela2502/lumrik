# Third-Party Notices

Lumrik is built using open-source software and can interoperate with external bioinformatics tools developed by other projects.

This document distinguishes between software dependencies incorporated into Lumrik builds and external programs that Lumrik can invoke at runtime.

## Rust dependencies

Lumrik uses third-party Rust crates distributed under their respective licenses.

The authoritative list of dependencies for a particular Lumrik version is defined by the workspace `Cargo.toml` files and the corresponding `Cargo.lock`.

Each dependency remains subject to the license terms specified by its respective authors and distributors.

The Lumrik license does not replace, modify, or supersede those licenses.

For redistribution of compiled Lumrik binaries, distributors should review the licenses of the exact dependency versions contained in the corresponding `Cargo.lock`.

Tools such as `cargo-license` may be useful for generating a dependency-specific license inventory for a release.

## External bioinformatics software

Lumrik can interact with external sequence alignment and bioinformatics programs, including:

* STAR
* minimap2
* BWA

These programs are independent projects and are not part of Lumrik merely because Lumrik can launch or communicate with them.

Users are responsible for installing external tools required by their selected workflow and for complying with the licenses applicable to those tools.

Unless explicitly stated for a particular distribution, Lumrik does not redistribute these external programs.

## Reference data and genomic resources

Lumrik workflows may use externally obtained resources such as:

* reference genome sequences;
* genome annotations;
* mapper indices;
* VCF files;
* barcode whitelists;
* primer sequences;
* feature reference sequences.

Such data may be subject to separate licenses, terms of use, attribution requirements, or redistribution restrictions.

Users are responsible for ensuring that they have the necessary rights to use and redistribute reference data employed in their analyses.

## Sequencing platforms and trademarks

Names of sequencing technologies, platforms, assays, software, and organisations may appear in Lumrik documentation to describe compatibility or supported data formats.

Such names and trademarks remain the property of their respective owners.

Mention of a product, platform, or organisation does not imply endorsement, sponsorship, or affiliation unless explicitly stated.

## Reporting licensing issues

If you believe a required attribution, copyright notice, or third-party license has been omitted or represented incorrectly, please report the issue so that it can be corrected.

This file is intended as an overview. The license files and notices distributed with individual third-party components remain authoritative.


## Lumrik + STAR runtime container

The optional Norn runtime container published from this repository includes
STAR in addition to the Lumrik executables. STAR is redistributed through the
pinned BioContainers base image under STAR's MIT license and remains subject to
its upstream copyright and notice terms. The container image does not relicense STAR under the Lumrik
license.

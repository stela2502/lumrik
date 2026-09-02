# Runtime image for Norn/Lumrik.
# STAR comes from the pinned BioContainers image; Lumrik binaries are static
# x86_64-musl release binaries built before this image is assembled.
FROM quay.io/biocontainers/star:2.7.11b--h5ca1c30_8

LABEL org.opencontainers.image.title="Lumrik + STAR runtime"
LABEL org.opencontainers.image.description="Runtime image for Lumrik/Norn single-cell processing"
LABEL org.opencontainers.image.source="https://github.com/stela2502/lumrik"
LABEL org.opencontainers.image.licenses="AGPL-3.0-or-later AND MIT"

COPY target/x86_64-unknown-linux-musl/release/nelrune /usr/local/bin/nelrune
COPY target/x86_64-unknown-linux-musl/release/nelrune-vdj /usr/local/bin/nelrune-vdj
COPY target/x86_64-unknown-linux-musl/release/vdj-index /usr/local/bin/vdj-index
COPY target/x86_64-unknown-linux-musl/release/gtf-splice-index /usr/local/bin/gtf-splice-index
COPY target/x86_64-unknown-linux-musl/release/lumrik-guides /usr/local/bin/lumrik-guides

RUN chmod 0755 \
      /usr/local/bin/nelrune \
      /usr/local/bin/nelrune-vdj \
      /usr/local/bin/vdj-index \
      /usr/local/bin/gtf-splice-index \
      /usr/local/bin/lumrik-guides \
    && STAR --version \
    && nelrune --help >/dev/null \
    && nelrune-vdj --help >/dev/null \
    && vdj-index --help >/dev/null \
    && gtf-splice-index --help >/dev/null \
    && lumrik-guides --help >/dev/null

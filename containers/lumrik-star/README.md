# Lumrik + STAR runtime image

This image is the runtime environment intended for Norn. It contains a pinned
STAR installation plus the static x86_64-musl Lumrik command-line programs used
by the workflow. It deliberately does not contain Rust, Cargo, the Lumrik source
tree, or build dependencies.

The Docker build context is the Lumrik repository root because the Dockerfile
copies release binaries from:

    target/x86_64-unknown-linux-musl/release/

Build the required binaries first:

    rustup target add x86_64-unknown-linux-musl
    cargo build --locked --release --target x86_64-unknown-linux-musl \
      --bin nelrune \
      --bin nelrune-vdj \
      --bin vdj-index \
      --bin gtf-splice-index \
      --bin lumrik-guides

Then build locally from the repository root:

    docker build \
      -f containers/lumrik-star/Dockerfile \
      -t ghcr.io/stela2502/lumrik:dev .

The image redistributes STAR as part of the selected BioContainers base image.
STAR remains governed by its own license. See the project third-party notices
and the upstream image/package metadata.

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

echo "Lumrik root: $ROOT"

# Ensure root Cargo.toml has the workspace dependency.
if ! grep -qE '^\[workspace\.dependencies\]' Cargo.toml; then
    cat >> Cargo.toml <<'EOF'

[workspace.dependencies]
rust-htslib = { path = "crates/rust-htslib" }
EOF
elif ! grep -qE '^rust-htslib\s*=' Cargo.toml; then
    perl -0pi -e '
        s/(\[workspace\.dependencies\]\s*)/$1rust-htslib = { path = "crates\/rust-htslib" }\n/
    ' Cargo.toml
else
    perl -0pi -e '
        s/^rust-htslib\s*=.*$/rust-htslib = { path = "crates\/rust-htslib" }/m
    ' Cargo.toml
fi

echo
echo "Normalizing crate dependencies..."

find crates \
    -mindepth 2 \
    -maxdepth 2 \
    -name Cargo.toml \
    -print0 |
while IFS= read -r -d '' manifest; do
    # Do not rewrite rust-htslib's own Cargo.toml.
    if [[ "$manifest" == "crates/rust-htslib/Cargo.toml" ]]; then
        continue
    fi

    if grep -qE '^rust-htslib\s*=' "$manifest"; then
        echo "  fixing $manifest"

        perl -0pi -e '
            s/^rust-htslib\s*=.*$/rust-htslib.workspace = true/m
        ' "$manifest"
    fi
done

echo
echo "Remaining rust-htslib declarations:"
grep -Rn \
    --include='Cargo.toml' \
    'rust-htslib' \
    Cargo.toml crates || true

echo
echo "Updating Cargo lockfile..."
cargo update

echo
echo "rust-htslib instances in dependency graph:"
cargo tree --workspace | grep 'rust-htslib' | sort -u || true

echo
echo "Done."
echo
echo "Review changes with:"
echo "  git diff -- Cargo.toml 'crates/*/Cargo.toml' Cargo.lock"

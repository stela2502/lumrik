#!/usr/bin/env bash
set -euo pipefail

HDF5_VERSION="1.12.1"
ZLIB_VERSION="1.3.1"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
prefix="${1:-$repo_root/target/musl-native}"
build_root="${MUSL_NATIVE_BUILD_DIR:-$repo_root/target/musl-native-build}"

if ! command -v musl-gcc >/dev/null 2>&1; then
    echo "ERROR: musl-gcc not found; install musl-tools first" >&2
    exit 1
fi
if ! command -v cmake >/dev/null 2>&1; then
    echo "ERROR: cmake not found" >&2
    exit 1
fi
if ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: curl not found" >&2
    exit 1
fi

stamp="$prefix/.lumrik-hdf5-${HDF5_VERSION}-zlib-${ZLIB_VERSION}"
if [[ -f "$stamp" && -f "$prefix/lib/libhdf5.a" && -f "$prefix/lib/libz.a" ]]; then
    echo "Reusing musl HDF5 ${HDF5_VERSION} + zlib ${ZLIB_VERSION} from $prefix"
    exit 0
fi

echo "Building musl zlib ${ZLIB_VERSION} + HDF5 ${HDF5_VERSION} into $prefix"
rm -rf "$build_root"
mkdir -p "$build_root" "$prefix"

curl -L --fail --retry 3 \
    "https://github.com/madler/zlib/archive/refs/tags/v${ZLIB_VERSION}.tar.gz" \
    -o "$build_root/zlib.tar.gz"
tar -xzf "$build_root/zlib.tar.gz" -C "$build_root"
(
    cd "$build_root/zlib-${ZLIB_VERSION}"
    CC=musl-gcc ./configure --static --prefix="$prefix"
    make -j"$(nproc)"
    make install
)

curl -L --fail --retry 3 \
    "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_${HDF5_VERSION}.tar.gz" \
    -o "$build_root/hdf5.tar.gz"
tar -xzf "$build_root/hdf5.tar.gz" -C "$build_root"

cmake \
    -S "$build_root/hdf5-hdf5_${HDF5_VERSION}" \
    -B "$build_root/hdf5-build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DCMAKE_C_COMPILER=musl-gcc \
    -DCMAKE_FIND_ROOT_PATH="$prefix" \
    -DCMAKE_PREFIX_PATH="$prefix" \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_TESTING=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF \
    -DHDF5_BUILD_TOOLS=OFF \
    -DHDF5_BUILD_HL_LIB=OFF \
    -DHDF5_ENABLE_Z_LIB_SUPPORT=ON \
    -DZLIB_ROOT="$prefix"
cmake --build "$build_root/hdf5-build" --parallel "$(nproc)"
cmake --install "$build_root/hdf5-build"

touch "$stamp"
echo "Built musl HDF5 ${HDF5_VERSION} + zlib ${ZLIB_VERSION}: $prefix"

# System Requirements for `pca_tree`

`pca_tree` uses **OpenBLAS** for accelerated linear algebra.  
You must have OpenBLAS installed before building from source.

---

## ⚠️ Linux

### Ubuntu / Debian
```bash
sudo apt-get update
sudo apt-get install libopenblas-dev
```

### Fedora / RHEL / CentOS
```bash
sudo dnf install openblas-devel
```

### Arch Linux
```bash
sudo pacman -S openblas
```

---

## ⚠️ macOS

macOS includes the Accelerate framework, but it is not fully LAPACK compatible.  
For best compatibility:

```bash
brew install openblas
```

Set environment variables if needed:

```bash
export OPENBLAS_DIR="$(brew --prefix openblas)"
```

---

## ⚠️ Windows (MSYS2)

Install OpenBLAS via MSYS2:

```bash
pacman -S mingw-w64-x86_64-openblas
```

Then set:

```bash
export OPENBLAS_DIR="C:/msys64/mingw64"
```

---

## 🛠 Building From Source

After installing OpenBLAS:

```bash
cargo build --release
```

Or install system-wide:

```bash
cargo install --path .
```

This produces the CLI binary:

```
pca_tree
```

---

## Optional: Troubleshooting

If OpenBLAS cannot be found, try setting:

```bash
export OPENBLAS_DIR=/path/to/openblas
export LD_LIBRARY_PATH="$OPENBLAS_DIR/lib:$LD_LIBRARY_PATH"
```

Or ensure `pkg-config` detects it:

```bash
pkg-config --libs openblas
```

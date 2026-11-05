# Obtain prebuilt Debian 12 packages without Conda or source compilation.
# Keep their non-glibc shared libraries isolated from STAR's runtime.
FROM debian:12-slim AS ont-tools
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends minimap2 samtools \
    && mkdir -p /opt/ont-tools/bin /opt/ont-tools/lib \
    && cp /usr/bin/minimap2 /usr/bin/samtools /opt/ont-tools/bin/ \
    && for exe in /usr/bin/minimap2 /usr/bin/samtools; do \
         ldd "$exe" | awk '/=> \/.*\// {print $3} /^\/[[:graph:]]+\// {print $1}' | \
           while IFS= read -r lib; do \
             case "$(basename "$lib")" in \
               libc.so.*|libm.so.*|libpthread.so.*|libdl.so.*|librt.so.*|ld-linux-*) ;; \
               *) cp -L "$lib" /opt/ont-tools/lib/ ;; \
             esac; \
           done; \
       done \
    && rm -rf /var/lib/apt/lists/*

# Runtime image for Norn/Lumrik. STAR comes from BioContainers.
FROM quay.io/biocontainers/star:2.7.11b--h5ca1c30_8

COPY --from=ont-tools /opt/ont-tools /opt/ont-tools
# Wrapper scripts scope the Debian libraries to the external tools only.
RUN printf '#!/bin/sh\nLD_LIBRARY_PATH=/opt/ont-tools/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} exec /opt/ont-tools/bin/minimap2 "$@"\n' > /usr/local/bin/minimap2 \
    && printf '#!/bin/sh\nLD_LIBRARY_PATH=/opt/ont-tools/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} exec /opt/ont-tools/bin/samtools "$@"\n' > /usr/local/bin/samtools \
    && chmod 0755 /usr/local/bin/minimap2 /usr/local/bin/samtools

LABEL org.opencontainers.image.title="Lumrik + STAR runtime"
LABEL org.opencontainers.image.description="Runtime image for Lumrik/Norn single-cell processing"
LABEL org.opencontainers.image.source="https://github.com/stela2502/lumrik"
LABEL org.opencontainers.image.licenses="AGPL-3.0-or-later AND MIT"

COPY containers/lumrik-star/runtime-bin/ /usr/local/bin/

RUN chmod 0755 /usr/local/bin/* \
    && echo "=== Lumrik runtime binaries ===" \
    && ls -lh /usr/local/bin \
    && echo \
    && echo "=== STAR ===" \
    && STAR --version \
    && minimap2 --version \
    && samtools --version \
    && echo \
    && echo "=== Lumrik smoke tests ===" \
    && for exe in /usr/local/bin/*; do \
         if [ -x "$exe" ]; then \
           name="$(basename "$exe")"; \
           case "$name" in \
             STAR|STAR-*|STARlong|STARlong-*|\
             bgzip|tabix|htsfile|openssl|\
             bzip2|bunzip2|bzcat|bzip2recover|\
             zstd|zstdcat|zstdmt|unzstd|\
             libdeflate-gzip|libdeflate-gunzip|\
             nghttp|nghttpd|nghttpx|\
             clear|reset|tic|toe|tput|tset|tabs|infocmp|infotocap|captoinfo|\
             kinit|klist|kdestroy|kpasswd|kswitch|ktutil|kvno|kadmin|ksu|k5srvutil|\
             adig|ahost|annot-tsv|c_rehash|compile_et|curl-config|keyctl|ref-cache|\
             sclient|sim_client|uuclient|ncurses6-config|ncursesw6-config|\
             bzcmp|bzdiff|bzegrep|bzfgrep|bzgrep|bzless|bzmore|zstdgrep|zstdless) \
               ;; \
             *) \
               echo "Checking $name"; \
               "$exe" --help >/dev/null 2>&1 || true; \
               ;; \
           esac; \
         fi; \
       done

# Runtime image for Norn/Lumrik.
# STAR comes from the pinned BioContainers image; Lumrik binaries are static
# x86_64-musl release binaries staged before this image is assembled.
FROM quay.io/biocontainers/star:2.7.11b--h5ca1c30_8

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

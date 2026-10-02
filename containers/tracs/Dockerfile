# Portable TRACS image, built from source by .github/workflows/tracs-image.yml.
#
# TRACS's setup.py hard-codes -march=native, so the bioconda image only runs on
# CPUs at least as new as bioconda's build host and SIGILLs (exit 132) elsewhere.
# Here the flag is patched to a fixed baseline so one image runs on any host.
#   x86-64-v2 = SSE4.2/POPCNT, every x86-64 server CPU since ~2009 (default)
#   x86-64-v3 = AVX2, ~2013+
# setuptools appends extra_compile_args last, so CFLAGS cannot override it.
#
# cxx-compiler is required: without it pybind11 fails late with the misleading
# "Unsupported compiler -- at least C++11 support is needed!". conda-forge's
# compiler activation scripts do not run here, so CC/CXX are set explicitly.
FROM mambaorg/micromamba:1.5.8

ARG TRACS_VERSION=1.1.1
ARG TRACS_MARCH=x86-64-v2

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 pip \
        cxx-compiler pybind11 boost-cpp zlib \
        numpy scipy pandas plotly matplotlib-base pyfastx tqdm joblib \
        sourmash ncbi-genome-download \
        samtools minimap2 htsbox \
        git procps-ng \
    && micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8

RUN export CC=$(ls /opt/conda/bin/*-linux-gnu-gcc | head -1) \
           CXX=$(ls /opt/conda/bin/*-linux-gnu-g++ | head -1) \
    && git clone --depth 1 --branch "v${TRACS_VERSION}" https://github.com/gtonkinhill/tracs.git /tmp/tracs \
    && sed -i "s/-march=native/-march=${TRACS_MARCH}/" /tmp/tracs/setup.py \
    && grep -q -- "-march=${TRACS_MARCH}\"" /tmp/tracs/setup.py \
    && pip install --no-cache-dir --no-build-isolation /tmp/tracs \
    && rm -rf /tmp/tracs \
    && tracs --version \
    && for c in align combine distance cluster build-db; do tracs "$c" --help >/dev/null; done \
    && for t in samtools minimap2 htsbox sourmash ps; do command -v "$t" >/dev/null; done

LABEL org.opencontainers.image.source="https://github.com/julianzaugg/genome_metagenome_all_in_one" \
      org.opencontainers.image.description="TRACS ${TRACS_VERSION} built for -march=${TRACS_MARCH}" \
      tracs.march="${TRACS_MARCH}"

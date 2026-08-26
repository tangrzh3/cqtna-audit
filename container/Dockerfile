# Runnable environment for the CD4 dynamic-eQTL target-nomination audit.
#
# READ CONTAINER.md FIRST. This image reproduces the software VERSIONS the
# analyses ran under. It does not reproduce the MACHINE: the reported numbers
# were produced on Windows 11, x86_64-w64-mingw32, R 4.4.1 with the reference
# BLAS, and this image is Linux with OpenBLAS. Treat digit-level agreement as
# something to be measured here, not assumed.
#
# Build (from the repository root):
#   docker build -f container/Dockerfile -t cqtna-audit:0.3.0 .
# Run:
#   docker run --rm -v "$PWD:/repo" -w /repo cqtna-audit:0.3.0 \
#       python3 step127_audit_manuscript_numbers.py /repo
#
# Bioconductor 3.19 is the release that pairs with R 4.4 and with the Bioc
# package versions in 150a_environment.lock (annotate 1.82.0, org.Hs.eg.db
# 3.19.1). Pinning the base to a RELEASE_ tag rather than :devel is what keeps
# this rebuildable later.
FROM bioconductor/bioconductor_docker:RELEASE_3_19

ARG DEBIAN_FRONTEND=noninteractive

# hdf5 for hdf5r, libxml2 for XML-dependent Bioc packages, libglpk for igraph
# (Seurat), unzip/wget for the two external binaries below.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-dev libxml2-dev libglpk-dev libgsl-dev \
        libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
        libharfbuzz-dev libfribidi-dev libtiff5-dev \
        python3 python3-pip python3-venv \
        wget unzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# --- Python -------------------------------------------------------------
# The analyses ran on CPython 3.12.4. This installs the distribution's 3.12.x,
# which is the same minor version but not necessarily the same patch. See
# CONTAINER.md section 3 for why that is recorded rather than hidden.
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"
COPY container/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r /tmp/requirements.txt

# --- R ------------------------------------------------------------------
COPY 150a_environment.lock /opt/150a_environment.lock
COPY container/install_r_packages.R /tmp/install_r_packages.R
RUN CQTNA_LOCK=/opt/150a_environment.lock Rscript /tmp/install_r_packages.R

# --- External binaries --------------------------------------------------
# SMR 1.3.1 and PLINK 2.0.0-a.7.2, the versions named in the Methods. The
# repository's bin/ holds the Windows builds; these are the Linux ones.
# The upstream URLs are not content-addressed, so verify the checksums against
# the ones you record at build time -- do not assume a URL still serves the
# same bytes it served for us.
RUN mkdir -p /opt/bin \
 && wget -q -O /tmp/smr.zip https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip \
 && unzip -q -j /tmp/smr.zip -d /opt/bin \
 && wget -q -O /tmp/plink2.zip https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20241206.zip \
 && unzip -q -j /tmp/plink2.zip -d /opt/bin \
 && chmod +x /opt/bin/* \
 && rm -f /tmp/smr.zip /tmp/plink2.zip
ENV PATH="/opt/bin:${PATH}"

# Fail fast and loudly rather than half-working.
RUN python3 -c "import numpy, pandas, scipy, pyarrow, matplotlib; print('python deps ok')" \
 && Rscript -e "stopifnot(getRversion() >= '4.4.0'); cat('R', as.character(getRversion()), '\n')"

WORKDIR /repo
CMD ["bash"]

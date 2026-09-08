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
#
# bioconductor_docker:RELEASE_3_19 is Ubuntu 22.04 (jammy), confirmed by
# running the image rather than assumed -- jammy's own repos top out at
# Python 3.10, and requirements.txt pins scipy==1.18.0, which needs >=3.12.
# The deadsnakes PPA supplies it.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-dev libxml2-dev libglpk-dev libgsl-dev \
        libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
        libharfbuzz-dev libfribidi-dev libtiff5-dev \
        wget unzip ca-certificates curl gnupg \
        git \
    && rm -rf /var/lib/apt/lists/*
# add-apt-repository's own key fetch uses the keyserver protocol, which this
# network could not reach (gpg exit status 2, discovered building
# container/Dockerfile.acceptance, the slim image, first). Fetching the same
# key over plain HTTPS from keyserver.ubuntu.com and dearmoring it manually
# is a more reliable transport for the same key. Fingerprint taken from that
# failed run's own log rather than re-derived.
RUN curl -fsSL "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0xBA6932366A755776" \
      | gpg --dearmor -o /etc/apt/trusted.gpg.d/deadsnakes.gpg \
 && echo "deb https://ppa.launchpadcontent.net/deadsnakes/ppa/ubuntu jammy main" \
      > /etc/apt/sources.list.d/deadsnakes.list \
 && apt-get update && apt-get install -y --no-install-recommends \
        python3.12 python3.12-venv python3.12-dev \
    && rm -rf /var/lib/apt/lists/*

# --- Python -------------------------------------------------------------
# The analyses ran on CPython 3.12.4. This installs the distribution's 3.12.x,
# which is the same minor version but not necessarily the same patch. See
# CONTAINER.md section 3 for why that is recorded rather than hidden.
RUN python3.12 -m venv /opt/venv
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
#
# > THE PLINK URL IS DEAD. s3://plink2-assets/alpha7 now holds only 2026 dated
# builds; plink2_linux_x86_64_20241206.zip -- the 2.0.0-a.7.2 the Methods
# names -- has been deleted upstream. Verified by listing the bucket, not
# inferred from the 404. Substituting a current build would put a different
# binary behind a version number the paper states, so nothing is substituted.
#
# The build no longer dies on it. What is obtained is recorded with a checksum
# in /opt/binaries_provenance.txt, and what is not is stated in the same file,
# so the image can say what it holds instead of merely looking complete. No
# analysis on S54's must-run list calls either binary -- step145 names SMR and
# HEIDI in comments but computes in R, and it ran green in the slim acceptance
# image, which carries neither -- so this does not affect the acceptance. It
# does affect a reader trying to rerun everything, which is why it is recorded
# rather than quietly worked around.
RUN mkdir -p /opt/bin \
 && : > /opt/binaries_provenance.txt \
 && if wget -q -O /tmp/smr.zip https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip; then \
      unzip -q -j /tmp/smr.zip -d /opt/bin \
      && echo "smr-1.3.1 OK sha256=$(sha256sum /tmp/smr.zip | cut -c1-64)" >> /opt/binaries_provenance.txt; \
    else \
      echo 'smr-1.3.1 UNOBTAINABLE' >> /opt/binaries_provenance.txt; \
    fi \
 && if wget -q -O /tmp/plink2.zip https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20241206.zip; then \
      unzip -q -j /tmp/plink2.zip -d /opt/bin \
      && echo "plink2-20241206 OK sha256=$(sha256sum /tmp/plink2.zip | cut -c1-64)" >> /opt/binaries_provenance.txt; \
    else \
      echo 'plink2-20241206 UNOBTAINABLE (404; upstream deleted this dated build)' >> /opt/binaries_provenance.txt; \
    fi \
 && chmod +x /opt/bin/* 2>/dev/null || true \
 && rm -f /tmp/smr.zip /tmp/plink2.zip \
 && cat /opt/binaries_provenance.txt
ENV PATH="/opt/bin:${PATH}"

# The package under test, built from source in the repository. Sixteen R
# scripts call library(cqtna), including step140 and step141 on S54's must-run
# list, so without this the acceptance dies at the start of phase 2 -- hours
# into a build that looked like it was succeeding. The slim acceptance image
# had it from the start; this one never did, and nothing caught that because
# this image had never actually been built.
COPY cqtna_r /tmp/cqtna_r
RUN R CMD INSTALL /tmp/cqtna_r

# Fail fast and loudly rather than half-working.
RUN python3 -c "import numpy, pandas, scipy, pyarrow, matplotlib; print('python deps ok')" \
 && Rscript -e "stopifnot(getRversion() >= '4.4.0'); cat('R', as.character(getRversion()), '\n')" \
 && Rscript -e "library(cqtna); cat('cqtna', as.character(packageVersion('cqtna')), '\n')" \
 && git --version

WORKDIR /repo
CMD ["bash"]

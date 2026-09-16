FROM bitnami/minideb:bookworm

RUN apt-get update && apt-get install --no-install-recommends --yes \
    ca-certificates \
    wget \
    curl \
    pkg-config \
    procps \
    git \
    gcc \
    man-db \
    g++ \
    make \
    autoconf \
    unzip \
    bzip2 \
    zlib1g-dev \
    dpkg-dev \
    build-essential \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    libgdbm-dev \
    libc6-dev \
    libsqlite3-dev \
    libssl-dev \
    openssl \
    libffi-dev \
    parallel \
    time \
    cmake \
    libclang-dev \
    grabix

# install htslib
ENV HTSLIB_VERSION="1.20"
RUN wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& rm htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& cd htslib-${HTSLIB_VERSION} \
	&& ./configure --enable-libcurl --enable-gcs --enable-s3 \
	&& make \
	&& make install \
	&& cd / \
	&& rm -rf htslib-${HTSLIB_VERSION}

# install samtools
ENV SAMTOOLS_VERSION="1.20"
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure --enable-libcurl --enable-gcs --enable-s3 \
	&& make \
	&& make install \
	&& cd / \
	&& rm -rf samtools-${SAMTOOLS_VERSION}

# install bedtools
ENV BED_TOOLS_VERSION="2.31.0"
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BED_TOOLS_VERSION}/bedtools.static \
	&& chmod 755 bedtools.static \
	&& mv bedtools.static /usr/bin/bedtools

# install python3 (bookworm marks its system python externally managed, which blocks pip)
RUN apt-get update && apt-get install --no-install-recommends --yes python3 python3-dev python3-distutils \
	&& rm -f /usr/lib/python3*/EXTERNALLY-MANAGED \
	&& curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
	&& python3 get-pip.py \
	&& rm get-pip.py

# install gcloud
ENV GCLOUD_SDK_VERSION="531.0.0"
RUN wget -q https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERSION}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERSION}-linux-x86_64.tar.gz \
    && rm google-cloud-sdk-${GCLOUD_SDK_VERSION}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

# install aws-cli (v1 from bookworm apt) for fast multipart S3 downloads via
# `aws s3 cp --no-sign-request`. Used by the HPRC R2 BAM-localization path.
RUN apt-get update && apt-get install --no-install-recommends --yes awscli \
    && rm -rf /var/lib/apt/lists/*

# install TRF
RUN wget -q https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 \
    && chmod 755 trf409.linux64 \
    && mv trf409.linux64 /usr/bin/trf

RUN python3 -m pip install --no-cache-dir git+https://github.com/bw2/trviz@9068d9abc50ec9b8ba478505d3d19496e1f431b9 pyarrow


# install Rust/cargo, build TRGT from source, then remove Rust/cargo.
#
# TRGT is pinned to an exact commit rather than tracking the fork's main branch, so that two builds
# of this Dockerfile produce the same genotyper and a delivered VCF can be traced to the binary that
# made it. Changing this value busts the cache for the build layer below, which is the point.
#
# The clone is deliberately NOT --depth 1: TRGT's build.rs runs `git describe` through vergen to
# embed the version (e.g. 5.1.0-9171d21), and that needs history and tags. A shallow clone would
# still build, but the binary would report a bare 5.1.0 with no commit, losing the provenance this
# pin exists to provide. The repo is small, so a full clone costs little.
ENV TRGT_COMMIT="9171d21009337bc92343fc148135722338902d9d"
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y \
    && . "$HOME/.cargo/env" \
    && git clone https://github.com/bw2/trgt.git \
    && cd trgt \
    && git checkout ${TRGT_COMMIT} \
    && cargo build --release \
    && mv target/release/trgt /usr/bin/trgt \
    && chmod 755 /usr/bin/trgt \
    && cd .. \
    && rm -rf trgt \
    && rustup self uninstall -y

# install trgt-lps
RUN wget -q https://github.com/PacificBiosciences/trgt-lps/releases/download/0.10.0/trgt-lps-v0.10.0.gz \
    && gunzip trgt-lps-v0.10.0.gz \
    && chmod 755 trgt-lps-v0.10.0 \
    && mv trgt-lps-v0.10.0 /usr/bin/trgt-lps

ENV HTSLIB_CONFIGURE_OPTIONS="--enable-s3 --enable-gcs --enable-libcurl --with-libdeflate"
# str-analysis version is injected at build time (CI passes the tested commit sha; defaults to main
# for a plain `docker build`). Changing the arg value busts the pip-install layer cache.
ARG STR_ANALYSIS_COMMIT=main
RUN python3 -m pip install --upgrade --no-cache-dir git+https://github.com/broadinstitute/str-analysis@${STR_ANALYSIS_COMMIT}
LABEL org.opencontainers.image.revision="${STR_ANALYSIS_COMMIT}"

ENV PYTHONPATH=/:.:$PYTHONPATH

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
    cmake \
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
    time

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

# install python3 (bookworm marks its system python externally managed, which blocks pip).
# python3-dev is needed because some ATaRVa dependencies (eg. pyabpoa) have no prebuilt wheel for
# this python version and are compiled from source.
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

# install ATaRVa (https://github.com/sowpatilab/ATaRVa), pinned so that rebuilds produce the same genotyper.
# This provides the `atarva` command.
ENV ATARVA_VERSION="0.7.1"
RUN python3 -m pip install --no-cache-dir ATaRVa==${ATARVA_VERSION} \
    && atarva --help > /dev/null

ENV HTSLIB_CONFIGURE_OPTIONS="--enable-s3 --enable-gcs --enable-libcurl --with-libdeflate"
# str-analysis version is injected at build time (CI passes the tested commit sha; defaults to main
# for a plain `docker build`). Changing the arg value busts the pip-install layer cache.
ARG STR_ANALYSIS_COMMIT=main
RUN python3 -m pip install --upgrade --no-cache-dir git+https://github.com/broadinstitute/str-analysis@${STR_ANALYSIS_COMMIT}
LABEL org.opencontainers.image.revision="${STR_ANALYSIS_COMMIT}"

ENV PYTHONPATH=/:.:$PYTHONPATH

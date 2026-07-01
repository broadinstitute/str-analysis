FROM bitnami/minideb:bullseye

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
    time

# install htslib (provides bgzip + tabix)
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

# install bcftools (used by ensembletr_pipeline.py to concat/sort/reheader the per-caller VCFs)
ENV BCFTOOLS_VERSION="1.20"
RUN wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& cd bcftools-${BCFTOOLS_VERSION} \
	&& ./configure --enable-libcurl --enable-gcs --enable-s3 \
	&& make \
	&& make install \
	&& cd / \
	&& rm -rf bcftools-${BCFTOOLS_VERSION}

# install python3.9 (python3.9-dev provides Python.h, needed to build EnsembleTR's pgenlib C extension from source)
RUN apt-get update && apt-get install --no-install-recommends --yes python3.9 python3.9-dev python3-distutils \
	&& curl https://bootstrap.pypa.io/pip/3.9/get-pip.py -o get-pip.py \
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

# install EnsembleTR from https://github.com/gymrek-lab/EnsembleTR (provides the `EnsembleTR` console script). Pinned to
# a specific commit for reproducibility; pip also pulls its TRTools/cyvcf2/numpy/networkx dependencies.
ENV ENSEMBLETR_COMMIT="dbda15bcf92c2416fb011d1337ddfead1b4d77b6"
RUN python3 -m pip install --no-cache-dir --upgrade pip crcmod \
	&& python3 -m pip install --no-cache-dir git+https://github.com/gymrek-lab/EnsembleTR@${ENSEMBLETR_COMMIT}

# install str-analysis (provides convert_expansion_hunter_json_to_vcf, convert_ensembletr_vcf_to_expansion_hunter_json,
# and combine_str_json_to_tsv --include-extra-ensembletr-fields, all used by ensembletr_pipeline.py). Pinned to the
# commit that adds those EnsembleTR converters + combine flag, for reproducibility (like the sibling images).
ENV HTSLIB_CONFIGURE_OPTIONS="--enable-s3 --enable-gcs --enable-libcurl --with-libdeflate"
RUN python3 -m pip install --upgrade --no-cache-dir git+https://github.com/broadinstitute/str-analysis@8c59663

ENV PYTHONPATH=/:.:$PYTHONPATH

FROM bitnami/minideb:bullseye

RUN apt-get update && apt-get install --no-install-recommends --yes \
    ca-certificates \
    wget -q \
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
	zlib1g-dev \
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
	&& make clean

# install samtools
ENV SAMTOOLS_VERSION="1.20"
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure --enable-libcurl --enable-gcs --enable-s3 \
	&& make \
	&& make install \
	&& make clean

# install bedtools
ENV BED_TOOLS_VERSION="2.31.0"
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BED_TOOLS_VERSION}/bedtools.static \
	&& chmod 777 bedtools.static \
	&& mv bedtools.static /usr/bin/bedtools

# install python3.9
RUN apt-get update && apt-get install --no-install-recommends --yes python3.9 python3-distutils \
	&& curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
	&& python3 get-pip.py

# install gcloud
ENV GCLOUD_SDK_VERISON="531.0.0"
RUN wget -q https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

# install HipSTR
RUN git clone https://github.com/HipSTR-Tool/HipSTR \
    && cd HipSTR \
    && make \
    && cp /HipSTR/HipSTR /bin/HipSTR \
    && make clean

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache

RUN python3 -m pip install --upgrade --no-cache-dir git+https://github.com/broadinstitute/str-analysis

ENV PYTHONPATH=/:.:$PYTHONPATH

FROM bitnami/minideb:bullseye

RUN apt-get update && apt-get install --no-install-recommends --yes \
    ca-certificates \
    wget -q \
    curl \
    pkg-config \
    procps \
    gcc \
    make \
    autoconf \
    unzip \
    bzip2 \
    zlib1g-dev \
    dpkg-dev \
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

# install python3.9
RUN apt-get update && apt-get install --no-install-recommends --yes python3.9 python3-distutils \
	&& curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
	&& python3 get-pip.py


ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache

RUN python3 -m pip install --no-cache-dir git+https://github.com/broadinstitute/str-analysis

ENV PYTHONPATH=/:.:$PYTHONPATH

FROM bitnami/minideb:stretch

RUN install_packages \
    ca-certificates \
    wget -q \
    curl \
    pkg-config \
    procps \
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
    liblzma-dev

ENV SAMTOOLS_VERSION="1.9"
RUN wget -q https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

# install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

# install python3.7 (adds 800Mb)
RUN install_packages python-smbus libncursesw5-dev libgdbm-dev libc6-dev zlib1g-dev libsqlite3-dev libssl-dev openssl libffi-dev

ENV PYTHON_VERSION="3.7.1"
RUN wget -q https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz \
	&& tar xvf Python-${PYTHON_VERSION}.tar.xz \
	&& rm Python-${PYTHON_VERSION}.tar.xz \
	&& cd Python-${PYTHON_VERSION} \
	&& ./configure --enable-optimizations \
	&& make install \
	&& make clean


# install ExpansionHunterDenovo binary from https://github.com/Illumina/ExpansionHunterDenovo
ENV VERSION=v0.9.0
RUN wget -q https://github.com/Illumina/ExpansionHunterDenovo/releases/download/${VERSION}/ExpansionHunterDenovo-${VERSION}-linux_x86_64.tar.gz \
    && tar xzf ExpansionHunterDenovo-${VERSION}-linux_x86_64.tar.gz \
    && mv ExpansionHunterDenovo-${VERSION}-linux_x86_64/bin/ExpansionHunterDenovo /usr/local/bin/ExpansionHunterDenovo \
    && rm -r ExpansionHunterDenovo-${VERSION}-linux_x86_64 ExpansionHunterDenovo-${VERSION}-linux_x86_64.tar.gz

# install ExpansionHunter binary
ENV EH_VERSION=v5.0.0
RUN wget -q https://github.com/Illumina/ExpansionHunter/releases/download/${EH_VERSION}/ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz \
    && tar xzf ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz \
    && mv /ExpansionHunter-${EH_VERSION}-linux_x86_64/bin/ExpansionHunter /usr/local/bin/ \
    && rm -r ExpansionHunter-${EH_VERSION}-linux_x86_64 ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz

# install REViewer
ENV REVIEWER_VERSION=v0.2.4
RUN wget -q https://github.com/Illumina/REViewer/releases/download/${REVIEWER_VERSION}/REViewer-${REVIEWER_VERSION}-linux_x86_64.gz \
    && gunzip REViewer-${REVIEWER_VERSION}-linux_x86_64.gz \
    && mv /REViewer-${REVIEWER_VERSION}-linux_x86_64 /usr/local/bin/REViewer \
    && chmod 777 /usr/local/bin/REViewer

# install gcloud
ENV GCLOUD_SDK_VERISON="354.0.0"
RUN wget -q https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

# python3 packages
RUN python3 -m pip install --upgrade pip \
	&& python3 -m pip install --no-cache-dir "intervaltree>=3.1.0" "numpy>=1.20.3" "pandas>=1.1.4" "pysam>=0.16.0.1"

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache

RUN python3 -m pip install --no-cache-dir str-analysis==0.7.7 \
    && python3 -m pip freeze \
    && sleep 2


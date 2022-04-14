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
    liblzma-dev \
    parallel \
    time

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

# install gcloud
ENV GCLOUD_SDK_VERISON="365.0.0"
RUN wget -q https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

# python3 packages
RUN python3 -m pip install --upgrade pip \
	&& python3 -m pip install --no-cache-dir "intervaltree>=3.1.0" "numpy>=1.20.3" "pandas>=1.1.4" "pysam>=0.16.0.1"

# install GangSTR from https://github.com/gymreklab/GangSTR
ENV GANGSTR_VERSION="2.5.0-Source"
ENV HTSLIB_CFLAGS="-I./htslib-1.9"
ENV HTSLIB_LIBS=" -L$./htslib-1.9 -lhts"
ENV NLOPT_CFLAGS="-I./nlopt-2.4.2"
ENV NLOPT_LIBS=" -L$./nlopt-2.4.2 -lnlopt"
ENV CPPFLAGS="-std=c++98"

ENV LD_LIBRARY_PATH=/usr/local/lib/

RUN wget -q https://github.com/gymreklab/GangSTR/releases/download/v2.5/GangSTR-${GANGSTR_VERSION}.tar.gz \
    && tar xzf GangSTR-${GANGSTR_VERSION}.tar.gz \
    && rm GangSTR-${GANGSTR_VERSION}.tar.gz

WORKDIR /GangSTR-${GANGSTR_VERSION}

RUN python3 -m pip install --upgrade cmake pip crcmod \
    && install_packages git autoconf automake libtool

RUN mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && cmake --install . \
    && cp GangSTR /usr/local/bin/ \
    && cp thirdparty/htslib/lib/libhts.so.3 /usr/local/lib/ \
    && rm -rf /GangSTR-${GANGSTR_VERSION}

WORKDIR /

# STRLing
ENV STRLING_VERSION="v0.5.1"
RUN curl -L https://github.com/quinlan-lab/STRling/releases/download/${STRLING_VERSION}/strling -o /usr/local/bin/strling \
    && chmod 777 /usr/local/bin/strling

# install ExpansionHunterDenovo binary from https://github.com/Illumina/ExpansionHunterDenovo
ENV EH_DENOVO_VERSION="v0.9.0"
RUN wget -q https://github.com/Illumina/ExpansionHunterDenovo/releases/download/${EH_DENOVO_VERSION}/ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz \
    && tar xzf ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz \
    && mv ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64/bin/ExpansionHunterDenovo /usr/local/bin/ExpansionHunterDenovo \
    && rm -r ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64 ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz

# install ExpansionHunter binary
#ENV EH_VERSION="v5.0.0"
#RUN wget -q https://github.com/Illumina/ExpansionHunter/releases/download/${EH_VERSION}/ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz \
#    && tar xzf ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz \
#    && mv /ExpansionHunter-${EH_VERSION}-linux_x86_64/bin/ExpansionHunter /usr/local/bin/ \
#    && rm -r ExpansionHunter-${EH_VERSION}-linux_x86_64 ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz

RUN wget https://raw.githubusercontent.com/bw2/ExpansionHunter/master/bin/linux/ExpansionHunter -O /usr/local/bin/ExpansionHunter \
    && chmod 777 /usr/local/bin/ExpansionHunter

# install REViewer
ENV REVIEWER_VERSION="v0.2.7"
RUN wget -q https://github.com/Illumina/REViewer/releases/download/${REVIEWER_VERSION}/REViewer-${REVIEWER_VERSION}-linux_x86_64.gz \
    && gunzip REViewer-${REVIEWER_VERSION}-linux_x86_64.gz \
    && mv /REViewer-${REVIEWER_VERSION}-linux_x86_64 /usr/local/bin/REViewer \
    && chmod 777 /usr/local/bin/REViewer

#RUN python3 -m pip install --no-cache-dir str-analysis==0.9.1 \
#    && python3 -m pip freeze \
#    && sleep 2

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache

COPY str_analysis ./str_analysis
COPY README.md .
COPY setup.py .
COPY requirements.txt .
RUN python3 setup.py install


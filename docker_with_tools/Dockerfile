FROM bitnami/minideb:bullseye

RUN apt-get update && apt-get install --no-install-recommends --yes \
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

RUN apt-get update && apt-get install --no-install-recommends --yes \
	libncursesw5-dev \
	libgdbm-dev \
	libc6-dev \
	zlib1g-dev \
	libsqlite3-dev \
	libssl-dev \
	openssl \
	libffi-dev

# install python3.9
RUN apt-get update && apt-get install --no-install-recommends --yes python3.9 python3-distutils \
	&& curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
	&& python3 get-pip.py

# python3 packages
RUN python3 -m pip install --upgrade --no-cache-dir "intervaltree>=3.1.0" "numpy>=1.20.3" "pandas>=1.1.4" "pysam>=0.16.0.1"

# install gcloud
ENV GCLOUD_SDK_VERISON="481.0.0"
RUN wget -q https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

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

RUN python3 -m pip install --no-cache-dir --upgrade cmake pip crcmod \
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

RUN git clone https://github.com/oscarlr-TRs/EH-expansion-evaluator.git

# STRLing
ENV STRLING_VERSION="v0.5.2"
RUN curl -L https://github.com/quinlan-lab/STRling/releases/download/${STRLING_VERSION}/strling -o /usr/local/bin/strling \
    && chmod 777 /usr/local/bin/strling

# install ExpansionHunterDenovo binary from https://github.com/Illumina/ExpansionHunterDenovo
ENV EH_DENOVO_VERSION="v0.9.0"
RUN wget -q https://github.com/Illumina/ExpansionHunterDenovo/releases/download/${EH_DENOVO_VERSION}/ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz \
    && tar xzf ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz \
    && mv ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64/bin/ExpansionHunterDenovo /usr/local/bin/ExpansionHunterDenovo \
    && rm -r ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64 ExpansionHunterDenovo-${EH_DENOVO_VERSION}-linux_x86_64.tar.gz


# install ExpansionHunter fork with optimizations
RUN wget -q https://github.com/bw2/ExpansionHunter/archive/refs/heads/master.zip \
	&& unzip master.zip \
	&& cd ExpansionHunter-master \
	&& mkdir build \
	&& cd build \
	&& cmake .. \
	&& make \
	&& cp install/bin/ExpansionHunter /usr/local/bin/ \
	&& cd / \
	&& rm -r ExpansionHunter-master master.zip


# install REViewer fork
RUN wget -q https://github.com/bw2/REViewer/archive/refs/heads/master.zip \
    && unzip master.zip \
    && cd REViewer-master \
	&& mkdir build \
	&& cd build \
	&& cmake .. \
	&& make \
	&& cp install/bin/REViewer /usr/local/bin/ \
	&& cd / \
	&& rm -r REViewer-master master.zip

# install HipSTR
RUN git clone https://github.com/HipSTR-Tool/HipSTR \
    && cd HipSTR \
    && make \
    && cp /HipSTR/HipSTR /bin/HipSTR \
    && make clean

# install TRGT
ENV TRGT_VERSION="1.0.0"
RUN wget -q https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
    && tar xzf trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
    && chmod 777 trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt \
    && mv trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt /usr/bin/trgt

# install LongTR
RUN git clone https://github.com/gymrek-lab/LongTR.git \
    && cd LongTR \
    && sed -i 's/git@github.com:/https:\/\/github.com\//' Makefile \
    && make \
    && cp ./LongTR /bin/LongTR \
    && make clean


# install jq
RUN apt-get update && apt-get install --no-install-recommends --yes  jq  python3-dev

# install bedtools
ENV BED_TOOLS_VERSION="2.31.0"
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BED_TOOLS_VERSION}/bedtools.static \
	&& chmod 777 bedtools.static \
	&& mv bedtools.static /usr/bin/bedtools

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache

COPY str_analysis /str_analysis
COPY README.md .
COPY setup.py .
COPY requirements.txt .
RUN python3 -m pip install --no-cache-dir -r requirements.txt
RUN python3 setup.py install

ENV PYTHONPATH=/:.:$PYTHONPATH

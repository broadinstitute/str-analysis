FROM bitnami/minideb:bullseye

RUN apt-get update && apt-get install --no-install-recommends --yes \
    ca-certificates \
    wget \
    curl \
    pkg-config \
    procps \
    git \
    gcc \
    g++ \
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

# install python3.9
RUN apt-get update && apt-get install --no-install-recommends --yes python3.9 python3-distutils \
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

# install dependencies for ExpansionHunter build
RUN apt-get update && apt-get install --no-install-recommends --yes \
    build-essential \
    libncursesw5-dev

RUN python3 -m pip install --upgrade cmake

# install official Illumina ExpansionHunter v5.0.0 prebuilt binary as IlluminaExpansionHunter
# (installed before the fork below so this pinned, stable layer stays cached when the fork's master.zip changes)
ENV ILLUMINA_EXPANSION_HUNTER_VERSION="5.0.0"
RUN wget -q https://github.com/Illumina/ExpansionHunter/releases/download/v${ILLUMINA_EXPANSION_HUNTER_VERSION}/ExpansionHunter-v${ILLUMINA_EXPANSION_HUNTER_VERSION}-linux_x86_64.tar.gz \
	&& tar xzf ExpansionHunter-v${ILLUMINA_EXPANSION_HUNTER_VERSION}-linux_x86_64.tar.gz \
	&& cp ExpansionHunter-v${ILLUMINA_EXPANSION_HUNTER_VERSION}-linux_x86_64/bin/ExpansionHunter /usr/local/bin/IlluminaExpansionHunter \
	&& chmod 755 /usr/local/bin/IlluminaExpansionHunter \
	&& rm -r ExpansionHunter-v${ILLUMINA_EXPANSION_HUNTER_VERSION}-linux_x86_64 ExpansionHunter-v${ILLUMINA_EXPANSION_HUNTER_VERSION}-linux_x86_64.tar.gz

# install mimalloc (faster allocator; ~30% faster genotyping with lower peak RSS, byte-identical output).
# Built from source so its CMake config is installed and the ExpansionHunter fork's find_package(mimalloc) links it.
# (before the fork build below so this stable layer stays cached when the fork's master.zip changes)
ENV MIMALLOC_VERSION="2.1.7"
RUN wget -q https://github.com/microsoft/mimalloc/archive/refs/tags/v${MIMALLOC_VERSION}.tar.gz \
	&& tar xzf v${MIMALLOC_VERSION}.tar.gz \
	&& rm v${MIMALLOC_VERSION}.tar.gz \
	&& cd mimalloc-${MIMALLOC_VERSION} \
	&& mkdir build \
	&& cd build \
	&& cmake .. -DCMAKE_BUILD_TYPE=Release -DMI_BUILD_TESTS=OFF \
	&& make -j$(nproc) \
	&& make install \
	&& cd / \
	&& rm -r mimalloc-${MIMALLOC_VERSION}

# install ExpansionHunter fork
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
#RUN wget -q https://github.com/bw2/REViewer/archive/refs/heads/master.zip \
#    && unzip master.zip \
#    && cd REViewer-master \
#	&& mkdir build \
#	&& cd build \
#	&& cmake .. \
#	&& make \
#	&& cp install/bin/REViewer /usr/local/bin/ \
#	&& cd / \
#	&& rm -r REViewer-master master.zip

# install TRF
RUN wget -q https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 \
	&& chmod 755 trf409.linux64 \
	&& mv trf409.linux64 /usr/bin/trf


ENV HTSLIB_CONFIGURE_OPTIONS="--enable-s3 --enable-gcs --enable-libcurl --with-libdeflate"
RUN python3 -m pip install --upgrade --no-cache-dir git+https://github.com/broadinstitute/str-analysis@29d7e90969e83c9405de4bcb509f88932ed7f570

ENV PYTHONPATH=/:.:$PYTHONPATH

FROM ubuntu:18.04

MAINTAINER jmonlong@ucsc.edu

## Misc
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        wget \
        git \
        gcc g++ cmake parallel \
        libxml2-dev libssl-dev libmariadbclient-dev libcurl4-openssl-dev \
        apt-transport-https software-properties-common dirmngr gpg-agent \
        && rm -rf /var/lib/apt/lists/*

## Install pandoc
RUN apt-get update \
    && wget https://github.com/jgm/pandoc/releases/download/1.18/pandoc-1.18-1-amd64.deb \
    && dpkg -i pandoc-1.18-1-amd64.deb \
    && apt-get install -f \
    && rm pandoc-1.18-1-amd64.deb

## Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
        && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
        && apt-get update \
        && DEBIAN_FRONTEND=noninteractive apt-get install -y r-base r-base-dev

WORKDIR /home

## Install R packages
ADD install.R /home

RUN R -f /home/install.R

## Install sveval package
RUN git clone https://github.com/jmonlong/sveval.git /sveval

RUN R -e "devtools::install('/sveval')"

## Install Sniffles
RUN git clone https://github.com/fritzsedlazeck/Sniffles.git /sniffles && \
        cd /sniffles && \
        mkdir -p build/ && \
        cd build/ && \
        cmake .. && \
        make -j

ENV PATH=$PATH:/sniffles/bin/sniffles-core-1.0.12

## Install bcftools, samtools, tabix, snakemake
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        less \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-dev \
        bcftools \
        samtools \
        tabix \
        && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir requests snakemake pandas numpy

## Install R packages I forgot
RUN R -e "install.packages(c('DT'))"

## Install go and indexcov
RUN wget --quiet https://golang.org/dl/go1.15.5.linux-amd64.tar.gz && \
        tar -C /usr/local -xzf go1.15.5.linux-amd64.tar.gz

ENV PATH=$PATH:/usr/local/go/bin

RUN go get -u github.com/brentp/goleft/... && \
        GOBIN=/usr/local/bin/ go install github.com/brentp/goleft/cmd/goleft

# ENV PATH=$PATH:/root/go/bin

## shasta
RUN wget -O /bin/shasta https://github.com/chanzuckerberg/shasta/releases/download/0.7.0/shasta-Linux-0.7.0 && \
        chmod ugo+x /bin/shasta

## minimap2
WORKDIR /build
RUN wget --no-check-certificate https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 \
        && tar -jxvf minimap2-2.17_x64-linux.tar.bz2 \
        && rm minimap2-2.17_x64-linux.tar.bz2

ENV PATH=/build/minimap2-2.17_x64-linux/:$PATH

## python dependencies for the assembly module
RUN pip3 install --no-cache-dir biopython pysam

## SVIM-asm
RUN git clone https://github.com/eldariont/svim-asm.git && \
        cd svim-asm && \
        pip install .

## samplot
RUN git clone https://github.com/ryanlayer/samplot.git && \
        cd samplot && \
        pip install .

## mosdepth
RUN wget --no-check-certificate https://github.com/brentp/mosdepth/releases/download/v0.3.1/mosdepth && \
        chmod +x mosdepth && \
        mv mosdepth /bin/

## IGV soon

## Add scripts
ADD Snakefile /scripts/
ADD extract_tier1_svs.R /scripts/
ADD extract_reads.py /scripts/
ADD Nanopore-Sep2020.conf /scripts/
ADD merge_contigs.py /scripts/
ADD sv-report.Rmd /scripts/
ADD sv_annotation_database.RData /scripts/

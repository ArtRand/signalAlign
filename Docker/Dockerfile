# Basics
FROM gcc:5.4
FROM opentable/anaconda

MAINTAINER Miten Jain, miten@soe.ucsc.edu

# apt-get installs
RUN apt-get update && apt-get install -y git make zlib1g-dev g++
WORKDIR /home/

# signalAlign
RUN git clone --recursive https://github.com/artrand/signalAlign.git
WORKDIR /home/signalAlign/
RUN make

# BWA
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /home/signalAlign/bwa/
RUN make
WORKDIR /home/signalAlign/

# set PATH variables
ENV ROOT /home/signalAlign
ENV PATH "$ROOT/bwa:$PATH"
ENV PATH "$ROOT/bin:$PATH"

# set signalAlign bin as workDir
WORKDIR /home/signalAlign/bin/



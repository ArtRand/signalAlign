# signalMachine dockerfile

# Basics
FROM gcc:5.4

MAINTAINER Art Rand, arand@soe.ucsc.edu

# apt-get installs
RUN apt-get update && apt-get install -y git make zlib1g-dev g++
WORKDIR /home/

# signalAlign
RUN git clone --recursive -b toil_refactor https://github.com/ArtRand/signalAlign.git
WORKDIR /home/signalAlign/
RUN make core

# prepare
RUN mkdir /opt/signalAlign && \
          cp  /home/signalAlign/bin/signalMachine /opt/signalAlign/

COPY wrapper.sh /opt/signalAlign

RUN mkdir /data
WORKDIR /data

ENV PATH "/opt/signalAlign/:$PATH"

ENTRYPOINT ["sh", "/opt/signalAlign/wrapper.sh"]

FROM ubuntu:20.04 AS builder

RUN apt update && apt install -y --no-install-recommends \
  autoconf \
  automake \
  gcc \
  libglib2.0-dev \
  make \
  && rm -rf /var/lib/apt/lists/*
  
WORKDIR /usr/src/app
COPY . .
RUN autoreconf -i \
  && ./configure \
  && make -j \
  && make test \
  && make install \
  && rm -rf /usr/src/app

FROM ubuntu:20.04

RUN apt update && apt install -y --no-install-recommends \
  libglib2.0-0
COPY --from=builder /usr/local/bin /usr/local/bin
COPY --from=builder /usr/local/share/man /usr/local/share/man

ENTRYPOINT ["/usr/local/bin/exonerate"]

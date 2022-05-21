FROM ubuntu:20.04 AS builder

RUN apt update && apt install -y --no-install-recommends \
  autoconf \
  automake \
  ccache \
  gcc \
  libglib2.0-dev \
  make \
  && ln -s $(command -v ccache) /usr/local/bin/gcc \
  && rm -rf /var/lib/apt/lists/*
  
WORKDIR /usr/src/app
COPY . .
RUN --mount=type=cache,target=/ccache \
  export CCACHE_DIR=/ccache \
  && autoreconf -i \
  && ./configure \
  && make -j \
  && make check \
  && make install \
  && rm -rf /usr/src/app

FROM ubuntu:20.04

RUN apt update && apt install -y --no-install-recommends \
  libglib2.0-0
COPY --from=builder /usr/local/bin /usr/local/bin
COPY --from=builder /usr/local/share/man /usr/local/share/man

ENTRYPOINT ["/usr/local/bin/exonerate"]

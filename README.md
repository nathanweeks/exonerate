exonerate
=========

A continuation of exonerate, a generic tool for sequence alignment by Guy St.
C. Slater et al. (http://www.ebi.ac.uk/~guy/exonerate/).

Exonerate is apparently unmaintained. The goal of this fork is mainly bug
fixes to the final version (2.4.0).

To check out/build the final (stable) version:

```
$ git clone https://github.com/nathanweeks/exonerate.git
$ cd exonerate
$ git checkout v2.4.0
$ ./configure [YOUR_CONFIGURE_OPTIONS]
$ make
$ make check
$ sudo make install
```

**NOTE**: do not use the `exonerate -c/--cores` option, as the multi-threading code contains data races.
Multi-threading has been disabled after v2.4.0.

To build the master branch:

```
$ git clone https://github.com/nathanweeks/exonerate.git
$ cd exonerate
$ autoreconf -i
$ ./configure [YOUR_CONFIGURE_OPTIONS]
$ make
$ make check
$ sudo make install
```

This repository contains a Dockerfile that formalizes the process to build and install exonerate on Ubuntu.
To build a [Docker](https://www.docker.com/) image and execute exonerate (assuming input sequences are in the current working directory or a subdirectory of on the host):

```
$ git clone https://github.com/nathanweeks/exonerate.git
$ cd exonerate
$ docker build -t exonerate .
$ docker run -it --rm -v $PWD:/mnt -w /mnt exonerate --version
exonerate from exonerate version 2.4.0
Using glib version 2.64.3
Built on Jan 23 2021
```

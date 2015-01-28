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
$ make install
```

To build the development version (at your own risk!):

```
$ git clone https://github.com/nathanweeks/exonerate.git
$ cd exonerate
$ git checkout develop
$ autoreconf -i
$ ./configure [YOUR_CONFIGURE_OPTIONS]
$ make
$ make check
$ make install
```

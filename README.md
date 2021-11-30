# The Linbox Library

[![Build Status](https://ci.inria.fr/linbox/buildStatus/icon?job=LinBox)](https://ci.inria.fr/linbox/job/LinBox/)

## Purpose

The Linbox library provides functionality for exact linear algebra.
See doc/mainpage.doxy for more info.

## Quick Install

Download [linbox-auto-install.sh](https://github.com/linbox-team/linbox/raw/master/linbox-auto-install.sh) and run it!

For instance, on a machine with an installation of GMP and OpenBLAS in the standard search paths:
```
./linbox-auto-install.sh --stable=yes --make-flags="-j 3" --with-blas-libs="-lopenblas"
```
This script will install stable versions of Givaro, fflas-ffpack and then LinBox.

Here is another example fetching and installing the latest versions of GMP, Givaro, OpenBLAS, fflas-ffpack and then LinBox.
```
./linbox-auto-install.sh --enable-openblas=yes --enable-gmp=yes
```

## Installation

In brief: ```./configure <options> && make && make install```

See INSTALL for generic installation information.

## Availability

From github.com/linbox-team


## Requirements
- fflas-ffpack
and by transitivity:
- any BLAS (Fortran or C): e.g. ATLAS, OpenBLAS, ...
- Givaro
- GMP



## Optional Dependencies 
- NTL, 
- IML, 
- FLINT, 

See  doc/install*html for details.

This library requires the GNU C++ compiler (gcc-4.3 or newer) or any 
compiler supporting advanced template features.

## AUTHORS

The LinBox group (see AUTHORS file for a list of contributors).

## Citing LinBox

If your research depends on the LinBox library, please consider citing the project as

```
@manual{linbox,
title = {{LinBox}},
author = {The {LinBox} group},
edition = {v1.6.3},
year = {2019},
url = {http://github.com/linbox-team/linbox}
}
```
## Contact and discussions

Corrections, suggestions and comments to linbox-use@googlegroups.com

 


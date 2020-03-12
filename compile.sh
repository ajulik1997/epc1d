#!/bin/bash

PYVER=3.8
INCLUDE_PYLIB=$(python$PYVER -c "from distutils import sysconfig; print(sysconfig.get_python_inc())")
LIBPYTHON=$(python$PYVER -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
INCLUDE_NUMPY=$(python$PYVER -c "import numpy; print(numpy.get_include())")

cython -3 -a -v -f --embed -o epc1d.c epc1d.pyx
gcc epc1d.c -pthread -lm -lutil -ldl -fPIC -fwrapv -O3 -Wall -fno-strict-aliasing -I$INCLUDE_PYLIB -I$INCLUDE_NUMPY -L$LIBPYTHON -lpython$PYVER -oepc1d

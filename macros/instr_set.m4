# Copyright(c)'1994-2017 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
#
# Author Clement Pernet
#/

AC_DEFUN([INSTR_SET],
[
        SIMD_CFLAGS=""
        AS_ECHO("Detecting SIMD instruction set")

        AC_ARG_ENABLE(sse,[AC_HELP_STRING([--disable-sse], [ disable SSE instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(sse2,[AC_HELP_STRING([--disable-sse2], [ disable SSE2 instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(sse3,[AC_HELP_STRING([--disable-sse3], [ disable SSE3 instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(ssse3,[AC_HELP_STRING([--disable-ssse3], [ disable SSSE3 instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(sse41,[AC_HELP_STRING([--disable-sse41], [ disable SSE4.1 instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(sse42,[AC_HELP_STRING([--disable-sse42], [ disable SSE4.2 instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(avx,[AC_HELP_STRING([--disable-avx], [ disable AVX instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(avx2,[AC_HELP_STRING([--disable-avx2], [ disable AVX2 instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(fma,[AC_HELP_STRING([--disable-fma], [ disable FMA instruction set (enabled by default when available)])],[],[])
        AC_ARG_ENABLE(fma4,[AC_HELP_STRING([--disable-fma4], [ disable FMA4 instruction set (enabled by default when available)])],[],[])

        AC_TRY_RUN([
                        #include "macros/CodeChunk/instrset_detect.cpp"
                        // increment by one to distinguish from compilation failure error code
                        int main(){return instrset_detect()+1;}
                ],[AS_ECHO("Using 80386 instruction set")],[
                iset=$?
                AS_IF([ test "$iset" -ge "2" -a "x$enable_sse" != "xno" ], [
                        AS_ECHO("SSE enabled")
                        SIMD_CFLAGS="${SIMD_CFLAGS} -msse"
                        HAVE_SSE="yes"
                ],[AS_ECHO("SSE disabled")])
                AS_IF([ test "$iset" -ge "3" -a "x$enable_sse2" != "xno" ], [
                        AS_ECHO("SSE2 enabled")
                        SIMD_CFLAGS="${SIMD_CFLAGS} -msse2"
                ],[AS_ECHO("SSE2 disabled")])
                AS_IF([ test "$iset" -ge "4" -a "x$enable_sse3" != "xno" ], [
                        AS_ECHO("SSE3 enabled")
                        SIMD_CFLAGS="${SIMD_CFLAGS} -msse3"
                ],[AS_ECHO("SSE3 disabled")])
                AS_IF([ test "$iset" -ge "5" -a "x$enable_ssse3" != "xno" ], [
                        AS_ECHO("SSSE3 enabled")
                        SIMD_CFLAGS="${SIMD_CFLAGS} -mssse3"
                ],[AS_ECHO("SSSE3 disabled")])
                AS_IF([ test "$iset" -ge "6" -a "x$enable_sse41" != "xno" ], [
                        AS_ECHO("SSE4.1 enabled")
                        SIMD_CFLAGS="${SIMD_CFLAGS} -msse4.1"
                ],[AS_ECHO("SSE4.1 disabled")])
                AS_IF([ test "$iset" -ge "7" -a "x$enable_sse42" != "xno" ], [
                        AS_ECHO("SSE4.2 enabled")
                        SIMD_CFLAGS="${SIMD_CFLAGS} -msse4.2"
                ],[AS_ECHO("SSE4.2 disabled")])
                AS_IF([ test "$iset" -ge "8" -a "x$enable_avx" != "xno" ], [
                        AS_ECHO("AVX enabled")
                        SIMD_CFLAGS="${SIMD_CFLAGS} -mavx"
                ],[AS_ECHO("AVX disabled")])
                AS_IF([ test "$iset" -ge "9" -a "x$enable_avx2" != "xno" ], [
                        AS_ECHO("AVX2 enabled")
                        SIMD_CFLAGS="${SIMD_CFLAGS} -mavx2"
                ],[AS_ECHO("AVX2 disabled")])
        ])
        AC_TRY_RUN([
                        #include "macros/CodeChunk/instrset_detect.cpp"
                        int main(){return !hasFMA3();}
                   ],[
                        AS_IF([ test "x$enable_fma" != "xno" ], [
                                AS_ECHO("FMA3 enabled")
                                SIMD_CFLAGS="${SIMD_CFLAGS} -mfma"
                              ],[AS_ECHO("FMA3 disabled")])
                ],[AS_ECHO("FMA3 disabled")])
        AC_TRY_RUN([
                        #include "macros/CodeChunk/instrset_detect.cpp"
                        int main(){return !hasFMA4();}
                   ],[
                        AS_IF([ test "x$enable_fma4" != "xno" ], [
                                AS_ECHO("FMA4 enabled")
                                SIMD_CFLAGS="${SIMD_CFLAGS} -mfma4"
                              ],[AS_ECHO("FMA4 disabled")])
                ],[AS_ECHO("FMA4 disabled")])
])

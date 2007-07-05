INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/gmp++ ${GMP_VERSION}

AM_CFLAGS=-Wall -O2

noinst_LTLIBRARIES=libgmpxx.la

libgmpxx_la_SOURCES=\
	gmp++_int_mul.C		gmp++_int_mod.C		\
	gmp++_int_cstor.C	gmp++_int_add.C		\
	gmp++_int_div.C		gmp++_int_io.C		\
	gmp++_int_pow.C		gmp++_int_gcd.C		\
	gmp++_int_misc.C	gmp++_int_sub.C		\
	gmp++_int_compare.C	
#gmp++_int.C

pkgincludesubdir=$(pkgincludedir)/util/gmp++

pkgincludesub_HEADERS=\
	gmp++_int_mul.C		gmp++_int_mod.C		\
	gmp++_int_cstor.C	gmp++_int_add.C		\
	gmp++_int_div.C		gmp++_int_io.C		\
	gmp++_int_pow.C		gmp++_int_gcd.C		\
	gmp++_int_misc.C	gmp++_int_sub.C		\
	gmp++_int_compare.C	gmp++_int.C

## Process this file with automake to produce Makefile.in

EXTRA_DIST = 

pkgdata_DATA =

SUBDIRS = macros intl po 

INCLUDES = $(GTK_CFLAGS) -Wall -ansi	\
	$(GNOME_INCLUDEDIR)

#sysdir   = $(datadir)/gnome/apps/Applications
#sys_DATA = 

LDADD =	$(GNOME_XML_LIB)	\
	$(GNOMEUI_LIBS)		\
	$(GTK_LIBS)		\
	$(INTLLIBS)

bin_PROGRAMS = erdos

erdos_SOURCES = 		\
	color.c			\
	erdos.c			\
	erdos-expression.c	\
	erdos-factory.c		\
	erdos-font.c		\
	erdos-mfrac.c		\
	erdos-mn.c		\
	erdos-mo.c		\
	erdos-mrow.c		\
	erdos-mstyle.c		\
	erdos-presentation.c	\
	erdos-renderer.c	\
	erdos-renderer-x.c	\
	erdos-style.c		\
	erdos-utils.c		\
	item-erdos.c		\
	main.c			\
	color.h			\
	erdos.h			\
	erdos-expression.h	\
	erdos-factory.h		\
	erdos-font.h		\
	erdos-mfrac.h		\
	erdos-mn.h		\
	erdos-mo.h		\
	erdos-mrow.h		\
	erdos-mstyle.h		\
	erdos-presentation.h	\
	erdos-renderer.h	\
	erdos-renderer-x.h	\
	erdos-style.h		\
	erdos-utils.h		\
	item-erdos.h		\
	main.h

erdos_LDFLAGS = $(GNOME_LIBDIR)

CPPFLAGS = \
	-DGNOMELOCALEDIR=\""$(datadir)/locale"\" 	\
	-DBINDIR=\""$(bindir)"\"	\
	-DDATADIR=\""$(datadir)"\"

## to automatically rebuild aclocal.m4 if any of the macros in
## `macros/' change
@MAINT@include macros/macros.dep
@MAINT@macros/macros.dep: macros/Makefile.am
@MAINT@	cd macros && $(MAKE) macros.dep

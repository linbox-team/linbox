oafdir = $(datadir)/oaf
oaf_DATA = equation-editor.oafinfo

EXTRA_DIST = $(oaf_DATA)

INCLUDES = \
	-DGNOMELOCALEDIR=\""$(datadir)/locale"\" 		\
	-DGLADEDIR=\""$(datadir)/equation-editor-data"\"	\
	-I$(top_srcdir)/intl -I$(top_builddir)/intl		\
	-I$(includedir) $(GNOME_INCLUDEDIR)			\
	$(LIBGLADE_CFLAGS)					\
	-DPIXMAPS_DIR=\""$(appicondir)"\"			\
	$(GDK_PIXBUF_CFLAGS)					\
	$(EXTRA_GNOME_CFLAGS)

bin_PROGRAMS = equation-editor

equation_editor_SOURCES = \
	main.c \
	math-expression.c	math-expression.h		\
	math-object.c		math-object.h			\
	identifier.c		identifier.h			\
	number.c		number.h			\
	symbol.c		symbol.h			\
	block.c			block.h				\
	row-block.c		row-block.h			\
	fraction-block.c	fraction-block.h		\
	matrix-block.c		matrix-block.h			\
	math-expression-view.c	math-expression-view.h		\
	layout.c		layout.h			\
	block-layout.c		block-layout.h			\
	glyph-layout.c		glyph-layout.h			\
	row-block-layout.c	row-block-layout.h		\
	fraction-block-layout.c	fraction-block-layout.h		\
	matrix-block-layout.c	matrix-block-layout.h		\
	renderer.c		renderer.h			\
	canvas-renderer.c	canvas-renderer.h		\
	print-renderer.c	print-renderer.h		\
	math-unit.c		math-unit.h			\
	unit-layout.c		unit-layout.h			\
	grouping.c		grouping.h			

equation_editor_LDADD = \
	$(GNOME_LIBDIR)					\
	$(GNOMEUI_LIBS)					\
	$(INTLLIBS)					\
	$(LIBGLADE_LIBS)				\
	$(GDK_PIXBUF_LIBS)				\
	$(EXTRA_GNOME_LIBS)

noinst_PROGRAMS = test

test_SOURCES = \
	test.c							\
	math-expression.c	math-expression.h		\
	math-object.c		math-object.h			\
	identifier.c		identifier.h			\
	number.c		number.h			\
	symbol.c		symbol.h			\
	block.c			block.h				\
	row-block.c		row-block.h			\
	fraction-block.c	fraction-block.h		\
	matrix-block.c		matrix-block.h			\
	math-expression-view.c	math-expression-view.h		\
	layout.c		layout.h			\
	block-layout.c		block-layout.h			\
	glyph-layout.c		glyph-layout.h			\
	row-block-layout.c	row-block-layout.h		\
	fraction-block-layout.c	fraction-block-layout.h		\
	matrix-block-layout.c	matrix-block-layout.h		\
	renderer.c		renderer.h			\
	canvas-renderer.c	canvas-renderer.h		\
	print-renderer.c	print-renderer.h		\
	math-unit.c		math-unit.h			\
	unit-layout.c		unit-layout.h			\
	grouping.c		grouping.h			\

test_LDADD = \
	$(GNOME_LIBDIR)					\
	$(GNOMEUI_LIBS)					\
	$(INTLLIBS)					\
	$(LIBGLADE_LIBS)				\
	$(GDK_PIXBUF_LIBS)

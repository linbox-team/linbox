/* -*- mode: c; style: linux -*- */

/* glyph-layout.c
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "glyph-layout.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _GlyphLayoutPrivate 
{
	/* Private data members */
};

static LayoutClass *parent_class;

static void glyph_layout_init        (GlyphLayout *glyph_layout);
static void glyph_layout_class_init  (GlyphLayoutClass *class);

static void glyph_layout_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void glyph_layout_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void glyph_layout_finalize    (GtkObject *object);

guint
glyph_layout_get_type (void)
{
	static guint glyph_layout_type = 0;

	if (!glyph_layout_type) {
		GtkTypeInfo glyph_layout_info = {
			"GlyphLayout",
			sizeof (GlyphLayout),
			sizeof (GlyphLayoutClass),
			(GtkClassInitFunc) glyph_layout_class_init,
			(GtkObjectInitFunc) glyph_layout_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		glyph_layout_type = 
			gtk_type_unique (layout_get_type (), 
					 &glyph_layout_info);
	}

	return glyph_layout_type;
}

static void
glyph_layout_init (GlyphLayout *glyph_layout)
{
	glyph_layout->p = g_new0 (GlyphLayoutPrivate, 1);
}

static void
glyph_layout_class_init (GlyphLayoutClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("GlyphLayout::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = glyph_layout_finalize;
	object_class->set_arg = glyph_layout_set_arg;
	object_class->get_arg = glyph_layout_get_arg;

	parent_class = LAYOUT_CLASS
		(gtk_type_class (layout_get_type ()));
}

static void
glyph_layout_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	GlyphLayout *glyph_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GLYPH_LAYOUT (object));

	glyph_layout = GLYPH_LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
glyph_layout_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	GlyphLayout *glyph_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GLYPH_LAYOUT (object));

	glyph_layout = GLYPH_LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
glyph_layout_finalize (GtkObject *object) 
{
	GlyphLayout *glyph_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GLYPH_LAYOUT (object));

	glyph_layout = GLYPH_LAYOUT (object);

	g_free (glyph_layout->p);
}

GtkObject *
glyph_layout_new (void) 
{
	return gtk_object_new (glyph_layout_get_type (),
			       NULL);
}

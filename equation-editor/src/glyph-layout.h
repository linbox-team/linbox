/* -*- mode: c; style: linux -*- */

/* glyph-layout.h
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

#ifndef __GLYPH_LAYOUT_H
#define __GLYPH_LAYOUT_H

#include <gnome.h>

#include "layout.h"

BEGIN_GNOME_DECLS

#define GLYPH_LAYOUT(obj)          GTK_CHECK_CAST (obj, glyph_layout_get_type (), GlyphLayout)
#define GLYPH_LAYOUT_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, glyph_layout_get_type (), GlyphLayoutClass)
#define IS_GLYPH_LAYOUT(obj)       GTK_CHECK_TYPE (obj, glyph_layout_get_type ())

typedef struct _GlyphLayout GlyphLayout;
typedef struct _GlyphLayoutClass GlyphLayoutClass;

struct _GlyphLayout 
{
	Layout parent;
};

struct _GlyphLayoutClass 
{
	LayoutClass layout_class;
};

guint glyph_layout_get_type         (void);

GtkObject *glyph_layout_new         (void);

END_GNOME_DECLS

#endif /* __GLYPH_LAYOUT_H */

/* -*- mode: c; style: linux -*- */

/* fraction-block-layout.h
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

#ifndef __FRACTION_BLOCK_LAYOUT_H
#define __FRACTION_BLOCK_LAYOUT_H

#include <gnome.h>

#include "block-layout.h"

BEGIN_GNOME_DECLS

#define FRACTION_BLOCK_LAYOUT(obj)          GTK_CHECK_CAST (obj, fraction_block_layout_get_type (), FractionBlockLayout)
#define FRACTION_BLOCK_LAYOUT_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, fraction_block_layout_get_type (), FractionBlockLayoutClass)
#define IS_FRACTION_BLOCK_LAYOUT(obj)       GTK_CHECK_TYPE (obj, fraction_block_layout_get_type ())

typedef struct _FractionBlockLayout FractionBlockLayout;
typedef struct _FractionBlockLayoutClass FractionBlockLayoutClass;

struct _FractionBlockLayout 
{
	BlockLayout parent;
};

struct _FractionBlockLayoutClass 
{
	BlockLayoutClass block_layout_class;
};

guint fraction_block_layout_get_type         (void);

GtkObject *fraction_block_layout_new         (void);

END_GNOME_DECLS

#endif /* __FRACTION_BLOCK_LAYOUT_H */

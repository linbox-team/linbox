/* -*- mode: c; style: linux -*- */

/* print-renderer.h
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

#ifndef __PRINT_RENDERER_H
#define __PRINT_RENDERER_H

#include <gnome.h>

#include "renderer.h"

BEGIN_GNOME_DECLS

#define PRINT_RENDERER(obj)          GTK_CHECK_CAST (obj, print_renderer_get_type (), PrintRenderer)
#define PRINT_RENDERER_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, print_renderer_get_type (), PrintRendererClass)
#define IS_PRINT_RENDERER(obj)       GTK_CHECK_TYPE (obj, print_renderer_get_type ())

typedef struct _PrintRenderer PrintRenderer;
typedef struct _PrintRendererClass PrintRendererClass;
typedef struct _PrintRendererPrivate PrintRendererPrivate;

struct _PrintRenderer 
{
	Renderer parent;

	PrintRendererPrivate *p;
};

struct _PrintRendererClass 
{
	RendererClass renderer_class;
};

guint print_renderer_get_type         (void);

GtkObject *print_renderer_new         (void);

END_GNOME_DECLS

#endif /* __PRINT_RENDERER_H */

/* -*- mode: c; style: linux -*- */

/* canvas-renderer.h
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

#ifndef __CANVAS_RENDERER_H
#define __CANVAS_RENDERER_H

#include <gnome.h>

#include "renderer.h"

BEGIN_GNOME_DECLS

#define CANVAS_RENDERER(obj)          GTK_CHECK_CAST (obj, canvas_renderer_get_type (), CanvasRenderer)
#define CANVAS_RENDERER_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, canvas_renderer_get_type (), CanvasRendererClass)
#define IS_CANVAS_RENDERER(obj)       GTK_CHECK_TYPE (obj, canvas_renderer_get_type ())

typedef struct _CanvasRenderer CanvasRenderer;
typedef struct _CanvasRendererClass CanvasRendererClass;
typedef struct _CanvasRendererPrivate CanvasRendererPrivate;

struct _CanvasRenderer 
{
	Renderer parent;

	CanvasRendererPrivate *p;
};

struct _CanvasRendererClass 
{
	RendererClass renderer_class;
};

guint canvas_renderer_get_type         (void);

GtkObject *canvas_renderer_new         (void);

END_GNOME_DECLS

#endif /* __CANVAS_RENDERER_H */

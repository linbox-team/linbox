/* -*- mode: c; style: linux -*- */

/* renderer.h
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

#ifndef __RENDERER_H
#define __RENDERER_H

#include <gnome.h>


BEGIN_GNOME_DECLS

#define RENDERER(obj)          GTK_CHECK_CAST (obj, renderer_get_type (), Renderer)
#define RENDERER_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, renderer_get_type (), RendererClass)
#define IS_RENDERER(obj)       GTK_CHECK_TYPE (obj, renderer_get_type ())

typedef struct _Renderer Renderer;
typedef struct _RendererClass RendererClass;
typedef struct _RendererPrivate RendererPrivate;

struct _Renderer 
{
	GtkObject parent;

	RendererPrivate *p;
};

struct _RendererClass 
{
	GtkObjectClass gtk_object_class;
};

guint renderer_get_type         (void);

GtkObject *renderer_new         (void);

END_GNOME_DECLS

#endif /* __RENDERER_H */

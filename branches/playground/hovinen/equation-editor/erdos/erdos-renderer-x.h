/* erdos-renderer.h
 * Copyright (C) 1999  Chris Lahey <clahey@umich.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#ifndef __ERDOS_RENDERER_X_H__
#define __ERDOS_RENDERER_X_H__

#include "erdos-renderer.h"
#include <gdk/gdk.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

typedef struct _ErdosRendererX       ErdosRendererX;
typedef struct _ErdosRendererXClass  ErdosRendererXClass;

struct _ErdosRendererX
{
  ErdosRenderer renderer;
  /* Put your own, widget-specific fields here */
  GdkDrawable *drawable;
  GdkGC *gc;
};

struct _ErdosRendererXClass
{
  ErdosRendererClass parent_class;
  /* Signals go here */
};

void erdos_renderer_x_class_init	(ErdosRendererXClass	 *klass);
ErdosObjectClass *erdos_renderer_x_get_class(void);
ErdosRenderer *erdos_renderer_x_new(GdkDrawable *drawable, GdkGC *gc);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ERDOS_RENDERER_X_H__ */

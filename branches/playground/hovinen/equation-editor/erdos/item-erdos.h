/* item-erdos.h
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
#ifndef __ITEM_ERDOS_H__
#define __ITEM_ERDOS_H__

#include <gnome.h>
#include "erdos-expression.h"
#include "erdos-renderer-x.h"

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

#define ITEM_TYPE_ERDOS			(item_erdos_get_type ())
#define ITEM_ERDOS(obj)			(GTK_CHECK_CAST ((obj), ITEM_TYPE_ERDOS, ItemErdos))
#define ITEM_ERDOS_CLASS(klass)		(GTK_CHECK_CLASS_CAST ((klass), ITEM_TYPE_ERDOS, ItemErdosClass))
#define ITEM_IS_ERDOS(obj)			(GTK_CHECK_TYPE ((obj), ITEM_TYPE_ERDOS))
#define ITEM_IS_ERDOS_CLASS(klass)		(GTK_CHECK_CLASS_TYPE ((obj), ITEM_TYPE_ERDOS))


typedef struct _ItemErdos       ItemErdos;
typedef struct _ItemErdosClass  ItemErdosClass;

struct _ItemErdos
{
  GnomeCanvasItem parent;

  /* Put your own, widget-specific fields here */
  ErdosExpression *expression;
  ErdosStyle *style;
  ErdosRenderer *renderer;
  double height;
  double ascent;
  GdkGC *gc;			/* GC for drawing points */
  GdkGC *bg_gc;
  
  GdkColor   background;
  GdkColor   foreground;
};

struct _ItemErdosClass
{
  GnomeCanvasItemClass parent_class;
};


GtkType    item_erdos_get_type (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ITEM_ERDOS_H__ */

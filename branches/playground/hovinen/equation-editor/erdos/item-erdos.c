/* item-erdos.c
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

#include <gnome.h>
#include "item-erdos.h"
#include "color.h"
static void item_erdos_init		(ItemErdos		 *erdos);
static void item_erdos_class_init	(ItemErdosClass	 *klass);
static void item_erdos_update (GnomeCanvasItem *item, double *affine, ArtSVP *clip_path, int flags);
static void item_erdos_draw (GnomeCanvasItem *item, GdkDrawable *drawable, int x, int y, int width, int height);
static double item_erdos_point (GnomeCanvasItem *item, double x, double y, int cx, int cy, GnomeCanvasItem **actual_item);
static void item_erdos_translate (GnomeCanvasItem *item, double dx, double dy);
static void item_erdos_set_arg (GtkObject *o, GtkArg *arg, guint arg_id);
static void item_erdos_get_arg (GtkObject *object, GtkArg *arg, guint arg_id);
static void item_erdos_realize (GnomeCanvasItem *item);
static void item_erdos_unrealize (GnomeCanvasItem *item);


static GnomeCanvasItemClass *parent_class = NULL;


/* The arguments we take */
enum {
	ARG_0,
	ARG_EXPRESSION,
	/*	ARG_HEIGHT*/
};

GtkType
item_erdos_get_type (void)
{
  static GtkType erdos_type = 0;

  if (!erdos_type)
    {
      static const GtkTypeInfo erdos_info =
      {
        "ItemErdos",
        sizeof (ItemErdos),
        sizeof (ItemErdosClass),
        (GtkClassInitFunc) item_erdos_class_init,
        (GtkObjectInitFunc) item_erdos_init,
        /* reserved_1 */ NULL,
        /* reserved_2 */ NULL,
        (GtkClassInitFunc) NULL,
      };

      erdos_type = gtk_type_unique (gnome_canvas_item_get_type (), &erdos_info);
    }

  return erdos_type;
}

static void
item_erdos_class_init (ItemErdosClass *klass)
{
  GtkObjectClass *object_class;
  GnomeCanvasItemClass *item_class;

  object_class = (GtkObjectClass*) klass;
  item_class = (GnomeCanvasItemClass *) klass;

  parent_class = gtk_type_class (gnome_canvas_item_get_type ());
  
  gtk_object_add_arg_type ("ItemErdos::expression", GTK_TYPE_POINTER, 
			   GTK_ARG_WRITABLE, ARG_EXPRESSION); 
  /*  gtk_object_add_arg_type ("ItemErdos::height", GTK_TYPE_DOUBLE, 
      GTK_ARG_READWRITE, ARG_HEIGHT);*/
 
  object_class->set_arg = item_erdos_set_arg;
  object_class->get_arg = item_erdos_get_arg;
  /*  object_class->destroy = item_erdos_destroy; */

  /* GnomeCanvasItem method overrides */
  item_class->update      = item_erdos_update;
  item_class->realize     = item_erdos_realize;
  item_class->unrealize   = item_erdos_unrealize;
  item_class->draw        = item_erdos_draw;
  item_class->point       = item_erdos_point;
  item_class->translate   = item_erdos_translate;
}

static void
item_erdos_init (ItemErdos *erdos)
{
  erdos->expression = NULL;
}

static void
item_erdos_update (GnomeCanvasItem *item, double *affine, ArtSVP *clip_path, int flags)
{
  ItemErdos *item_erdos = ITEM_ERDOS (item);
  if (GNOME_CANVAS_ITEM_CLASS (parent_class)->update)
    (* GNOME_CANVAS_ITEM_CLASS (parent_class)->update) (item, affine, clip_path, flags);
  
  gnome_canvas_group_child_bounds (GNOME_CANVAS_GROUP (item->parent), item);
}

static void
item_erdos_draw (GnomeCanvasItem *item, GdkDrawable *drawable, int x, int y, int width, int height)
{
  ItemErdos *item_erdos = ITEM_ERDOS (item);
  ((ErdosRendererX *)item_erdos->renderer)->drawable = drawable;
  ((ErdosRendererX *)item_erdos->renderer)->gc = item_erdos->gc;
  gdk_draw_rectangle( drawable, item_erdos->bg_gc, TRUE, 0, 0, width, height );
  erdos_expression_draw( item_erdos->expression, item_erdos->renderer, -x, -y + item_erdos->ascent );
}

static double
item_erdos_point (GnomeCanvasItem *item, double x, double y, int cx, int cy,
		 GnomeCanvasItem **actual_item)
{
  *actual_item = item;
  return 0.0;
}

static void
item_erdos_translate (GnomeCanvasItem *item, double dx, double dy)
{
  printf ("item_erdos_translate %g, %g\n", dx, dy);
}

static void
item_erdos_set_arg (GtkObject *o, GtkArg *arg, guint arg_id)
{
	GnomeCanvasItem *item;
	ItemErdos *item_erdos;

	item = GNOME_CANVAS_ITEM (o);
	item_erdos = ITEM_ERDOS (o);
	
	switch (arg_id){
	case ARG_EXPRESSION:
	  item_erdos->expression = GTK_VALUE_POINTER (*arg);
	  gnome_canvas_item_request_update (item);
	  break;
	  /*	case ARG_HEIGHT:
		item_erdos->height = GTK_VALUE_DOUBLE (*arg);
		break*/
	}
}

static void
item_erdos_get_arg (GtkObject *object, GtkArg *arg, guint arg_id)
{
	ItemErdos *item_erdos;

	item_erdos = ITEM_ERDOS (object);

	switch (arg_id) {
	  /*	case ARG_WIDTH:
		GTK_VALUE_DOUBLE (*arg) = item_erdos->width;
		break;*/
	default:
		arg->type = GTK_TYPE_INVALID;
		break;
	}
}

static void
item_erdos_realize (GnomeCanvasItem *item)
{
        double ascent, descent;
	ItemErdos *item_erdos;

	item_erdos = ITEM_ERDOS (item);

	if (parent_class->realize)
	  (* parent_class->realize) (item);

	if (!item->canvas->aa)
	  {
	    item_erdos->gc = gdk_gc_new (item->canvas->layout.bin_window);
	    item_erdos->bg_gc = gdk_gc_new (item->canvas->layout.bin_window);
	    
	    /* Allocate the default colors */
	    item_erdos->background = ia_white;
	    item_erdos->foreground = ia_black;
	    
	    gdk_gc_set_foreground (item_erdos->gc, &item_erdos->foreground);
	    gdk_gc_set_background (item_erdos->gc, &item_erdos->background);
	    
	    gdk_gc_set_foreground (item_erdos->bg_gc, &item_erdos->background);
	    gdk_gc_set_background (item_erdos->bg_gc, &item_erdos->background);

	    item_erdos->renderer = erdos_renderer_x_new(NULL, item_erdos->gc);
	    erdos_object_ref( (ErdosObject *) item_erdos->renderer );
	    item_erdos->style = erdos_style_new();
	    erdos_style_set_font( item_erdos->style, erdos_font_new() );
	    item_erdos->style->font->fontsize = 40;
	    item_erdos->style->font->fontfamily = "adobe-times";
	    erdos_object_ref( (ErdosObject *) item_erdos->style );
	    erdos_expression_set_style(item_erdos->expression, item_erdos->style);
	    erdos_expression_get_size(item_erdos->expression, item_erdos->renderer, &item->x2, &ascent, &descent);
	    
	    item->x1 = 0;
	    item->y1 = 0;
	    
	    item_erdos->ascent = ascent;
	    item->y2 = ascent + descent;
	  }
}

static void
item_erdos_unrealize (GnomeCanvasItem *item)
{
  ItemErdos *item_erdos;

  item_erdos = ITEM_ERDOS (item);

  if (!item->canvas->aa)
    {
      gdk_gc_unref (item_erdos->gc);
      gdk_gc_unref (item_erdos->bg_gc);
      item_erdos->gc = 0;
      item_erdos->bg_gc = 0;

      erdos_object_unref( (ErdosObject *) item_erdos->renderer );
      erdos_object_unref( (ErdosObject *) item_erdos->style );
    }

  if (parent_class->unrealize)
    (* parent_class->unrealize) (item);
}

/* erdos-renderer-x.c
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

#include "erdos-renderer-x.h"
#include <glib.h>
#include <stdio.h>

static void erdos_renderer_x_init		(ErdosRendererX		 *renderer);
static void erdos_renderer_x_string_extents (ErdosRenderer *renderer, ErdosStyle *style, char *string, double *lbearing, double *rbearing, double *width, double *ascent, double *descent);
static void erdos_renderer_x_render             (ErdosRenderer *renderer, ErdosStyle *style, char *string, double x, double y);
static void erdos_renderer_x_line             (ErdosRenderer *renderer, ErdosStyle *style, double x0, double y0, double x1, double y1, double thickness);

static void erdos_renderer_x_finalize (ErdosObject *object);


static ErdosObjectClass *parent_class = NULL;

static GHashTable *font_hash = NULL;

static guint
erdos_font_hash( gconstpointer v )
{
  ErdosFont *font = (ErdosFont *) v;
  return (font->fontstyle << 16) + g_str_hash( font->fontfamily ) + (int) font->fontsize;
}

static gint
erdos_font_equal( gconstpointer v1,
		  gconstpointer v2 )
{
  ErdosFont *font1 = (ErdosFont *) v1;
  ErdosFont *font2 = (ErdosFont *) v2;
  
  return font1->fontstyle == font2->fontstyle &&
    font1->fontsize == font2->fontsize &&
    strcmp( font1->fontfamily, font2->fontfamily ) == 0;
}

static GdkFont *
get_font( ErdosFont *font )
{
  GdkFont *gdkfont = (GdkFont *) g_hash_table_lookup( font_hash, (gconstpointer) font );
  if ( gdkfont == NULL )
    {
      char fontname[100];
      sprintf(fontname, "-%s-medium-r-normal--%d-0-0-0-p-0-iso8859-1", font->fontfamily, (int) font->fontsize);
      gdkfont = gdk_font_load(fontname);
      erdos_object_ref( (ErdosObject *) font );
      g_hash_table_insert( font_hash, font, gdkfont );
    }
  return gdkfont;
}


ErdosObjectClass *
erdos_renderer_x_get_class (void)
{
  static ErdosObjectClass *renderer_type = NULL;

  if (!renderer_type)
    {
      renderer_type = g_malloc(sizeof(ErdosRendererXClass));
      erdos_renderer_x_class_init( (ErdosRendererXClass *) renderer_type );
    }

  return renderer_type;
}

void
erdos_renderer_x_class_init (ErdosRendererXClass *klass)
{
  ErdosObjectClass *object_class;
  ErdosRendererClass *renderer_class;

  object_class = (ErdosObjectClass*) klass;
  renderer_class = (ErdosRendererClass*) klass;

  parent_class = erdos_renderer_get_class();

  erdos_renderer_class_init (renderer_class);

  renderer_class->string_extents = erdos_renderer_x_string_extents;
  renderer_class->render = erdos_renderer_x_render;
  renderer_class->line = erdos_renderer_x_line;
  object_class->finalize = erdos_renderer_x_finalize;

  font_hash = g_hash_table_new( erdos_font_hash, erdos_font_equal );
}


static void
erdos_renderer_x_init (ErdosRendererX *renderer)
{
  erdos_renderer_init( (ErdosRenderer *) renderer);
  ((ErdosRendererXClass *)(((ErdosObject *)renderer)->klass)) = erdos_renderer_x_get_class();
  renderer->drawable = NULL;
  renderer->gc = NULL;
}

void erdos_renderer_x_string_extents(ErdosRenderer *renderer, ErdosStyle *style, char *string, double *lbearing, double *rbearing, double *width, double *ascent, double *descent)
{
  int n_lbearing;
  int n_rbearing;
  int n_width;
  int n_ascent;
  int n_descent;
  GdkFont *font = get_font( style->font );
  gdk_string_extents( font, string, &n_lbearing, &n_rbearing, &n_width, &n_ascent, &n_descent );
  if ( lbearing )
    *lbearing = n_lbearing;
  if ( rbearing )
    *rbearing = n_rbearing;
  if ( width )
    *width = n_width;
  if ( ascent )
    *ascent = n_ascent;
  if ( descent )
    *descent = n_descent;
}

void
erdos_renderer_x_render(ErdosRenderer *renderer, ErdosStyle *style, char *string, double x, double y)
{
  ErdosRendererX *renderer_x = (ErdosRendererX *) renderer;
  GdkFont *font = get_font( style->font );
  gdk_draw_string( renderer_x->drawable, font, renderer_x->gc, x, y, string );
}

static void
erdos_renderer_x_line             (ErdosRenderer *renderer, ErdosStyle *style, double x0, double y0, double x1, double y1, double thickness)
{
  ErdosRendererX *renderer_x = (ErdosRendererX *) renderer;
  gdk_gc_set_line_attributes( renderer_x->gc, thickness, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER );
  gdk_draw_line( renderer_x->drawable, renderer_x->gc, x0, y0, x1, y1 );
}

ErdosRenderer *
erdos_renderer_x_new(GdkDrawable *drawable, GdkGC *gc)
{
  ErdosRendererX *return_val = g_malloc(sizeof(ErdosRendererX));
  erdos_renderer_x_init(return_val);
  return_val->drawable = drawable;
  return_val->gc = gc;
  /*  gdk_gc_ref( gc ); */
  return (ErdosRenderer *) return_val;
}

static void
erdos_renderer_x_finalize(ErdosObject *object)
{
  /*  ErdosRendererX *renderer = (ErdosRendererX *) object;
      gdk_gc_unref( renderer->gc ); */
}

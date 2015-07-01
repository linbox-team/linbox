/* erdos-style.c
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

#include "erdos-style.h"
#include "erdos-utils.h"

#include <glib.h>

static void erdos_style_init		(ErdosStyle		 *style);
static void erdos_style_real_object_finalize (ErdosObject *object);
static ErdosObject *erdos_style_real_object_clone (ErdosObject *object);
/*static void erdos_real_style_output	(ErdosStyle		 *style); */


static ErdosObjectClass *parent_class = NULL;


ErdosObjectClass *
erdos_style_get_class (void)
{
  static ErdosObjectClass *style_type = NULL;

  if (!style_type)
    {
      style_type = g_malloc(sizeof(ErdosStyleClass));
      erdos_style_class_init( (ErdosStyleClass *) style_type );
    }

  return style_type;
}

void
erdos_style_class_init (ErdosStyleClass *klass)
{
  ErdosObjectClass *object_class;

  object_class = (ErdosObjectClass*) klass;

  parent_class = erdos_object_get_class();

  erdos_object_class_init (object_class);

  object_class->finalize = erdos_style_real_object_finalize;
  object_class->clone = erdos_style_real_object_clone;
}

static void
erdos_style_init (ErdosStyle *style)
{
  erdos_object_init( (ErdosObject *) style );
  ((ErdosStyleClass *)(((ErdosObject *)style)->klass)) = erdos_style_get_class();
  style->font = NULL;
  style->form = ERDOS_INFIX;
  style->thickness = NULL;
  style->script_level = 0;
  style->display_style = 1;
  style->scriptminsize = 8;
  style->scriptsizemultiplier = .71;
}

static void
erdos_style_real_object_finalize (ErdosObject *object)
{
  ErdosStyle *style = ( (ErdosStyle *) object );
  erdos_object_unref( (ErdosObject *) style->font );
}

static ErdosObject *
erdos_style_real_object_clone( ErdosObject *object )
{
  ErdosStyle *style = ( (ErdosStyle *) object );
  ErdosStyle *newstyle = erdos_style_new();
  *newstyle = *style;
  erdos_object_ref( (ErdosObject *) newstyle->font );
  if ( style->thickness )
    newstyle->thickness = erdos_line_thickness_clone(style->thickness);
  return (ErdosObject *) newstyle;
}

ErdosStyle *
erdos_style_new()
{
  ErdosStyle *return_val = g_malloc(sizeof(ErdosStyle));
  erdos_style_init(return_val);
  return (ErdosStyle *) return_val;
}

void
erdos_style_set_font(ErdosStyle *style, ErdosFont *font)
{
  erdos_object_unref( (ErdosObject *) style->font );
  style->font = font;
  erdos_object_ref( (ErdosObject *) style->font );
}

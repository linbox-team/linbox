/* erdos-font.c
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

#include "erdos-font.h"
#include <glib.h>

static void erdos_font_class_init	(ErdosFontClass	 *klass);
static void erdos_font_init	(ErdosFont		 *font);
static void erdos_font_real_object_finalize (ErdosObject *object);
static ErdosObject *erdos_font_real_object_clone (ErdosObject *object);

static ErdosObjectClass *parent_class = NULL;

ErdosObjectClass *
erdos_font_get_class (void)
{
  static ErdosObjectClass *font_type = NULL;

  if (!font_type)
    {
      font_type = g_malloc(sizeof(ErdosFontClass));
      erdos_font_class_init( (ErdosFontClass *) font_type );
    }

  return font_type;
}

static void
erdos_font_class_init (ErdosFontClass *klass)
{
  ErdosObjectClass *object_class;

  object_class = (ErdosObjectClass *) klass;

  parent_class = erdos_object_get_class();

  erdos_object_class_init( (ErdosObjectClass *) klass );

  object_class->finalize = erdos_font_real_object_finalize;
  object_class->clone = erdos_font_real_object_clone;
}

static void
erdos_font_init (ErdosFont *font)
{
  erdos_object_init( (ErdosObject *) font );
  ((ErdosFontClass *)(((ErdosObject *)font)->klass)) = erdos_font_get_class();
}

ErdosFont *
erdos_font_new(void)
{
  ErdosFont *return_val = g_malloc(sizeof(ErdosFont));
  erdos_font_init(return_val);
  return (ErdosFont *) return_val;
}

static void
erdos_font_real_object_finalize (ErdosObject *object)
{
  ErdosFont *font = ( (ErdosFont *) object );
  g_free( font->fontfamily );
}

static ErdosObject *
erdos_font_real_object_clone( ErdosObject *object )
{
  ErdosFont *font = ( (ErdosFont *) object );
  ErdosFont *newfont = erdos_font_new();
  *newfont = *font;
  newfont->fontfamily = g_strdup( font->fontfamily );
  return (ErdosObject *) newfont;
}

/* erdos-renderer.c
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

#include "erdos-renderer.h"

#include <glib.h>

static double erdos_renderer_real_em    (ErdosRenderer *renderer,
					 ErdosStyle *style);
static double erdos_renderer_real_ex    (ErdosRenderer *renderer,
					 ErdosStyle *style);
static double erdos_renderer_real_value_lookup    (ErdosRenderer *renderer,
						   ErdosStyle *style,
						   int which);
static double erdos_renderer_real_unit_convert (ErdosRenderer *renderer,
						ErdosStyle *style,
						double value,
						char unit[3]);
/*static void erdos_real_renderer_output	(ErdosRenderer		 *renderer); */


static ErdosObjectClass *parent_class = NULL;


ErdosObjectClass *
erdos_renderer_get_class (void)
{
  static ErdosObjectClass *renderer_type = NULL;

  if (!renderer_type)
    {
      renderer_type = g_malloc(sizeof(ErdosRendererClass));
      erdos_renderer_class_init( (ErdosRendererClass *) renderer_type );
    }

  return renderer_type;
}

void
erdos_renderer_class_init (ErdosRendererClass *klass)
{
  ErdosObjectClass *object_class;

  object_class = (ErdosObjectClass*) klass;

  parent_class = erdos_object_get_class();

  erdos_object_class_init (object_class);

  klass->string_extents = NULL;
  klass->render = NULL;
  klass->line = NULL;
  klass->em = erdos_renderer_real_em;
  klass->ex = erdos_renderer_real_ex;
  klass->value_lookup = erdos_renderer_real_value_lookup;
  klass->unit_convert = erdos_renderer_real_unit_convert;
}


void
erdos_renderer_init (ErdosRenderer *renderer)
{
  erdos_object_init( (ErdosObject *) renderer );
  ((ErdosRendererClass *)(((ErdosObject *)renderer)->klass)) = erdos_renderer_get_class();
}

void erdos_renderer_string_extents(ErdosRenderer *renderer, ErdosStyle *style, char *string, double *lbearing, double *rbearing, double *width, double *ascent, double *descent)
{
  if ( ERDOS_RENDERER_CLASS(renderer)->string_extents )
    {
      ERDOS_RENDERER_CLASS(renderer)->string_extents( renderer, style, string, lbearing, rbearing, width, ascent, descent );
    }
}

void erdos_renderer_render(ErdosRenderer *renderer, ErdosStyle *style, char *string, double x, double y)
{
  if ( ERDOS_RENDERER_CLASS(renderer)->render )
    {
      ERDOS_RENDERER_CLASS(renderer)->render( renderer, style, string, x, y );
    }
}

void erdos_renderer_line(ErdosRenderer *renderer, ErdosStyle *style, double x0, double y0, double x1, double y1, double thickness)
{
  if ( ERDOS_RENDERER_CLASS(renderer)->line )
    {
      ERDOS_RENDERER_CLASS(renderer)->line( renderer, style, x0, y0, x1, y1, thickness );
    }
}

double erdos_renderer_em(ErdosRenderer *renderer, ErdosStyle *style)
{
  if ( ERDOS_RENDERER_CLASS(renderer)->em )
    {
      return ERDOS_RENDERER_CLASS(renderer)->em( renderer, style );
    }
  else
    return 0.0;
}

double erdos_renderer_ex(ErdosRenderer *renderer, ErdosStyle *style)
{
  if ( ERDOS_RENDERER_CLASS(renderer)->ex )
    {
      return ERDOS_RENDERER_CLASS(renderer)->ex( renderer, style );
    }
  else
    return 0.0;
}

double erdos_renderer_value_lookup(ErdosRenderer *renderer, ErdosStyle *style, int which)
{
  if ( ERDOS_RENDERER_CLASS(renderer)->value_lookup )
    {
      return ERDOS_RENDERER_CLASS(renderer)->value_lookup( renderer, style, which );
    }
  else
    return 0.0;
}

double erdos_renderer_unit_convert(ErdosRenderer *renderer, ErdosStyle *style, double value, char *unit)
{
  if ( ERDOS_RENDERER_CLASS(renderer)->unit_convert )
    {
      return ERDOS_RENDERER_CLASS(renderer)->unit_convert( renderer, style, value, unit );
    }
  else
    return 0.0;
}

double erdos_renderer_real_em(ErdosRenderer *renderer, ErdosStyle *style)
{
  double x;
  erdos_renderer_string_extents(renderer, style, "m", NULL, NULL, &x, NULL, NULL);
  return x;
}

double erdos_renderer_real_ex(ErdosRenderer *renderer, ErdosStyle *style)
{
  double a;
  erdos_renderer_string_extents(renderer, style, "x", NULL, NULL, NULL, &a, NULL);
  return a;
}

double erdos_renderer_real_value_lookup(ErdosRenderer *renderer, ErdosStyle *style, int which)
{
  switch ( which )
    {
      /*
    case 8:
      return erdos_renderer_ex( renderer, style ) * 1.570606452;
    case 9:
      return erdos_renderer_ex( renderer, style ) * .913832258;
    case 10:
      return erdos_renderer_ex( renderer, style ) * 1.02996129;
    case 11:
      return erdos_renderer_ex( renderer, style ) * 1.593832258;
    case 12:
    return erdos_renderer_ex( renderer, style ) * 0.987380645;*/
    case 0x16:
      return erdos_renderer_ex( renderer, style ) * 0.58;
    case 0x88:
      return erdos_renderer_ex( renderer, style ) * 0.0929;
    }
  return 0;
}

double erdos_renderer_real_unit_convert(ErdosRenderer *renderer, ErdosStyle *style, double value, char *unit)
{
  if ( (!unit) || (!*unit) )
    return value;
  switch ( ( (*unit) << 8 ) + *(unit + 1) )
    {
    case ('p' << 8) + 't':
      return value;
    case ('p' << 8) + 'x':
      return value;
    case ('e' << 8) + 'm':
      return erdos_renderer_em(renderer, style) * value;
    case ('e' << 8) + 'x':
      return erdos_renderer_ex(renderer, style) * value;
    case ('i' << 8) + 'n':
      return value * 72.0;
    case ('c' << 8) + 'm':
      return value * (72.0 / 2.54);
    case ('m' << 8) + 'm':
      return value * (72.0 / 25.4);
    case ('p' << 8) + 'c':
      return value / 12.0;
    default:
      return 0;
    }
}

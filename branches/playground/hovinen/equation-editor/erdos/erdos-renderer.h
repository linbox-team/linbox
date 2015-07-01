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
#ifndef __ERDOS_RENDERER_H__
#define __ERDOS_RENDERER_H__

#include "erdos.h"
#include "erdos-style.h"

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

typedef struct _ErdosRenderer       ErdosRenderer;
typedef struct _ErdosRendererClass  ErdosRendererClass;

#define ERDOS_RENDERER_CLASS(object) ((ErdosRendererClass *)((ErdosObject *)(object))->klass)

struct _ErdosRenderer
{
  ErdosObject object;
  /* Put your own, widget-specific fields here */

};

struct _ErdosRendererClass
{
  ErdosObjectClass parent_class;
  /* Signals go here */

  void (*string_extents) (ErdosRenderer *renderer, ErdosStyle *style, char *string, double *lbearing, double *rbearing, double *width, double *ascent, double *descent);
  void (*render) (ErdosRenderer *renderer, ErdosStyle *style, char *string, double x, double y);
  void (*line) (ErdosRenderer *renderer, ErdosStyle *style, double x0, double y0, double x1, double y1, double thickness);
  double (*em) (ErdosRenderer *renderer, ErdosStyle *style);
  double (*ex) (ErdosRenderer *renderer, ErdosStyle *style);
  double (*value_lookup) (ErdosRenderer *renderer, ErdosStyle *style, int which);
  double (*unit_convert) (ErdosRenderer *renderer, ErdosStyle *style, double value, char *unit);
  /*  void (*output)	(ErdosRenderer *renderer); */
};

void erdos_renderer_string_extents(ErdosRenderer *renderer, ErdosStyle *style, char *string, double *lbearing, double *rbearing, double *width, double *ascent, double *descent);
void erdos_renderer_render(ErdosRenderer *renderer, ErdosStyle *style, char *string, double x, double y);
void erdos_renderer_line(ErdosRenderer *renderer, ErdosStyle *style, double x0, double y0, double x1, double y1, double thickness);
double erdos_renderer_em(ErdosRenderer *renderer, ErdosStyle *style);
double erdos_renderer_ex(ErdosRenderer *renderer, ErdosStyle *style);
double erdos_renderer_value_lookup(ErdosRenderer *renderer, ErdosStyle *style, int which);
double erdos_renderer_unit_convert(ErdosRenderer *renderer, ErdosStyle *style, double value, char *unit);

void erdos_renderer_class_init	(ErdosRendererClass	 *klass);
void erdos_renderer_init		(ErdosRenderer		 *renderer);
ErdosObjectClass *erdos_renderer_get_class(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ERDOS_RENDERER_H__ */

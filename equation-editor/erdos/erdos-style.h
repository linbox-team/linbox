/* erdos-style.h
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
#ifndef __ERDOS_STYLE_H__
#define __ERDOS_STYLE_H__

#include "erdos.h"
#include "erdos-types.h"
#include "erdos-font.h"

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

typedef struct _ErdosStyle       ErdosStyle;
typedef struct _ErdosStyleClass  ErdosStyleClass;

#define ERDOS_STYLE_CLASS(object) ((ErdosStyleClass *)((ErdosObject *)(object))->klass)

struct _ErdosStyle
{
  ErdosObject object;
  /* Put your own, widget-specific fields here */
  ErdosFont *font;
  ErdosForm form;
  ErdosLineThickness *thickness;
  int script_level;
  int display_style;
  double scriptminsize;
  double scriptsizemultiplier;
};

struct _ErdosStyleClass
{
  ErdosObjectClass parent_class;
  /* Signals go here */

  /*  void (*output)	(ErdosStyle *style); */
};

void erdos_style_class_init	(ErdosStyleClass	 *klass);
ErdosObjectClass *erdos_style_get_class(void);

void erdos_style_set_font (ErdosStyle *style, ErdosFont *font);
ErdosStyle *erdos_style_new(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ERDOS_STYLE_H__ */

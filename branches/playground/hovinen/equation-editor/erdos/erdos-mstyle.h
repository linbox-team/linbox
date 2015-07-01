/* erdos-mstyle.h
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
#ifndef __ERDOS_MSTYLE_H__
#define __ERDOS_MSTYLE_H__

#include "erdos-presentation.h"

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

#include <glib.h>

typedef struct _ErdosMstyle       ErdosMstyle;
typedef struct _ErdosMstyleClass  ErdosMstyleClass;

struct _ErdosMstyle
{
  ErdosPresentation presentation;
  ErdosExpression *child;
  ErdosLineThickness *thickness;
};

struct _ErdosMstyleClass
{
  ErdosPresentationClass parent_class;
  
  /* Signals go here */
};

ErdosObjectClass *erdos_mstyle_get_class(void);

ErdosExpression *erdos_mstyle_new(ErdosExpression *expression);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ERDOS_MSTYLE_H__ */

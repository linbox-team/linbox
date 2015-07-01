/* erdos-mfrac.h
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
#ifndef __ERDOS_MFRAC_H__
#define __ERDOS_MFRAC_H__

#include "erdos-presentation.h"

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

#include <glib.h>

typedef struct _ErdosMfrac       ErdosMfrac;
typedef struct _ErdosMfracClass  ErdosMfracClass;

struct _ErdosMfrac
{
  ErdosPresentation presentation;
  
  ErdosExpression *numerator;
  ErdosExpression *denominator;

  ErdosLineThickness *thickness;
  double actual_thickness;
};

struct _ErdosMfracClass
{
  ErdosPresentationClass parent_class;
  
  /* Signals go here */
};

ErdosObjectClass *erdos_mfrac_get_class(void);

ErdosExpression *erdos_mfrac_new(ErdosExpression *numerator, ErdosExpression *denominator);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ERDOS_MFRAC_H__ */

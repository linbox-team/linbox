/* erdos-mo.h
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
#ifndef __ERDOS_MO_H__
#define __ERDOS_MO_H__

#include "erdos-presentation.h"

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

typedef struct _ErdosMo       ErdosMo;
typedef struct _ErdosMoClass  ErdosMoClass;

typedef struct _ErdosOperator ErdosOperator;
typedef struct _ErdosOperatorInformation ErdosOperatorInformation;

struct _ErdosMo
{
  ErdosPresentation presentation;
  char *value;
  double lspace;
  double rspace;
};

struct _ErdosMoClass
{
  ErdosPresentationClass parent_class;

  /* Signals go here */
};

struct _ErdosOperator
{
  char *text;
  ErdosForm position;
};

struct _ErdosOperatorInformation
{
  double lspace;
  double rspace;
  int stretchy;
  int fence;
  int separator;
};

ErdosObjectClass *erdos_mo_get_class(void);
void erdos_mo_operator_add(ErdosOperator *operator, ErdosOperatorInformation *information);

ErdosExpression *erdos_mo_new(char *value);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ERDOS_MO_H__ */

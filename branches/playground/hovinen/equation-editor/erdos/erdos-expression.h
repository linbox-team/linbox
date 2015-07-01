/* erdos-expression.h
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
#ifndef __ERDOS_EXPRESSION_H__
#define __ERDOS_EXPRESSION_H__

#include "erdos.h"
#include "erdos-style.h"
#include "erdos-renderer.h"

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

typedef struct _ErdosExpression       ErdosExpression;
typedef struct _ErdosExpressionClass  ErdosExpressionClass;

#define ERDOS_EXPRESSION_CLASS(object) ((ErdosExpressionClass *)((ErdosObject *)(object))->klass)

struct _ErdosExpression
{
  ErdosObject object;
  /* Put your own, widget-specific fields here */
  ErdosStyle *parent_style;
};

struct _ErdosExpressionClass
{
  ErdosObjectClass parent_class;
  /* Signals go here */
  void (*output)	(ErdosExpression *expression);
  void (*draw)	(ErdosExpression *expression, ErdosRenderer *renderer, double x, double y);
  void (*get_size) (ErdosExpression *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent);
  void (*set_style) (ErdosExpression *expression, ErdosStyle *style);
  /* nonzero = failure */
  int (*set_attribute) (ErdosExpression *expression, char *key, char *value);
};

void erdos_expression_output(ErdosExpression *expression);
void erdos_expression_draw(ErdosExpression *expression, ErdosRenderer *renderer, double x, double y);
void erdos_expression_get_size(ErdosExpression *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent);
void erdos_expression_set_style(ErdosExpression *expression, ErdosStyle *style);
int erdos_expression_set_attribute(ErdosExpression *expression, char *key, char *value);

void erdos_expression_class_init	(ErdosExpressionClass	 *klass);
void erdos_expression_init		(ErdosExpression		 *expression);
ErdosObjectClass *erdos_expression_get_class(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ERDOS_EXPRESSION_H__ */

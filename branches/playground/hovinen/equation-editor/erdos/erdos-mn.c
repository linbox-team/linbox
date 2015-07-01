/* erdos-mn.c
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

#include "erdos-mn.h"

#include <glib.h>

static void erdos_mn_init	(ErdosMn		 *mn);
static void erdos_mn_class_init	(ErdosMnClass	 *klass);
static void erdos_mn_real_expression_draw	(ErdosExpression *expression,
						 ErdosRenderer *renderer,
						 double x,
						 double y );
static void erdos_mn_real_expression_get_size	(ErdosExpression		 *expression,
						 ErdosRenderer *renderer,
						 double *x,
						 double *ascent,
						 double *descent);
static void erdos_mn_real_expression_output (ErdosExpression *expression);
static void erdos_mn_real_object_finalize (ErdosObject *object);


static ErdosObjectClass *parent_class = NULL;

ErdosObjectClass *
erdos_mn_get_class (void)
{
  static ErdosObjectClass *mn_type = NULL;

  if (!mn_type)
    {
      mn_type = g_malloc(sizeof(ErdosMnClass));
      erdos_mn_class_init( (ErdosMnClass *) mn_type );
    }

  return mn_type;
}

static void
erdos_mn_class_init (ErdosMnClass *klass)
{
  ErdosObjectClass *object_class = (ErdosObjectClass *) klass;
  ErdosExpressionClass *expression_class = (ErdosExpressionClass *) klass;
  ErdosPresentationClass *presentation_class = (ErdosPresentationClass *) klass;

  parent_class = erdos_presentation_get_class();

  erdos_presentation_class_init( presentation_class );

  object_class->finalize = erdos_mn_real_object_finalize;
  expression_class->output = erdos_mn_real_expression_output;
  expression_class->draw = erdos_mn_real_expression_draw;
  expression_class->get_size = erdos_mn_real_expression_get_size;
}

static void
erdos_mn_init (ErdosMn *mn)
{
  erdos_presentation_init( (ErdosPresentation *) mn );
  ((ErdosMnClass *)(((ErdosObject *)mn)->klass)) = erdos_mn_get_class();
  mn->value = NULL;
}

static void
erdos_mn_real_expression_draw (ErdosExpression *expression, ErdosRenderer *renderer, double x, double y)
{
  ErdosMn *mn = (ErdosMn *) expression;
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  erdos_renderer_render( renderer, ((ErdosExpression *) mn)->parent_style, mn->value, x, y);
}

static void
erdos_mn_real_expression_get_size (ErdosExpression *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent )
{
  ErdosMn *mn = (ErdosMn *) expression;
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  erdos_renderer_string_extents( renderer, ((ErdosExpression *) mn)->parent_style, mn->value, NULL, NULL, width, ascent, descent );
}

static void
erdos_mn_real_expression_output (ErdosExpression *expression)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
}

void
erdos_mn_real_object_finalize( ErdosObject *object )
{
  ErdosMn *mn = (ErdosMn *) object;
  g_free( mn->value );
  if ( parent_class->finalize )
    {
      parent_class->finalize( object );
    }
}

ErdosExpression *
erdos_mn_new(char *value)
{
  ErdosMn *return_val = g_malloc(sizeof(ErdosMn));
  erdos_mn_init(return_val);
  return_val->value = g_strdup(value);
  return (ErdosExpression *) return_val;
}

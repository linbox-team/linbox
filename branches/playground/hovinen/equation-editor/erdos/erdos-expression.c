/* erdos-expression.c
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

#include "erdos-expression.h"
#include <glib.h>

static void erdos_real_expression_output	(ErdosExpression		 *expression);
static void erdos_real_expression_draw	(ErdosExpression		 *expression,
					 ErdosRenderer                   *renderer,
					 double x,
					 double y);
static void erdos_real_expression_get_size	(ErdosExpression		 *expression,
						 ErdosRenderer                   *renderer,
						 double *width,
						 double *ascent,
						 double *descent);
static void erdos_real_expression_set_style (ErdosExpression *expression,
					     ErdosStyle *style);
static void erdos_real_expression_finalize (ErdosObject *object);


static ErdosObjectClass *parent_class = NULL;


ErdosObjectClass *
erdos_expression_get_class (void)
{
  static ErdosObjectClass *expression_type = NULL;

  if (!expression_type)
    {
      expression_type = g_malloc(sizeof(ErdosExpressionClass));
      erdos_expression_class_init( (ErdosExpressionClass *) expression_type );
    }

  return expression_type;
}

void
erdos_expression_class_init (ErdosExpressionClass *klass)
{
  ErdosObjectClass *object_class;

  object_class = (ErdosObjectClass*) klass;

  parent_class = erdos_object_get_class();

  erdos_object_class_init (object_class);

  object_class->finalize = erdos_real_expression_finalize;
  klass->output = erdos_real_expression_output;
  klass->draw = erdos_real_expression_draw;
  klass->get_size = erdos_real_expression_get_size;
  klass->set_style = erdos_real_expression_set_style;
  klass->set_attribute = NULL;
}


void
erdos_expression_init (ErdosExpression *expression)
{
  erdos_object_init( (ErdosObject *) expression );
  ERDOS_EXPRESSION_CLASS(expression) = erdos_expression_get_class();
  expression->parent_style = NULL;
}


static void
erdos_real_expression_output (ErdosExpression *expression)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
}

static void
erdos_real_expression_draw (ErdosExpression *expression, ErdosRenderer *renderer, double x, double y)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
}

static void
erdos_real_expression_get_size (ErdosExpression *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
}

static void
erdos_real_expression_set_style (ErdosExpression *expression, ErdosStyle *style)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  if ( expression->parent_style )
    {
      erdos_object_unref( (ErdosObject *) expression->parent_style );
    }
  expression->parent_style = style;
  erdos_object_ref( (ErdosObject *) style );
}

static void
erdos_real_expression_finalize (ErdosObject *object)
{
  erdos_object_unref( (ErdosObject *) ((ErdosExpression *)object)->parent_style );
}

void
erdos_expression_output(ErdosExpression *expression)
{
  if ( ERDOS_EXPRESSION_CLASS(expression)->output )
    {
      ERDOS_EXPRESSION_CLASS(expression)->output( expression );
    }
}

void
erdos_expression_draw(ErdosExpression *expression, ErdosRenderer *renderer, double x, double y)
{
  if ( ERDOS_EXPRESSION_CLASS(expression)->draw )
    {
      ERDOS_EXPRESSION_CLASS(expression)->draw( expression, renderer, x, y );
    }
}

void
erdos_expression_get_size(ErdosExpression *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent)
{
  if ( ERDOS_EXPRESSION_CLASS(expression)->get_size )
    {
      ERDOS_EXPRESSION_CLASS(expression)->get_size( expression, renderer, width, ascent, descent );
    }
}

void
erdos_expression_set_style(ErdosExpression *expression, ErdosStyle *style)
{
  if ( ERDOS_EXPRESSION_CLASS(expression)->set_style )
    {
      ERDOS_EXPRESSION_CLASS(expression)->set_style( expression, style );
    }
}

int
erdos_expression_set_attribute(ErdosExpression *expression, char *key, char *value)
{
  if ( ERDOS_EXPRESSION_CLASS(expression)->set_attribute )
    {
      return ERDOS_EXPRESSION_CLASS(expression)->set_attribute( expression, key, value );
    }
  else
    return 0;
}

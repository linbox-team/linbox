/* erdos-mstyle.c
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

#include "erdos-mstyle.h"
#include "erdos-utils.h"

static void erdos_mstyle_init	(ErdosMstyle		 *mstyle);
static void erdos_mstyle_class_init	(ErdosMstyleClass	 *klass);
static void erdos_mstyle_real_object_finalize (ErdosObject *object);
static void erdos_mstyle_real_expression_draw	(ErdosExpression		 *expression, ErdosRenderer *renderer, double x, double y);
static void erdos_mstyle_real_expression_get_size	(ErdosExpression		 *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent);
static void erdos_mstyle_real_expression_output (ErdosExpression *expression);
static void erdos_mstyle_real_expression_set_style (ErdosExpression *expression, ErdosStyle *style);
static int erdos_mstyle_real_expression_set_attribute (ErdosExpression *expression, char *key, char *value);

static ErdosObjectClass *parent_class = NULL;

ErdosObjectClass *
erdos_mstyle_get_class (void)
{
  static ErdosObjectClass *mstyle_type = NULL;

  if (!mstyle_type)
    {
      mstyle_type = g_malloc(sizeof(ErdosMstyleClass));
      erdos_mstyle_class_init( (ErdosMstyleClass *) mstyle_type );
    }

  return mstyle_type;
}

static void
erdos_mstyle_class_init (ErdosMstyleClass *klass)
{
  ErdosObjectClass *object_class;
  ErdosExpressionClass *expression_class;
  ErdosPresentationClass *presentation_class;

  object_class = (ErdosObjectClass *) klass;
  expression_class = (ErdosExpressionClass *) klass;
  presentation_class = (ErdosPresentationClass *) klass;

  parent_class = erdos_presentation_get_class();

  erdos_presentation_class_init( (ErdosPresentationClass *) klass );

  object_class->finalize = erdos_mstyle_real_object_finalize;
  expression_class->output = erdos_mstyle_real_expression_output;
  expression_class->draw = erdos_mstyle_real_expression_draw;
  expression_class->get_size = erdos_mstyle_real_expression_get_size;
  expression_class->set_style = erdos_mstyle_real_expression_set_style;
  expression_class->set_attribute = erdos_mstyle_real_expression_set_attribute;
}

static void
erdos_mstyle_init (ErdosMstyle *mstyle)
{
  erdos_presentation_init( (ErdosPresentation *) mstyle );
  ((ErdosMstyleClass *)(((ErdosObject *)mstyle)->klass)) = erdos_mstyle_get_class();
  mstyle->child = NULL;
}

ErdosExpression *
erdos_mstyle_new(ErdosExpression *child)
{
  ErdosMstyle *return_val = g_malloc(sizeof(ErdosMstyle));
  erdos_mstyle_init(return_val);
  return_val->child = child;
  erdos_object_ref( (ErdosObject *) child );
  return (ErdosExpression *) return_val;
}

static void
erdos_mstyle_real_expression_draw (ErdosExpression *expression, ErdosRenderer *renderer, double x, double y)
{
  ErdosMstyle *mstyle = ((ErdosMstyle *)expression);
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */

  erdos_expression_draw( mstyle->child, renderer, x, y );
}

static void
erdos_mstyle_real_expression_get_size (ErdosExpression *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent)
{
  ErdosMstyle *mstyle = ((ErdosMstyle *)expression);
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  erdos_expression_get_size(mstyle->child, renderer, width, ascent, descent);
}

static void
erdos_mstyle_real_expression_output (ErdosExpression *expression)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
}

static void
erdos_mstyle_real_expression_set_style (ErdosExpression *expression, ErdosStyle *style)
{
  ErdosMstyle *mstyle = ((ErdosMstyle *)expression);
  g_return_if_fail (expression != NULL);
  if ( ((ErdosExpressionClass *)parent_class)->set_style )
    {
      ((ErdosExpressionClass *)parent_class)->set_style( expression, style );
    }
  
  erdos_object_ref( (ErdosObject *) style );
  
  if ( mstyle->thickness )
    {
      style = (ErdosStyle *) erdos_object_cow( (ErdosObject *) style );
      g_free( style->thickness );
      style->thickness = erdos_line_thickness_clone( mstyle->thickness );
    }
  
  erdos_expression_set_style( mstyle->child, style );
  erdos_object_unref( (ErdosObject *)style );
}

static void erdos_mstyle_real_object_finalize (ErdosObject *object)
{
  ErdosMstyle *mstyle = ( (ErdosMstyle *) object );
  g_return_if_fail (object != NULL);
  
  erdos_object_unref( (ErdosObject *) mstyle->child );
  
  if ( parent_class->finalize )
    {
      parent_class->finalize( object );
    }
}


static int
erdos_mstyle_real_expression_set_attribute(ErdosExpression *expression, char *key, char *value)
{
  ErdosMstyle *mstyle = ((ErdosMstyle *)expression);
  /*  g_return_if_fail_val (expression != NULL, -1); */
  if ( ! strcmp( key, "linethickness" ) )
    {
      mstyle->thickness = erdos_line_thickness_new( value );
      if ( mstyle->thickness == NULL )
	return -1;
      else
	return 0;
    }
  else if ( ((ErdosExpressionClass *)parent_class)->set_attribute )
    return ((ErdosExpressionClass *)parent_class)->set_attribute( expression, key, value );
  else
    return 0;
}
 

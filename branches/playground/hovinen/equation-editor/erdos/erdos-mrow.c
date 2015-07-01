/* erdos-mrow.c
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

#include "erdos-mrow.h"

static void erdos_mrow_init	(ErdosMrow		 *mrow);
static void erdos_mrow_class_init	(ErdosMrowClass	 *klass);
static void erdos_mrow_real_object_finalize (ErdosObject *object);
static void erdos_mrow_real_expression_draw	(ErdosExpression		 *expression, ErdosRenderer *renderer, double x, double y);
static void erdos_mrow_real_expression_get_size	(ErdosExpression		 *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent);
static void erdos_mrow_real_expression_output (ErdosExpression *expression);
static void erdos_mrow_real_expression_set_style (ErdosExpression *expression, ErdosStyle *style);

static ErdosObjectClass *parent_class = NULL;

ErdosObjectClass *
erdos_mrow_get_class (void)
{
  static ErdosObjectClass *mrow_type = NULL;

  if (!mrow_type)
    {
      mrow_type = g_malloc(sizeof(ErdosMrowClass));
      erdos_mrow_class_init( (ErdosMrowClass *) mrow_type );
    }

  return mrow_type;
}

static void
erdos_mrow_class_init (ErdosMrowClass *klass)
{
  ErdosObjectClass *object_class;
  ErdosExpressionClass *expression_class;
  ErdosPresentationClass *presentation_class;

  object_class = (ErdosObjectClass *) klass;
  expression_class = (ErdosExpressionClass *) klass;
  presentation_class = (ErdosPresentationClass *) klass;

  parent_class = erdos_presentation_get_class();

  erdos_presentation_class_init( (ErdosPresentationClass *) klass );

  object_class->finalize = erdos_mrow_real_object_finalize;
  expression_class->output = erdos_mrow_real_expression_output;
  expression_class->draw = erdos_mrow_real_expression_draw;
  expression_class->get_size = erdos_mrow_real_expression_get_size;
  expression_class->set_style = erdos_mrow_real_expression_set_style;
}

static void
erdos_mrow_init (ErdosMrow *mrow)
{
  erdos_presentation_init( (ErdosPresentation *) mrow );
  ((ErdosMrowClass *)(((ErdosObject *)mrow)->klass)) = erdos_mrow_get_class();
  mrow->children = NULL;
}

ErdosExpression *
erdos_mrow_new(void)
{
  ErdosMrow *return_val = g_malloc(sizeof(ErdosMrow));
  erdos_mrow_init(return_val);
  return (ErdosExpression *) return_val;
}

void
erdos_mrow_append(ErdosMrow *mrow, ErdosExpression *expression)
{
  mrow->children = g_list_append(mrow->children, expression);
  erdos_object_ref( (ErdosObject *) expression );
  if ( ((ErdosExpression *) mrow)->parent_style )
    {
      erdos_expression_set_style( expression, ((ErdosExpression *) mrow)->parent_style );
    }
}

static void
erdos_mrow_real_expression_draw (ErdosExpression *expression, ErdosRenderer *renderer, double x, double y)
{
  ErdosMrow *mrow = ((ErdosMrow *)expression);
  GList *children;
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
  for ( children = mrow->children; children; children = g_list_next(children) )
    {
      double thisx;
      erdos_expression_draw((ErdosExpression *)children->data, renderer, x, y);
      erdos_expression_get_size((ErdosExpression *)children->data, renderer, &thisx, NULL, NULL);
      x += thisx;
    }
}

static void
erdos_mrow_real_expression_get_size (ErdosExpression *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent)
{
  ErdosMrow *mrow = ((ErdosMrow *)expression);
  GList *children;
  double thisx, thisa, thisd;
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  if ( ascent )
    *ascent = 0;
  if ( descent )
    *descent = 0;
  if ( width )
    *width = 0;
  if ( width || ascent || descent )
    {
      for ( children = mrow->children; children; children = g_list_next(children) )
	{
	  erdos_expression_get_size((ErdosExpression *)children->data, renderer, &thisx, &thisa, &thisd);
	  if ( ascent && thisa > *ascent )
	    *ascent = thisa;
	  if ( descent && thisd > *descent )
	    *descent = thisd;
	  if ( width )
	    *width += thisx;
	}
    }
}

static void
erdos_mrow_real_expression_output (ErdosExpression *expression)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
}

static void
erdos_mrow_real_expression_set_style (ErdosExpression *expression, ErdosStyle *stylep)
{
  ErdosMrow *mrow = ((ErdosMrow *)expression);
  GList *children;
  ErdosForm form = ERDOS_PREFIX;
  ErdosStyle *style = stylep;
  g_return_if_fail (expression != NULL);
  if ( ((ErdosExpressionClass *)parent_class)->set_style )
    {
      ((ErdosExpressionClass *)parent_class)->set_style( expression, stylep );
    }
  for ( children = mrow->children; children; children = g_list_next(children) )
    {
      if ( (  children->next  &&   children->prev ) ||
	   ((!children->next) && (!children->prev)) )
	form = ERDOS_INFIX;
      else if ( children->prev )
	form = ERDOS_SUFFIX;
      if ( form != style->form )
	{
	  erdos_object_ref( (ErdosObject *) style );
	  style = (ErdosStyle *) erdos_object_cow( (ErdosObject *) style );
	  style->form = form;
	  erdos_expression_set_style((ErdosExpression *)children->data, style);
	  erdos_object_unref( (ErdosObject *) style );
	}
      else
	erdos_expression_set_style((ErdosExpression *)children->data, style);
    }
}

static void erdos_mrow_real_object_finalize (ErdosObject *object)
{
  ErdosMrow *mrow = ( (ErdosMrow *) object );
  GList *children;
  g_return_if_fail (object != NULL);
  for ( children = mrow->children; children; children = g_list_next(children) )
    {
      erdos_object_unref( (ErdosObject *) children->data );
    }
  g_list_free( mrow->children );
  if ( parent_class->finalize )
    {
      parent_class->finalize( object );
    }
}

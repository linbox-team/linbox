/* erdos-mo.c
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

#include "erdos-mo.h"

#include <glib.h>

static void erdos_mo_init	(ErdosMo		 *mo);
static void erdos_mo_class_init	(ErdosMoClass	 *klass);
static void erdos_mo_real_expression_draw	(ErdosExpression *expression,
						 ErdosRenderer *renderer,
						 double x,
						 double y );
static void erdos_mo_real_expression_get_size	(ErdosExpression		 *expression,
						 ErdosRenderer *renderer,
						 double *x,
						 double *ascent,
						 double *descent);
static void erdos_mo_real_expression_output (ErdosExpression *expression);
static void erdos_mo_real_object_finalize (ErdosObject *object);
static void erdos_mo_real_expression_set_style (ErdosExpression *expression,
						ErdosStyle *style);

static ErdosObjectClass *parent_class = NULL;

static GHashTable *operator_table = NULL;

static guint
op_hash( gconstpointer v )
{
  ErdosOperator *operator = (ErdosOperator *) v;
  guint ret_val = g_str_hash( operator->text );
  if ( operator->position == ERDOS_PREFIX )
    ret_val &= 0xFFFF0000;
  else if ( operator->position == ERDOS_SUFFIX )
    ret_val &= 0x0000FFFF;
  return ret_val;
}

static gint
op_equal( gconstpointer v1, gconstpointer v2 )
{  
  ErdosOperator *operator1 = (ErdosOperator *) v1;
  ErdosOperator *operator2 = (ErdosOperator *) v2;
  return strcmp( operator1->text, operator2->text ) == 0 &&
    operator1->position == operator2->position;
}

void
erdos_mo_operator_add(ErdosOperator *operator, ErdosOperatorInformation *information)
{
  ErdosOperator *newop = g_new( ErdosOperator, 1 );
  ErdosOperatorInformation *newinfo = g_new( ErdosOperatorInformation, 1 );
  
  *newop = *operator;
  newop->text = g_strdup(operator->text);
  *newinfo = *information;
  
  erdos_mo_get_class();
  g_hash_table_insert(operator_table, newop, newinfo);
}

ErdosObjectClass *
erdos_mo_get_class (void)
{
  static ErdosObjectClass *mo_type = NULL;

  if (!mo_type)
    {
      mo_type = g_malloc(sizeof(ErdosMoClass));
      erdos_mo_class_init( (ErdosMoClass *) mo_type );
    }

  return mo_type;
}

static void
erdos_mo_class_init (ErdosMoClass *klass)
{
  ErdosObjectClass *object_class = (ErdosObjectClass *) klass;
  ErdosExpressionClass *expression_class = (ErdosExpressionClass *) klass;
  ErdosPresentationClass *presentation_class = (ErdosPresentationClass *) klass;
  ErdosOperator op;
  ErdosOperatorInformation info;

  parent_class = erdos_presentation_get_class();

  erdos_presentation_class_init( presentation_class );

  object_class->finalize = erdos_mo_real_object_finalize;
  expression_class->output = erdos_mo_real_expression_output;
  expression_class->draw = erdos_mo_real_expression_draw;
  expression_class->get_size = erdos_mo_real_expression_get_size;
  expression_class->set_style = erdos_mo_real_expression_set_style;

  operator_table = g_hash_table_new(op_hash, op_equal);
  
  info.fence = 0;
  info.stretchy = 0;
  info.separator = 0;
  
  op.text = "+";
  op.position = ERDOS_INFIX;
  info.lspace = 0.22222;
  info.rspace = 0.22222;
  erdos_mo_operator_add(&op, &info);
  op.text = "-";
  erdos_mo_operator_add(&op, &info);
  op.text = "=";
  info.lspace = 0.27777;
  info.rspace = 0.27777;
  erdos_mo_operator_add(&op, &info);

  op.text = "+";
  op.position = ERDOS_PREFIX;
  info.lspace = 0.0;
  info.rspace = 0.05555;
  erdos_mo_operator_add(&op, &info);
  op.text = "-";
  erdos_mo_operator_add(&op, &info);

  op.text = "(";
  op.position = ERDOS_PREFIX;
  info.lspace = 0.0;
  info.rspace = 0.0;
  info.fence = 1;
  info.stretchy = 1;
  erdos_mo_operator_add(&op, &info);
  op.text = ")";
  op.position = ERDOS_SUFFIX;
  erdos_mo_operator_add(&op, &info);
}

static void
erdos_mo_init (ErdosMo *mo)
{
  erdos_presentation_init( (ErdosPresentation *) mo );
  ((ErdosMoClass *)(((ErdosObject *)mo)->klass)) = erdos_mo_get_class();
  mo->value = NULL;
  mo->lspace = 0.2777;
  mo->rspace = 0.2777;
}


static void
erdos_mo_real_expression_draw (ErdosExpression *expression, ErdosRenderer *renderer, double x, double y)
{
  ErdosMo *mo = (ErdosMo *) expression;
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  x += erdos_renderer_em( renderer, expression->parent_style ) * ( mo->lspace );
  erdos_renderer_render( renderer, expression->parent_style, mo->value, x, y);
}

static void
erdos_mo_real_expression_get_size (ErdosExpression *expression, ErdosRenderer *renderer, double *width, double *ascent, double *descent )
{
  ErdosMo *mo = (ErdosMo *) expression;
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  erdos_renderer_string_extents( renderer, expression->parent_style, mo->value, NULL, NULL, width, ascent, descent );
  if ( width )
    *width += erdos_renderer_em( renderer, expression->parent_style ) * ( mo->lspace + mo->rspace );
}

static void
erdos_mo_real_expression_output (ErdosExpression *expression)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
}

static void
erdos_mo_real_expression_set_style (ErdosExpression *expression, ErdosStyle *style)
{
  ErdosMo *mo = (ErdosMo *) expression;
  g_return_if_fail (expression != NULL);
  if ( ((ErdosExpressionClass *) parent_class)->set_style )
    {
      ((ErdosExpressionClass *) parent_class)->set_style(expression, style);
    }
  if ( expression->parent_style )
    {
      ErdosOperator op = { mo->value, expression->parent_style->form };
      ErdosOperatorInformation *info;
      info = g_hash_table_lookup( operator_table, &op );
      if ( info == NULL )
	{
	  op.position = ERDOS_INFIX;
	  info = g_hash_table_lookup( operator_table, &op );
	}
      if ( info == NULL )
	{
	  op.position = ERDOS_SUFFIX;
	  info = g_hash_table_lookup( operator_table, &op );
	}
      if ( info == NULL )
	{
	  op.position = ERDOS_PREFIX;
	  info = g_hash_table_lookup( operator_table, &op );
	}
      if ( info == NULL )
	{
	  mo->lspace = 0.27777;
	  mo->rspace = 0.27777;
	}
      else
	{
	  mo->lspace = info->lspace;
	  mo->rspace = info->rspace;
	}
    }
}

void
erdos_mo_real_object_finalize( ErdosObject *object )
{
  ErdosMo *mo = (ErdosMo *) object;
  g_free( mo->value );
  if ( parent_class->finalize )
    {
      parent_class->finalize( object );
    }
}

ErdosExpression *
erdos_mo_new(char *value)
{
  ErdosMo *return_val = g_malloc(sizeof(ErdosMo));
  erdos_mo_init(return_val);
  return_val->value = g_strdup(value);
  return (ErdosExpression *) return_val;
}

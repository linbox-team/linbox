/* erdos-mfrac.c
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

#include "erdos-mfrac.h"
#include "erdos-utils.h"

static void erdos_mfrac_init	(ErdosMfrac		 *mfrac);
static void erdos_mfrac_class_init	(ErdosMfracClass	 *klass);
static void erdos_mfrac_real_expression_draw	(ErdosExpression		 *expression, ErdosRenderer *renderer, double x, double y);
static void erdos_mfrac_real_expression_get_size	(ErdosExpression		 *expression, ErdosRenderer *renderer, double *x, double *ascent, double *descent);
static void erdos_mfrac_real_expression_set_style (ErdosExpression *expression, ErdosStyle *style);
static void erdos_mfrac_real_expression_output (ErdosExpression *expression);
static void erdos_mfrac_real_object_finalize (ErdosObject *object);
static int erdos_mfrac_real_expression_set_attribute (ErdosExpression *expression, char *key, char *value);

static ErdosObjectClass *parent_class = NULL;

#define FONT_HEIGHT (((expression)->parent_style->display_style) ? 3.5 : 1.5)

ErdosObjectClass *
erdos_mfrac_get_class (void)
{
  static ErdosObjectClass *mfrac_type = NULL;

  if (!mfrac_type)
    {
      mfrac_type = g_malloc(sizeof(ErdosMfracClass));
      erdos_mfrac_class_init( (ErdosMfracClass *) mfrac_type );
    }

  return mfrac_type;
}

static void
erdos_mfrac_class_init (ErdosMfracClass *klass)
{
  ErdosObjectClass *object_class;
  ErdosExpressionClass *expression_class;
  ErdosPresentationClass *presentation_class;

  object_class = (ErdosObjectClass *) klass;
  expression_class = (ErdosExpressionClass *) klass;
  presentation_class = (ErdosPresentationClass *) klass;

  parent_class = erdos_presentation_get_class();

  erdos_presentation_class_init( (ErdosPresentationClass *) klass );

  object_class->finalize = erdos_mfrac_real_object_finalize;
  expression_class->output = erdos_mfrac_real_expression_output;
  expression_class->draw = erdos_mfrac_real_expression_draw;
  expression_class->get_size = erdos_mfrac_real_expression_get_size;
  expression_class->set_style = erdos_mfrac_real_expression_set_style;
  expression_class->set_attribute = erdos_mfrac_real_expression_set_attribute;
}

static void
erdos_mfrac_init (ErdosMfrac *mfrac)
{
  erdos_presentation_init( (ErdosPresentation *) mfrac );
  ((ErdosMfracClass *)(((ErdosObject *)mfrac)->klass)) = erdos_mfrac_get_class();
  mfrac->numerator = NULL;
  mfrac->denominator = NULL;
  mfrac->thickness = NULL;
  mfrac->actual_thickness = 0.0;
}

ErdosExpression *
erdos_mfrac_new(ErdosExpression *numerator, ErdosExpression *denominator)
{
  ErdosMfrac *return_val = g_malloc(sizeof(ErdosMfrac));
  erdos_mfrac_init(return_val);
  return_val->numerator = numerator;
  return_val->denominator = denominator;
  return (ErdosExpression *) return_val;
}

static void
erdos_mfrac_real_expression_draw (ErdosExpression *expression, ErdosRenderer *renderer, double x, double y)
{
  ErdosMfrac *mfrac = ((ErdosMfrac *)expression);
  double numx, denx, numd, dena;
  double rule_thick;
  double u, v;
  int uref, vref;
  double phi, rho;
  double width;
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
  erdos_expression_get_size(mfrac->numerator, renderer, &numx, NULL, &numd);
  width = numx;
  erdos_expression_get_size(mfrac->denominator, renderer, &denx, &dena, NULL);      
  if ( width < denx )
    width = denx;
  
  if ( mfrac->thickness )
    {
      rule_thick = erdos_line_thickness_calculate_actual(mfrac->thickness, renderer, expression->parent_style);
    }
  else
    {
      rule_thick = erdos_line_thickness_calculate_actual(expression->parent_style->thickness, renderer, expression->parent_style);
    }	  
	  
  if ( (expression)->parent_style->display_style )
    {
      uref = 8;
      vref = 11;
    }
  else
    {
      if ( rule_thick == 0 )
	uref = 10;
      else
	uref = 9;
      vref = 12;
    }
  u = erdos_renderer_value_lookup( renderer, expression->parent_style, uref );
  v = erdos_renderer_value_lookup( renderer, expression->parent_style, vref );
	  
  if ( rule_thick == 0 )
    {
      if ( (expression)->parent_style->display_style )
	phi = 7 * erdos_renderer_value_lookup( renderer, expression->parent_style, 0x88 );
      else
	phi = 3 * erdos_renderer_value_lookup( renderer, expression->parent_style, 0x88 );
      rho = (u - numd) - (dena - v);
      if ( rho < phi )
	{
	  u += (phi - rho) / 2.0;
	  v += (phi - rho) / 2.0;
	}
    }
  else
    {
      double a = erdos_renderer_value_lookup( renderer, expression->parent_style, 22 );
      if ( (expression)->parent_style->display_style )
	phi = 3 * rule_thick;
      else
	phi = rule_thick;
      if ( (u - numd) - (a + phi / 2.0) < phi )
	u += phi - ((u - numd) - (a + phi / 2.0));
      if ( (a - phi / 2.0) - (dena - v) < phi )
	v += phi - ((a - phi / 2.0) - (dena - v));
      erdos_renderer_line(renderer, expression->parent_style, x, y - a, x + width, y - a, rule_thick );/*  - 0.5 * rule_thick, x + width, y - a - 0.5 * rule_thick, rule_thick );*/
    }
  erdos_expression_draw(mfrac->numerator, renderer, x + ( width - numx ) / 2.0, y - u );
  erdos_expression_draw(mfrac->denominator, renderer, x + ( width - denx ) / 2.0, y + v );
}

static void
erdos_mfrac_real_expression_get_size (ErdosExpression *expression, ErdosRenderer *renderer, double *x, double *ascent, double *descent)
{
  ErdosMfrac *mfrac = ((ErdosMfrac *)expression);
  double thisx, numa, numd, dena, dend;
  double rule_thick;
  double u, v;
  double phi, rho;
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
  if ( x )
    *x = 0;
  if ( x || ascent || descent )
    {
      erdos_expression_get_size(mfrac->numerator, renderer, &thisx, &numa, &numd);
      if ( x && thisx > *x )
	*x = thisx;
      erdos_expression_get_size(mfrac->denominator, renderer, &thisx, &dena, &dend);      
      if ( x && thisx > *x )
	*x = thisx;
      
      if ( ascent || descent )
	{
	  int uref, vref;
	  if ( mfrac->thickness )
	    {
	      rule_thick = erdos_line_thickness_calculate_actual(mfrac->thickness, renderer, expression->parent_style);
	    }
	  else
	    {
	      rule_thick = erdos_line_thickness_calculate_actual(expression->parent_style->thickness, renderer, expression->parent_style);
	    }	  
	  
	  if ( (expression)->parent_style->display_style )
	    {
	      uref = 8;
	      vref = 11;
	    }
	  else
	    {
	      if ( rule_thick == 0 )
		uref = 10;
	      else
		uref = 9;
	      vref = 12;
	    }
	  u = erdos_renderer_value_lookup( renderer, expression->parent_style, uref );
	  v = erdos_renderer_value_lookup( renderer, expression->parent_style, vref );
	  
	  if ( rule_thick == 0 )
	    {
	      if ( (expression)->parent_style->display_style )
		phi = 7 * erdos_renderer_value_lookup( renderer, expression->parent_style, 0x88 );
	      else
		phi = 3 * erdos_renderer_value_lookup( renderer, expression->parent_style, 0x88 );
	      rho = (u - numd) - (dena - v);
	      if ( rho < phi )
		{
		  u += (phi - rho) / 2.0;
		  v += (phi - rho) / 2.0;
		}
	    }
	  else
	    {
	      double a = erdos_renderer_value_lookup( renderer, expression->parent_style, 22 );
	      if ( (expression)->parent_style->display_style )
		phi = 3 * rule_thick;
	      else
		phi = rule_thick;
	      if ( (u - numd) - (a + phi / 2.0) < phi )
		u += phi - ((u - numd) - (a + phi / 2.0));
	      if ( (a - phi / 2.0) - (dena - v) < phi )
		v += phi - ((a - phi / 2.0) - (dena - v));
	    }
	  if ( ascent )
	    *ascent = u + numa;
	  if ( descent )
	    *descent = v + dend;
	}
    }
}

static void
erdos_mfrac_real_expression_output (ErdosExpression *expression)
{
  g_return_if_fail (expression != NULL);
  /*  g_return_if_fail (ERDOS_IS_EXPRESSION (expression)); */
  
}

static void
erdos_mfrac_real_expression_set_style(ErdosExpression *expression, ErdosStyle *style)
{
  ErdosMfrac *mfrac = ((ErdosMfrac *)expression);
  g_return_if_fail (expression != NULL);
  if ( ((ErdosExpressionClass *)parent_class)->set_style )
    {
      ((ErdosExpressionClass *)parent_class)->set_style(expression, style);
    }
  erdos_object_ref( (ErdosObject *) style );
  style = (ErdosStyle *) erdos_object_cow( (ErdosObject *) style );
  if ( style->display_style )
    style->display_style = 0;
  else
    {
      style->script_level ++;
      if ( style->font->fontsize > style->scriptminsize )
	{
	  style->font = (ErdosFont *) erdos_object_cow( (ErdosObject *) style->font );
	  style->font->fontsize = style->font->fontsize * style->scriptsizemultiplier;
	  if ( style->font->fontsize < style->scriptminsize )
	    style->font->fontsize = style->scriptminsize;
	}
    }
  
  erdos_expression_set_style(mfrac->numerator, style);
  erdos_expression_set_style(mfrac->denominator, style);
  erdos_object_unref( (ErdosObject *) style );
}

static void
erdos_mfrac_real_object_finalize(ErdosObject *object)
{
  ErdosMfrac *mfrac = ((ErdosMfrac *)object);
  g_return_if_fail (object != NULL);
  erdos_object_unref((ErdosObject *)mfrac->numerator);
  erdos_object_unref((ErdosObject *)mfrac->denominator);
}

static int
erdos_mfrac_real_expression_set_attribute(ErdosExpression *expression, char *key, char *value)
{
  ErdosMfrac *mfrac = ((ErdosMfrac *)expression);
  /*  g_return_if_fail_val (expression != NULL, -1); */
  if ( ! strcmp( key, "linethickness" ) )
    {
      mfrac->thickness = erdos_line_thickness_new( value );
      if ( mfrac->thickness == NULL )
	return -1;
      else
	return 0;
    }
  else if ( ((ErdosExpressionClass *)parent_class)->set_attribute )
    return ((ErdosExpressionClass *)parent_class)->set_attribute( expression, key, value );
  else
    return 0;
}

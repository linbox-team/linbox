/* erdos-presentation.c
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

#include "erdos-presentation.h"

#include <glib.h>

static ErdosObjectClass *parent_class = NULL;

ErdosObjectClass *
erdos_presentation_get_class (void)
{
  static ErdosObjectClass *presentation_type = NULL;

  if (!presentation_type)
    {
      presentation_type = g_malloc(sizeof(ErdosPresentationClass));
      erdos_presentation_class_init( (ErdosPresentationClass *) presentation_type );
    }

  return presentation_type;
}

void
erdos_presentation_class_init (ErdosPresentationClass *klass)
{
  ErdosExpressionClass *expression_class;

  expression_class = (ErdosExpressionClass *) klass;

  parent_class = erdos_expression_get_class();

  erdos_expression_class_init( expression_class );
}


void
erdos_presentation_init (ErdosPresentation *presentation)
{
  erdos_expression_init( (ErdosExpression *) presentation );
  ((ErdosPresentationClass *)(((ErdosObject *)presentation)->klass)) = erdos_presentation_get_class();
}

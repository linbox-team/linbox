/* -*- mode: c; style: linux -*- */

/* math-expression.h
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#ifndef __MATH_EXPRESSION_H
#define __MATH_EXPRESSION_H

#include <gnome.h>

#include "math-object.h"

BEGIN_GNOME_DECLS

#define MATH_EXPRESSION(obj)          GTK_CHECK_CAST (obj, math_expression_get_type (), MathExpression)
#define MATH_EXPRESSION_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, math_expression_get_type (), MathExpressionClass)
#define IS_MATH_EXPRESSION(obj)       GTK_CHECK_TYPE (obj, math_expression_get_type ())

typedef struct _MathExpression MathExpression;
typedef struct _MathExpressionClass MathExpressionClass;
typedef struct _MathExpressionPrivate MathExpressionPrivate;

struct _MathExpression 
{
	GtkObject parent;

	MathExpressionPrivate *p;
};

struct _MathExpressionClass 
{
	GtkObjectClass gtk_object_class;
};

guint math_expression_get_type            (void);

GtkObject *math_expression_new            (void);

MathObject *math_expression_get_toplevel  (MathExpression *expression);

END_GNOME_DECLS

#endif /* __MATH_EXPRESSION_H */

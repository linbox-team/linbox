/* -*- mode: c; style: linux -*- */

/* identifier.h
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

#ifndef __IDENTIFIER_H
#define __IDENTIFIER_H

#include <gnome.h>

#include "math-object.h"

BEGIN_GNOME_DECLS

#define IDENTIFIER(obj)          GTK_CHECK_CAST (obj, identifier_get_type (), Identifier)
#define IDENTIFIER_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, identifier_get_type (), IdentifierClass)
#define IS_IDENTIFIER(obj)       GTK_CHECK_TYPE (obj, identifier_get_type ())

typedef struct _Identifier Identifier;
typedef struct _IdentifierClass IdentifierClass;
typedef struct _IdentifierPrivate IdentifierPrivate;

struct _Identifier 
{
	MathObject parent;

	IdentifierPrivate *p;
};

struct _IdentifierClass 
{
	MathObjectClass math_object_class;
};

guint identifier_get_type         (void);

GtkObject *identifier_new         (void);

END_GNOME_DECLS

#endif /* __IDENTIFIER_H */

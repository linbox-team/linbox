/* -*- mode: c; style: linux -*- */

/* symbol.h
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

#ifndef __SYMBOL_H
#define __SYMBOL_H

#include <gnome.h>

#include "math-object.h"

BEGIN_GNOME_DECLS

#define SYMBOL(obj)          GTK_CHECK_CAST (obj, symbol_get_type (), Symbol)
#define SYMBOL_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, symbol_get_type (), SymbolClass)
#define IS_SYMBOL(obj)       GTK_CHECK_TYPE (obj, symbol_get_type ())

typedef struct _Symbol Symbol;
typedef struct _SymbolClass SymbolClass;
typedef struct _SymbolPrivate SymbolPrivate;

struct _Symbol 
{
	MathObject parent;

	SymbolPrivate *p;
};

struct _SymbolClass 
{
	MathObjectClass math_object_class;
};

guint symbol_get_type         (void);

GtkObject *symbol_new         (void);

END_GNOME_DECLS

#endif /* __SYMBOL_H */

/* -*- mode: c; style: linux -*- */

/* math-atom.h
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
 *
 * Caveat emptor! This class will silently fail if one tries to mix ISO-10346
 * 32-bit codes with 8-bit chars, since the `pos' argument in insert and
 * delete is used merely as an array index. So, please do not mix them!!!
 */

#ifndef __MATH_ATOM_H
#define __MATH_ATOM_H

#include <gnome.h>

#include "math-unit.h"

BEGIN_GNOME_DECLS

#define MATH_ATOM(obj)          GTK_CHECK_CAST (obj, math_atom_get_type (), MathAtom)
#define MATH_ATOM_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, math_atom_get_type (), MathAtomClass)
#define IS_MATH_ATOM(obj)       GTK_CHECK_TYPE (obj, math_atom_get_type ())

typedef struct _MathAtom MathAtom;
typedef struct _MathAtomClass MathAtomClass;
typedef struct _MathAtomPrivate MathAtomPrivate;

typedef enum _MathAtomType MathAtomType;

/* Atom types
 *
 * MATH_ATOM_NONE represents an uninitialized atom
 *
 * MATH_ATOM_DIVSTRING is a "divisible" string, i.e. the cursor can be placed
 * inside it, but it is still treated as indivisible by the formula
 * editor. Example: multi-digit numbers
 *
 * MATH_ATOM_STRING is an indivisible string, e.g. sin, cos, tan
 * 
 * MATH_ATOM_SINGLETON is a standalone variable or symbol, like x, y, the
 * integral sign, etc.
 */

enum _MathAtomType {
	MATH_ATOM_NONE, MATH_ATOM_DIVSTRING,
	MATH_ATOM_STRING, MATH_ATOM_SINGLETON
};

struct _MathAtom 
{
	MathUnit parent;

	MathAtomPrivate *p;
};

struct _MathAtomClass 
{
	MathUnitClass math_unit_class;
};

guint        math_atom_get_type       (void);

GtkObject   *math_atom_new            (MathAtomType type);

void         math_atom_set_code       (MathAtom *atom, gint code);

void         math_atom_append         (MathAtom *atom, gint code);
void         math_atom_insert         (MathAtom *atom, gint pos, gint code);
void         math_atom_delete         (MathAtom *atom, gint pos);

MathAtomType math_atom_get_atom_type  (const MathAtom *atom);
const gchar *math_atom_get_text       (const MathAtom *atom);
gint         math_atom_get_length     (const MathAtom *atom);

END_GNOME_DECLS

#endif /* __MATH_ATOM_H */

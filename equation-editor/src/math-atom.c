/* -*- mode: c; style: linux -*- */

/* math-atom.c
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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "math-atom.h"
#include "layout.h"
#include "glyph-layout.h"

enum {
	ARG_0,
	ARG_TYPE
};

struct _MathAtomPrivate 
{
	MathAtomType  type;
	gchar        *text;
	gint          length;
	gint          alloc_len;
};

static MathUnitClass *parent_class;

static void math_atom_init        (MathAtom *math_atom);
static void math_atom_class_init  (MathAtomClass *class);

static void math_atom_set_arg     (GtkObject *object, 
				   GtkArg *arg, 
				   guint arg_id);
static void math_atom_get_arg     (GtkObject *object, 
				   GtkArg *arg, 
				   guint arg_id);

static void math_atom_finalize    (GtkObject *object);

static void alloc_buffer          (MathAtom *atom, int req_len);

static Layout *math_atom_get_layout (MathObject *math_object);

guint
math_atom_get_type (void)
{
	static guint math_atom_type = 0;

	if (!math_atom_type) {
		GtkTypeInfo math_atom_info = {
			"MathAtom",
			sizeof (MathAtom),
			sizeof (MathAtomClass),
			(GtkClassInitFunc) math_atom_class_init,
			(GtkObjectInitFunc) math_atom_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		math_atom_type = 
			gtk_type_unique (math_unit_get_type (), 
					 &math_atom_info);
	}

	return math_atom_type;
}

static void
math_atom_init (MathAtom *math_atom)
{
	math_atom->p = g_new0 (MathAtomPrivate, 1);
	math_atom->p->type = MATH_ATOM_NONE;

	math_atom->p->alloc_len = 1;
	math_atom->p->text = g_new0 (gchar, 1);
}

static void
math_atom_class_init (MathAtomClass *class) 
{
	GtkObjectClass *object_class;
	MathObjectClass *math_object_class;

	gtk_object_add_arg_type ("MathAtom::type",
				 GTK_TYPE_INT,
				 GTK_ARG_READWRITE,
				 ARG_TYPE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = math_atom_finalize;
	object_class->set_arg = math_atom_set_arg;
	object_class->get_arg = math_atom_get_arg;

	math_object_class = MATH_OBJECT_CLASS (class);
	math_object_class->get_layout = math_atom_get_layout;

	parent_class = MATH_UNIT_CLASS
		(gtk_type_class (math_unit_get_type ()));
}

static void
math_atom_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathAtom *math_atom;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_ATOM (object));

	math_atom = MATH_ATOM (object);

	switch (arg_id) {
	case ARG_TYPE:
		g_return_if_fail (math_atom->p->type == MATH_ATOM_NONE);

		math_atom->p->type = GTK_VALUE_INT (*arg);
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
math_atom_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathAtom *math_atom;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_ATOM (object));

	math_atom = MATH_ATOM (object);

	switch (arg_id) {
	case ARG_TYPE:
		GTK_VALUE_INT (*arg) = math_atom->p->type;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
math_atom_finalize (GtkObject *object) 
{
	MathAtom *math_atom;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_ATOM (object));

	math_atom = MATH_ATOM (object);

	g_free (math_atom->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
math_atom_new (MathAtomType type) 
{
	return gtk_object_new (math_atom_get_type (),
			       "type", type,
			       NULL);
}

/**
 * math_atom_set_code:
 * @atom: 
 * @code: UTF-8 character code
 * 
 * Removes the existing text and creates a new string with just the code
 * given; used for changing the code in a singleton atom
 **/

void
math_atom_set_code (MathAtom *atom, gint code)
{
	g_return_if_fail (atom != NULL);
	g_return_if_fail (IS_MATH_ATOM (atom));

	while (atom->p->length > 0)
		math_atom_delete (atom, 0);

	math_atom_append (atom, code);

	gtk_signal_emit_by_name (GTK_OBJECT (atom), "changed", NULL);
}

/**
 * math_atom_append:
 * @atom: 
 * @code: 
 * 
 * Append the character with the given code to the end of the string
 **/

void
math_atom_append (MathAtom *atom, gint code)
{
	int req_len;

	g_return_if_fail (atom != NULL);
	g_return_if_fail (IS_MATH_ATOM (atom));

	req_len = (code & 0xfff0) ? 4 : 1;
	alloc_buffer (atom, req_len);

	if (req_len == 1)
		atom->p->text[atom->p->length] = (char) code;
	else
		*(int *) (atom->p->text + atom->p->length) = code;

	atom->p->length += req_len;
	atom->p->text[atom->p->length] = '\0';

	gtk_signal_emit_by_name (GTK_OBJECT (atom), "changed", NULL);
}

/**
 * math_atom_insert:
 * @atom: 
 * @pos: 
 * @code: 
 * 
 * Insert the character with the given code in the given position
 **/

void
math_atom_insert (MathAtom *atom, gint pos, gint code)
{
	gint req_len;

	g_return_if_fail (atom != NULL);
	g_return_if_fail (IS_MATH_ATOM (atom));

	req_len = (code & 0xfff0) ? 4 : 1;
	alloc_buffer (atom, req_len);

	memmove (atom->p->text + pos + req_len, atom->p->text + pos,
		 atom->p->length - pos);

	if (req_len == 1)
		atom->p->text[pos] = (char) code;
	else
		*(int *) (atom->p->text + pos) = code;

	atom->p->length += req_len;
	atom->p->text[atom->p->length] = '\0';

	gtk_signal_emit_by_name (GTK_OBJECT (atom), "changed", NULL);
}

/**
 * math_atom_delete:
 * @atom: 
 * @pos: 
 * 
 * Remove the character (assumed to be 8-bit) at the given position
 **/

void
math_atom_delete (MathAtom *atom, gint pos)
{
	g_return_if_fail (atom != NULL);
	g_return_if_fail (IS_MATH_ATOM (atom));
	g_return_if_fail (pos < atom->p->length);

	memmove (atom->p->text + pos, atom->p->text + pos + 1,
		 atom->p->length - pos);

	atom->p->length--;
	atom->p->text[atom->p->length] = '\0';

	gtk_signal_emit_by_name (GTK_OBJECT (atom), "changed", NULL);
}

/**
 * math_atom_get_atom_type:
 * @atom: 
 * 
 * Get the atom's type
 * 
 * Return value: Type, as defined in math-atom.h
 **/

MathAtomType 
math_atom_get_atom_type (const MathAtom *atom)
{
	g_return_val_if_fail (atom != NULL, 0);
	g_return_val_if_fail (IS_MATH_ATOM (atom), 0);

	return atom->p->type;
}

/**
 * math_atom_get_text:
 * @atom: 
 * 
 * Get the string of text associated with the atom
 * 
 * Return value: The text
 **/

const gchar *
math_atom_get_text (const MathAtom *atom)
{
	g_return_val_if_fail (atom != NULL, 0);
	g_return_val_if_fail (IS_MATH_ATOM (atom), 0);

	return atom->p->text;
}

/**
 * math_atom_get_length:
 * @atom: 
 * 
 * Get the length of the atom string
 * 
 * Return value: 
 **/

gint
math_atom_get_length (const MathAtom *atom)
{
	g_return_val_if_fail (atom != NULL, 0);
	g_return_val_if_fail (IS_MATH_ATOM (atom), 0);

	return atom->p->length;
}

static Layout *
math_atom_get_layout (MathObject *math_object)
{
	return LAYOUT (glyph_layout_new ());
}

/* Check buffer allocation to see if there is enough room for the current text
 * plus req_len and reallocate if necessary
 */

static void
alloc_buffer (MathAtom *atom, int req_len)
{
	if (atom->p->alloc_len < atom->p->length + req_len + 1) {
		atom->p->alloc_len *= 2;
		atom->p->text = 
			g_renew (gchar, atom->p->text, atom->p->alloc_len);
	}
}


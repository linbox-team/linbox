/* -*- mode: c; style: linux -*- */

/* number.h
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

#ifndef __NUMBER_H
#define __NUMBER_H

#include <gnome.h>

#include "unit.h"

BEGIN_GNOME_DECLS

#define NUMBER(obj)          GTK_CHECK_CAST (obj, number_get_type (), Number)
#define NUMBER_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, number_get_type (), NumberClass)
#define IS_NUMBER(obj)       GTK_CHECK_TYPE (obj, number_get_type ())

typedef struct _Number Number;
typedef struct _NumberClass NumberClass;
typedef struct _NumberPrivate NumberPrivate;

struct _Number 
{
	Unit parent;

	NumberPrivate *p;
};

struct _NumberClass 
{
	UnitClass unit_class;
};

guint      number_get_type    (void);

GtkObject *number_new         (void);

gdouble    number_get_value   (Number *number);
void       number_set_value   (Number *number, gdouble value);

END_GNOME_DECLS

#endif /* __NUMBER_H */

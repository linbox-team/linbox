/* -*- mode: c; style: linux -*- */

/* what.h
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

#ifndef __WHAT_H
#define __WHAT_H

#include <gnome.h>

#include ".h"

BEGIN_GNOME_DECLS

#define WHAT(obj)          GTK_CHECK_CAST (obj, what_get_type (), What)
#define WHAT_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, what_get_type (), WhatClass)
#define IS_WHAT(obj)       GTK_CHECK_TYPE (obj, what_get_type ())

typedef struct _What What;
typedef struct _WhatClass WhatClass;
typedef struct _WhatPrivate WhatPrivate;

struct _What 
{
	 parent;

	WhatPrivate *p;
};

struct _WhatClass 
{
	Class _class;
};

guint what_get_type         (void);

GtkObject *what_new         (void);

END_GNOME_DECLS

#endif /* __WHAT_H */

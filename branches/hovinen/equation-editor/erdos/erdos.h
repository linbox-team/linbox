/* erdos.h
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
#ifndef __ERDOS_H__
#define __ERDOS_H__

#ifdef __cplusplus
extern "C" {
#pragma }
#endif /* __cplusplus */

typedef struct _ErdosObject       ErdosObject;
typedef struct _ErdosObjectClass  ErdosObjectClass;

struct _ErdosObject
{
  ErdosObjectClass *klass;
  int ref_count;
};

struct _ErdosObjectClass
{
  void (*finalize) (ErdosObject *object);
  ErdosObject * (*clone) (ErdosObject *object);
};

void erdos_object_class_init	(ErdosObjectClass	 *klass);
void erdos_object_init		(ErdosObject		 *object);
ErdosObjectClass *erdos_object_get_class	(void);
void erdos_object_ref (ErdosObject *object);
void erdos_object_unref (ErdosObject *object);
ErdosObject *erdos_object_cow (ErdosObject *object);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __ERDOS_H__ */

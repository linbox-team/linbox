/* -*- mode: c; style: linux -*- */

/* cursor.c
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
 * The cursor maintains two stacks of selected math objects. The top of the
 * object stack is the object in which the cursor is located; the next lower
 * member is the object's parent, and so on. The top of the position stack is
 * the cursor's position in the current object; the next lower is the position
 * within the object's parent, and so on. It is an invariant that, for each
 * stack element except that at the bottom, the position of the stack element
 * (as returned by row_block_get_position_of) is equal to the position value
 * for the parent element.
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "cursor.h"
#include "row-block.h"
#include "fraction-block.h"
#include "math-atom.h"

enum {
	ARG_0,
	ARG_EXPR
};

enum {
	MOVED_SIGNAL,
	LAST_SIGNAL
};

struct _CursorPrivate 
{
	MathExpression *expr;

	GList *objects;
	GList *positions;
};

/* Evaluates to TRUE iff the cursor's given position is at the beginning of
 * the given object */

#define CURSOR_AT_BEGINNING(onode, pnode) \
            ((gint) pnode->data == 0)

static GtkObjectClass *parent_class;

static gint cursor_signals[LAST_SIGNAL] = { 0 };

static void cursor_init        (Cursor *cursor);
static void cursor_class_init  (CursorClass *class);

static void cursor_set_arg     (GtkObject *object, 
				GtkArg *arg, 
				guint arg_id);
static void cursor_get_arg     (GtkObject *object, 
				GtkArg *arg, 
				guint arg_id);

static void cursor_finalize    (GtkObject *object);

static gint get_object_length  (MathObject *object);
static void clean_objects_list (Cursor *cursor, 
				GList *node);
static GList *ascend           (Cursor *cursor,
				gboolean right,
				GList *onode,
				GList *pnode);
static gboolean descend        (Cursor *cursor,
				gboolean right,
				gboolean partial);
static gboolean is_navigable   (MathObject *object);
static gboolean cursor_at_end  (GList *onode, GList *pnode);

guint
cursor_get_type (void)
{
	static guint cursor_type = 0;

	if (!cursor_type) {
		GtkTypeInfo cursor_info = {
			"Cursor",
			sizeof (Cursor),
			sizeof (CursorClass),
			(GtkClassInitFunc) cursor_class_init,
			(GtkObjectInitFunc) cursor_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		cursor_type = 
			gtk_type_unique (gtk_object_get_type (), 
					 &cursor_info);
	}

	return cursor_type;
}

static void
cursor_init (Cursor *cursor)
{
	cursor->p = g_new0 (CursorPrivate, 1);
}

static void
cursor_class_init (CursorClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Cursor::expression",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_EXPR);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = cursor_finalize;
	object_class->set_arg = cursor_set_arg;
	object_class->get_arg = cursor_get_arg;

	cursor_signals[MOVED_SIGNAL] =
		gtk_signal_new ("moved", GTK_RUN_FIRST,
				object_class->type,
				GTK_SIGNAL_OFFSET (CursorClass, moved),
				gtk_signal_default_marshaller,
				GTK_TYPE_NONE, 0);

	gtk_object_class_add_signals (object_class, cursor_signals,
				      LAST_SIGNAL);

	parent_class = GTK_OBJECT_CLASS
		(gtk_type_class (gtk_object_get_type ()));
}

static void
cursor_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Cursor *cursor;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CURSOR (object));

	cursor = CURSOR (object);

	switch (arg_id) {
	case ARG_EXPR:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) != NULL);
		g_return_if_fail (IS_MATH_EXPRESSION 
				  (GTK_VALUE_POINTER (*arg)));

		if (cursor->p->expr != NULL)
			gtk_object_unref (GTK_OBJECT (cursor->p->expr));

		cursor->p->expr = GTK_VALUE_POINTER (*arg);

		if (cursor->p->expr != NULL) {
			gtk_object_ref (GTK_OBJECT (cursor->p->expr));
			cursor_move_to_beginning (cursor);
			g_assert (cursor->p->objects != NULL);
			g_assert (cursor->p->objects->data != NULL);
			g_assert (is_navigable 
				  (MATH_OBJECT (cursor->p->objects->data)));
			gtk_object_ref (GTK_OBJECT (cursor->p->objects->data));

			cursor->p->positions =
				g_list_prepend (NULL, (gpointer) 0);
		}
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
cursor_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Cursor *cursor;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CURSOR (object));

	cursor = CURSOR (object);

	switch (arg_id) {
	case ARG_EXPR:
		GTK_VALUE_POINTER (*arg) = cursor->p->expr;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
cursor_finalize (GtkObject *object) 
{
	Cursor *cursor;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CURSOR (object));

	cursor = CURSOR (object);

	g_list_foreach (cursor->p->objects, (GFunc) gtk_object_unref, NULL);
	g_list_free (cursor->p->objects);
	g_list_free (cursor->p->positions);
	g_free (cursor->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

/**
 * cursor_new:
 * @expr: 
 * 
 * Factory method
 * 
 * Return value: A new cursor object
 **/

GtkObject *
cursor_new (MathExpression *expr) 
{
	return gtk_object_new (cursor_get_type (),
			       "expression", expr,
			       NULL);
}

/**
 * cursor_get_current_object:
 * @cursor: object
 * 
 * Get the row block in which the user is currently navigating
 * 
 * Return value: The row block
 **/

MathObject *
cursor_get_current_object (Cursor *cursor)
{
	g_return_val_if_fail (cursor != NULL, NULL);
	g_return_val_if_fail (IS_CURSOR (cursor), NULL);

	gtk_object_ref (GTK_OBJECT (cursor->p->objects->data));
	return MATH_OBJECT (cursor->p->objects->data);
}

/**
 * cursor_get_insertion_point:
 * @cursor: 
 * 
 * Get the position of the insertion point in the current row block
 * 
 * Return value: The position of the insertion point
 **/

gint
cursor_get_insertion_point (Cursor *cursor)
{
	g_return_val_if_fail (cursor != NULL, 0);
	g_return_val_if_fail (IS_CURSOR (cursor), 0);

	return (gint) cursor->p->positions->data;
}

/**
 * cursor_get_object_at_insertion_point:
 * @cursor: 
 * 
 * Get the math object at which the insertion point is located; that is the
 * one to its immediate left
 * 
 * Return value: The math object
 **/

MathObject *
cursor_get_object_at_insertion_point (Cursor *cursor)
{
	MathObject *object;

	g_return_val_if_fail (cursor != NULL, NULL);
	g_return_val_if_fail (IS_CURSOR (cursor), NULL);

	if (IS_ROW_BLOCK (cursor->p->objects->data))
		object = row_block_get_object_at
			(ROW_BLOCK (cursor->p->objects->data),
			 (gint) cursor->p->positions->data);
	else
		object = cursor->p->objects->data;

	gtk_object_ref (GTK_OBJECT (object));
	return object;
}


/**
 * cursor_currently_in:
 * @cursor: 
 * @object: 
 * 
 * Check if the cursor is in the object given
 * 
 * Return value: TRUE if it is, FALSE otherwise
 **/

gboolean
cursor_currently_in (Cursor *cursor, MathObject *object)
{
	GList *node;

	g_return_val_if_fail (cursor != NULL, FALSE);
	g_return_val_if_fail (IS_CURSOR (cursor), FALSE);

	node = g_list_find (cursor->p->objects, object);
	if (node == NULL) return 0;

	while (node) {
		node = node->prev;
		if (node == NULL) break;
		if (IS_ROW_BLOCK (node->data)) return 0;
	}

	return 1;
}

/**
 * cursor_move_left:
 * @cursor: 
 * 
 * Move the cursor left one position, leaving it unchanged if it is at the
 * beginning of the math expression
 **/

void
cursor_move_left (Cursor *cursor)
{
	GList *node;

	g_return_if_fail (cursor != NULL);
	g_return_if_fail (IS_CURSOR (cursor));

	if (!descend (cursor, TRUE, FALSE)) {
		node = ascend (cursor, FALSE, 
			       cursor->p->objects, cursor->p->positions);
		if (node == NULL) return;

		if (node == cursor->p->objects || !IS_ROW_BLOCK (node->data)) {
			clean_objects_list (cursor, node);
			((gint) cursor->p->positions->data)--;
			descend (cursor, TRUE, TRUE);
		} else {
			clean_objects_list (cursor, node);
		}

		gtk_signal_emit (GTK_OBJECT (cursor), 
				 cursor_signals[MOVED_SIGNAL], NULL);
	}
}

/**
 * cursor_move_right:
 * @cursor: 
 * 
 * Move the cursor right one position, leaving it unchanged if it is at the
 * end of the math expression
 **/

void
cursor_move_right (Cursor *cursor)
{
	GList *node;

	g_return_if_fail (cursor != NULL);
	g_return_if_fail (IS_CURSOR (cursor));

	if (!descend (cursor, FALSE, FALSE)) {
		node = ascend (cursor, TRUE, 
			       cursor->p->objects, cursor->p->positions);
		if (node == NULL) return;

		if (node == cursor->p->objects || !IS_ROW_BLOCK (node->data)) {
			clean_objects_list (cursor, node);
			((gint) cursor->p->positions->data)++;
			descend (cursor, FALSE, TRUE);
		} else {
			clean_objects_list (cursor, node);
		}

		gtk_signal_emit (GTK_OBJECT (cursor), 
				 cursor_signals[MOVED_SIGNAL], NULL);
	}
}

/**
 * cursor_move_to_beginning:
 * @cursor: 
 * 
 * Positions the cursor at the beginning of the expression
 **/

void
cursor_move_to_beginning (Cursor *cursor)
{
	g_return_if_fail (cursor != NULL);
	g_return_if_fail (IS_CURSOR (cursor));

	clean_objects_list (cursor, NULL);
	cursor->p->objects = g_list_prepend
		(NULL, math_expression_get_toplevel (cursor->p->expr));
	cursor->p->positions =
		g_list_prepend (NULL, (gpointer) 0);
	descend (cursor, FALSE, FALSE);

	gtk_signal_emit (GTK_OBJECT (cursor), 
			 cursor_signals[MOVED_SIGNAL], NULL);
}

/**
 * cursor_move_to_end:
 * @cursor: 
 * 
 * Positions the cursor at the end of the expression
 **/

void
cursor_move_to_end (Cursor *cursor)
{
	g_return_if_fail (cursor != NULL);
	g_return_if_fail (IS_CURSOR (cursor));

	clean_objects_list (cursor, NULL);
	cursor->p->objects = 
		g_list_prepend (NULL,
				math_expression_get_toplevel
				(cursor->p->expr));
	cursor->p->positions =
		g_list_prepend (NULL, (gpointer) get_object_length 
				(MATH_OBJECT (cursor->p->objects->data)));
	descend (cursor, TRUE, FALSE);

	gtk_signal_emit (GTK_OBJECT (cursor), 
			 cursor_signals[MOVED_SIGNAL], NULL);
}

/**
 * get_object_length:
 * @math_object: 
 * 
 * Get the length of the object the cursor is currently in
 * 
 * Return value: 
 **/

static gint
get_object_length (MathObject *math_object) 
{
	g_return_val_if_fail (math_object != NULL, 0);
	g_return_val_if_fail (IS_MATH_OBJECT (math_object), 0);

	if (IS_ROW_BLOCK (math_object))
		return row_block_get_length (ROW_BLOCK (math_object));
/*  	else if (IS_MATH_ATOM (math_object)) */
/*  		return math_atom_get_length (MATH_ATOM (math_object)); */
	else if (IS_FRACTION_BLOCK (math_object))
		return 2;
	else
		g_assert_not_reached ();

	return 0;
}

/**
 * clean_objects_list:
 * @cursor: 
 * @node: The record up to which to clear the object and positions lists
 * 
 * Clear out the object and positions lists up to the record given; if node is
 * NULL, clears out the whole stack
 **/

static void
clean_objects_list (Cursor *cursor, GList *node) 
{
	GList *tmp;

	while (cursor->p->objects != node) {
  		gtk_object_unref (GTK_OBJECT (cursor->p->objects->data));
		tmp = cursor->p->objects;
		cursor->p->objects = cursor->p->objects->next;
  		g_list_free_1 (tmp);
		g_assert (cursor->p->objects != NULL);
		cursor->p->objects->prev = NULL;

		tmp = cursor->p->positions;
		cursor->p->positions = cursor->p->positions->next;
  		g_list_free_1 (tmp);
		g_assert (cursor->p->positions != NULL);
		cursor->p->positions->prev = NULL;
	}
}

/**
 * ascend:
 * @cursor: object
 * @right: TRUE if it should ascend the tree until the cursor may move to the
 * right; false if it should ascend the tree until the cursor may move to the
 * left
 * @onode: Internal data passed through the stack; initialize to the current
 * top of the objects stack
 * @pnode: Internal data passed through the stack; initialize to the current
 * top of the positions stack
 * 
 * Ascends the tree until the cursor can move left or right, depending on the
 * value of the flag @right
 *
 * Return value: Returns the new intended top of the stack, or NULL if the
 * cursor cannot be moved to the right
 **/

static GList *
ascend (Cursor *cursor, gboolean right, GList *onode, GList *pnode) 
{
	if (IS_ROW_BLOCK (onode->data) && onode != cursor->p->objects)
		return onode;

	/* If we can move now, do so and return */
	if (is_navigable (MATH_OBJECT (onode->data)) &&
	    ((right && !cursor_at_end (onode, pnode)) ||
	     (!right && !CURSOR_AT_BEGINNING (onode, pnode))))
		return onode;

	/* Otherwise, go up a level and try again */
	if (onode->next == NULL) return NULL; /* Can't move left */
	g_assert (pnode->next != NULL);
	return ascend (cursor, right, onode->next, pnode->next);
}

/**
 * descend:
 * @cursor: 
 * @right: TRUE if it should descend to the last possible object, FALSE
 * otherwise
 * @partial: TRUE iff it should only descend enough to make the system
 * consistent again
 * 
 * Descends to the lowest possible math object that it can, either on the left
 * side or on the right side
 *
 * Return value: TRUE if the descent was possible, FALSE otherwise
 **/

static gboolean
descend (Cursor *cursor, gboolean right, gboolean partial) 
{
	MathObject *ins_point;

	g_return_val_if_fail (cursor != NULL, FALSE);
	g_return_val_if_fail (IS_CURSOR (cursor), FALSE);

	if (!partial && IS_ROW_BLOCK (cursor->p->objects->data)) {
		ins_point = row_block_get_object_at
			(ROW_BLOCK (cursor->p->objects->data),
			 (gint) cursor->p->positions->data);

		if (is_navigable (ins_point))
			cursor->p->objects =
				g_list_prepend (cursor->p->objects, ins_point);
		else
			return FALSE;
	}
	else if (IS_FRACTION_BLOCK (cursor->p->objects->data)) {
		cursor->p->objects =
			g_list_prepend (cursor->p->objects,
					(gpointer) 
					cursor->p->positions->data == 0 ?
					fraction_block_get_numerator
					(FRACTION_BLOCK 
					 (cursor->p->objects->data)) :
					fraction_block_get_denominator
					(FRACTION_BLOCK
					 (cursor->p->objects->data)));
	} else {
		return FALSE;
	}

	if (is_navigable (MATH_OBJECT (cursor->p->objects->data)))
		cursor->p->positions =
			g_list_prepend (cursor->p->positions,
					(gpointer) (right ? 
					get_object_length
					(cursor->p->objects->data) - 1 : 0));
	else
		cursor->p->positions =
			g_list_prepend (cursor->p->positions, (gpointer) 0);

	if (!IS_ROW_BLOCK (cursor->p->objects->data))
		descend (cursor, right, partial);

	return TRUE;
}

/**
 * is_navigable:
 * @node: Object to check
 * 
 * Determine if the given object is internally navigable
 * 
 * Return value: TRUE if it is, FALSE otherwise
 **/

static gboolean
is_navigable (MathObject *object) 
{
#if 0
	MathAtomType type;
#endif

	if (IS_ROW_BLOCK (object)) return TRUE;
	if (IS_FRACTION_BLOCK (object)) return TRUE;

#if 0
	if (IS_MATH_ATOM (object)) {
		type = math_atom_get_atom_type (MATH_ATOM (object));
		if (type == MATH_ATOM_DIVSTRING) return TRUE;
	}
#endif

	return FALSE;
}

/**
 * cursor_at_end:
 * @onode: 
 * @pnode: 
 * 
 * Returns TRUE iff the cursor is in the last valid position of the block
 * 
 * Return value: 
 **/

static gboolean
cursor_at_end (GList *onode, GList *pnode) 
{
	if (IS_ROW_BLOCK (onode->data))
		return (gint) pnode->data == 
			get_object_length (MATH_OBJECT (onode->data));
	else
		return (gint) pnode->data == 
			get_object_length (MATH_OBJECT (onode->data)) - 1;
}

/* File: src/examples/field/field_categories.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _FIELD_CATEGORIES_
#define _FIELD_CATEGORIES_

/** List of vector categories.
 * This structure contains nineteen structures: one relating each field type
 * to be tested plus one for all fields combined.  
 * These types allow us to use template specialization to use different 
 * code for different fields.
 */
struct field_categories
{
  struct all_tag {};
  struct double_tag {};
  struct float_tag {};
  struct modular_tag {};
  struct fuzzy_tag {};
  struct param_modular_tag {};
  struct param_fuzzy_tag {};
  struct double_envelope_tag {};
  struct float_envelope_tag {};
  struct modular_envelope_tag {};
  struct fuzzy_envelope_tag {};
  struct param_modular_envelope_tag {};
  struct param_fuzzy_envelope_tag {};
  struct abstract_double_tag {};
  struct abstract_float_tag {};
  struct abstract_modular_tag {};
  struct abstract_fuzzy_tag {};
  struct abstract_param_modular_tag {};
  struct abstract_param_fuzzy_tag {};
}; // struct field_categories

#endif // _FIELD_CATEGORIES_


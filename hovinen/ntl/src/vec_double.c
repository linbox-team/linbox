
#include <NTL/vec_double.h>

inline void BlockConstruct(double *, double) { }
inline void BlockDestroy(double *, double) { }

NTL_vector_impl_plain(double,vec_double)

NTL_io_vector_impl(double,vec_double)

NTL_eq_vector_impl(double,vec_double)


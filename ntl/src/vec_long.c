
#include <NTL/vec_long.h>

inline void BlockConstruct(long *, long) { }
inline void BlockDestroy(long *, long) { }

NTL_vector_impl_plain(long,vec_long)

NTL_io_vector_impl(long,vec_long)

NTL_eq_vector_impl(long,vec_long)


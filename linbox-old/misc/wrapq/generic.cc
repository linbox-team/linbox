/* separately compile this one */
#include "base.h"

domainbase::element& genericalg(const domainbase& D, domainbase::element& b, const domainbase::element& a)
{ return D.mul(b, a, a); }


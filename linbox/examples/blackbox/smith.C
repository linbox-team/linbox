/** @name examples/blackbox/ju.C
 * @author bds
 *
 * @doc
 */
//@{

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <utility>

#include "linbox/util/timer.h"
#include "linbox/field/ntl-ZZ.h"

using std::cout;
using std::endl;

template<class I1, class Lp>
void distinct (I1 a, I1 b, Lp& c);

template <class I>
void display(I b, I e);

int main (int argc, char** argv) {

	typedef LinBox::NTL_ZZ Field;

	Field F;

	int n = atoi (argv[1]);

	std::vector<Field::Element> l (n);

	std::vector<Field::Element>::iterator p1, p2;

	int i = 1;

	for (p1 = l. begin(); p1 != l. end(); ++ p1, ++ i) 

		*p1 = i;


	Field::Element tmp1, tmp2;

	for (p1 = l .begin(); p1 != l. end(); ++ p1) 

		for (p2 = p1 + 1; p2 != l. end(); ++ p2) {
			
			F. gcd (tmp1, *p2, *p1);

			F. lcm (tmp2, *p2, *p1);

			*p1 = tmp1;

			*p2 = tmp2;

		}
	

	std::list<std::pair<Field::Element, int> >L;

	distinct (l. begin(), l. end(), L);

	std::cout << "Total distinct invariant factors: " << L. size() << '\n';

	display (L. begin(), L. end());

	return 0;
}

template<class I1, class Lp>
void distinct (I1 a, I1 b, Lp& c)
{ typename I1::value_type e;
  size_t count = 0;
  if (a != b) {e = *a; ++a; count = 1;}
  else return;
  while (a != b)
  {  if (*a == e) ++count;
     else
     { c.push_back(typename Lp::value_type(e, count));
       e = *a; count = 1;
     }
     ++a;
  }
  c.push_back(typename Lp::value_type(e, count));
  return;
}

template <class I>
void display(I b, I e)
{ cout << "(";
  for (I p = b; p != e; ++p) cout << p->first << " " << p->second << ", ";
  cout << ")" << endl;
}

	

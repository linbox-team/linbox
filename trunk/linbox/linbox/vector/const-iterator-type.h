#ifndef CONSTITERATORTYPE_H
#define CONSTITERATORTYPE_H
#include <vector>

namespace LinBox
{
  template<class Iterator>
    class Subiterator;

  template<class Iterator>
    class ConstIteratorType
    {
    public:
      typedef Iterator const_iterator;
    };
  
  template<class Iterator>
    class ConstIteratorType<Subiterator<Iterator> >
    {
    public:
      typedef Subiterator<ConstIteratorType<Iterator>::const_iterator> const_iterator;
    };

  template<class T>
    class ConstIteratorType<T* >
    {
    public:
      typedef const T* const_iterator;
    };

  template<>
    class ConstIteratorType<vector<bool>::iterator >
    {
    public:
      typedef vector<bool>::const_iterator const_iterator;
    };
}
#include <linbox/vector/subiterator.h>
#endif

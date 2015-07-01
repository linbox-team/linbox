/* File: src/library/objects/auxillary/cekstv_switch.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _CEKSTV_SWITCH_
#define _CEKSTV_SWITCH_

#ifdef TRACE
#include <iostream>
#endif // TRACE

#include <vector>

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Butterfly switch object from preconditioner paper.
   * This is a switch predicate object that is applied
   * to two references to elements to switch them as needed 
   * by the \Ref{Butterfly Switching Network BlackBox Matrix Object}
   * following the exchange matrix introduced in "Efficient Matrix
   * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
   * Kaltofen, Saunders, Turner, and Villard.
   * This class is templatized by the field in which the arithmetic
   * is done.
   */
  template <class Field>
  class cekstv_switch
  {
  public:

    /// Typedef
    typedef typename Field::element element;

    /** Constructor from a field and random field element generator.
     * The switch is applied using the vector of field elements.
     * If the current switch is marked by the field element a,
     * and the two elements are x and y, the output from the switch
     * is x' = x + a*y and y' = y + (x + a*y).
     * The generator is used to create random field elements for setting the 
     * switches.  Both apply function and applyTranspose function 
     * use the random elements in the order they are generated.  Neither 
     * can re-use random elements.  This means each time the matrix (or
     * its transpose) is applied it will be a different matrix.
     * @param F field in which arithmetic is done
     * @param R random field element generator
     */
    cekstv_switch(const Field& F, const typename Field::randIter& R);

    /** Constructor from a field and STL vector of field elements.
     * The switch is applied using the vector of field elements.
     * If the current switch is marked by the field element a,
     * and the two elements are x and y, the output from the switch
     * is x' = x + a*y and y' = y + (x + a*y).
     * The apply function starts at the beginning of the vector moving 
     * forward through it, and applyTranspose function starts at the end
     * moving backwards.  Both repeat the vector after they pass through it.
     * @param F field in which arithmetic is done
     * @param switches vector of switches
     */
    cekstv_switch(const Field& F, const std::vector<element>& switches);

    /** Destructor.
     */
    ~cekstv_switch(void) {}

    /** Apply switch function.
     * Switches the elements in references according to the
     * exchange matrix introduced in "Efficient Matrix
     * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
     * Kaltofen, Saunders, Turner, and Villard and the current field element
     * specified in the switch object.
     * @return bool true if swapped, false otherwise
     * @param x reference to first element to be switched
     * @param y reference to second element to be switched
     */
    bool apply(element& x, element& y);

    /** Apply switch transpose function.
     * Switches the elements in references according to the
     * transpose of the exchange matrix introduced in "Efficient Matrix
     * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
     * Kaltofen, Saunders, Turner, and Villard and the current field element
     * specified in the switch object.
     * @return bool true if swapped, false otherwise
     * @param x reference to first element to be switched
     * @param y reference to second element to be switched
     */
    bool applyTranspose(element& x, element& y);

  private:

    // Field in which arithemetic is done
    Field _F;

    // Random field element generator;
    typename Field::randIter _R;

    // STL vector of boolean flags for switches
    std::vector<element> _switches;

    // STL vector iterator and reverse iterator pointing to current switch
    // and its transpose
    std::vector<element>::iterator _iter;
    std::vector<element>::reverse_iterator _riter;

    // temporary field element used in arithmetic
    element _temp;
    
  }; // class cekstv_switch

  template <class Field>
  inline cekstv_switch<Field>
  ::cekstv_switch(const Field& F, const typename Field::randIter& R)
  : _F(F), _R(R)
  { 
#ifdef TRACE
    clog
      << "Called cekstv_switch constructor from random field element generator."
      << endl;
#endif // TRACE

    _iter = _switches.begin();
    _F.init(_temp, 0);
    
  } // cekstv_switch::cekstv_switch(F, R)

  template <class Field>
  inline cekstv_switch<Field>
  ::cekstv_switch(const Field& F,
		  const std::vector<typename Field::element>& switches)
  : _F(F), _R(_F), _switches(switches)
  { 
#ifdef TRACE
    clog
      << "Called cekstv_switch constructor from STL vector." << endl
      << "constructucted switch vector:" << endl
      << "    i        switches[i]" << endl
      << "    --------------------" << endl;
    
    for (size_t i = 0; i < _switches.size(); i++)
    {
      clog << "    " << i << "    ";
      _F.write(clog, _switches[i]);
      clog << endl;
    } // for (size_t i = 0; i <= _switches.size(); i++)
      
#endif // TRACE

    _iter = _switches.begin(); 
    _F.init(_temp, 0);

  } // cekstv_switch::cekstv_switch(F, switches)

  template <class Field> 
  inline bool cekstv_switch<Field>::apply(typename Field::element& x,
					  typename Field::element& y)
  {
    if (_switches.empty())
    {
#ifdef TRACE
      typename Field::element oldx(x), oldy(y), s(x);
#endif // TRACE

      _F.addin(x, _F.mul(_temp, 
#ifdef TRACE
s =
#endif // TRACE
			 _R(), y));
      _F.addin(y, x);
      
#ifdef TRACE
      clog << "Switched ";
      _F.write(clog, oldx);
      clog << " and ";
      _F.write(clog, oldy);
      clog << " with switch value ";
      _F.write(clog, s);
      clog << " to obtain ";
      _F.write(clog, x);
      clog << " and ";
      _F.write(clog, y);
      clog << endl;
#endif // TRACE

    } // if (_switches.empty())
    else
    {    
      // If at end of vector, repeat it
      if (_iter == _switches.end()) _iter = _switches.begin();

#ifdef TRACE
      typename Field::element oldx(x), oldy(y), s(x);
#endif // TRACE
      
      _F.addin(x, _F.mul(_temp, 
#ifdef TRACE
s =
#endif // TRACE
			 *_iter++, y));
      _F.addin(y, x);
      
#ifdef TRACE
      clog << "Switched ";
      _F.write(clog, oldx);
      clog << " and ";
      _F.write(clog, oldy);
      clog << " with switch value ";
      _F.write(clog, s);
      clog << " to obtain ";
      _F.write(clog, x);
      clog << " and ";
      _F.write(clog, y);
      clog << endl;
#endif // TRACE

    } // else
    
    // return with swap flag
    return true;

  } // bool cekstv_switch::apply(Element& x, Element& y)
  
  template <class Field> 
  inline bool cekstv_switch<Field>::applyTranspose(typename Field::element& x,
						   typename Field::element& y)
  {
    if (_switches.empty())
    {
#ifdef TRACE
      typename Field::element oldx(x), oldy(y), s(x);
#endif // TRACE
      
      _F.addin(x, _F.mul(_temp, 
#ifdef TRACE
s =
#endif // TRACE
			 _R(), y));
      _F.addin(y, x);
      
#ifdef TRACE
      clog << "Switched ";
      _F.write(clog, oldx);
      clog << " and ";
      _F.write(clog, oldy);
      clog << " with switch value ";
      _F.write(clog, s);
      clog << " to obtain ";
      _F.write(clog, x);
      clog << " and ";
      _F.write(clog, y);
      clog << endl;
#endif // TRACE

    } // if (_switches.empty())
    else
    {    
      // If at end of vector, extend it
      if (_riter == _switches.rend()) _riter = _switches.rbegin();

#ifdef TRACE
      typename Field::element oldx(x), oldy(y), s(x);
#endif // TRACE
      
      _F.addin(x, _F.mul(_temp, 
#ifdef TRACE
s =
#endif // TRACE
			 *_riter++, y));
      _F.addin(y, x);
      
#ifdef TRACE
      clog << "Switched ";
      _F.write(clog, oldx);
      clog << " and ";
      _F.write(clog, oldy);
      clog << " with switch value ";
      _F.write(clog, s);
      clog << " to obtain ";
      _F.write(clog, x);
      clog << " and ";
      _F.write(clog, y);
      clog << endl;
#endif // TRACE

    } // else
    
    // return with swap flag
    return true;

  } // bool cekstv_switch::applyTranspose(Element& x, Element& y)

 } // namespace LinBox

#endif // _CEKSTV_SWITCH_

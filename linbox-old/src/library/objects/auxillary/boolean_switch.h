/* File: src/library/objects/auxillary/boolean_switch.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _BOOLEAN_SWITCH_
#define _BOOLEAN_SWITCH_

#ifdef TRACE
#include <iostream>
#endif // TRACE

#include <vector>

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Boolean switch object.
   * This is a switch predicate object that is applied
   * to two references to elements to switch them as needed 
   * by the \Ref{Butterfly Switching Network BlackBox Matrix Object}.
   */
  class boolean_switch
  {
  public:

    /** Constructor from an STL vector of booleans.
     * The switch is applied using the vector of booleans.
     * A true value means to swap the two elements, and a false
     * value means not to.
     * The apply function starts at the beginning of the vector moving 
     * forward through it, and applyTranspose function starts at the end
     * moving backwards.  Both repeat the vector after they pass through it.
     * @param switches vector of switches
     */
    boolean_switch(const std::vector<bool>& switches);

    /** Destructor.
     */
    ~boolean_switch(void) {}

    /** Apply switch function.
     * Switches the elements in references according to current boolean
     * value.  Swaps the elements if boolean is true, otherwise does nothing.
     * It is templatized by the element type to be swapped.
     * @return bool true if swapped, false otherwise
     * @param x reference to first element to be switched
     * @param y reference to second element to be switched
     */
    template <class Element>
    bool apply(Element& x, Element& y);

    /** Apply switch transpose function.
     * Switches the elements in references according to current boolean
     * value.  Swaps the elements if boolean is true, otherwise does nothing.
     * It is templatized by the element type to be swapped.
     * @return bool true if swapped, false otherwise
     * @param x reference to first element to be switched
     * @param y reference to second element to be switched
     */
    template <class Element>
    bool applyTranspose(Element& x, Element& y);

  private:

    // STL vector of boolean flags for switches
    std::vector<bool> _switches;

    // STL vector iterator and reverse iterator pointing to current switch
    // and its transpose
    std::vector<bool>::iterator _iter;
    std::vector<bool>::reverse_iterator _riter;

  }; // class boolean_switch

  inline boolean_switch::boolean_switch(const std::vector<bool>& switches)
  : _switches(switches)
  { 
#ifdef TRACE
    clog
      << "Called boolean_switch constructor from STL vector." << endl
      << "constructucted switch vector:" << endl
      << "    i        switches[i]" << endl
      << "    --------------------" << endl;
    
    for (size_t i = 0; i < _switches.size(); i++)
      clog << "    " << i << "    " << _switches[i] << endl;
      
#endif // TRACE

    _iter = _switches.begin(); 
    _riter = _switches.rbegin(); 
  } // boolean_switch::boolean_switch(const std::vector<bool>& switches)

  template <class Element> 
  inline bool boolean_switch::apply(Element& x, Element& y)
  {
    // If at end of vector, extend it
    if (_iter == _switches.end())
    {
      if (_switches.size() == 0)
	return false;
      else   
	_iter = _switches.begin();
    } // if (_iter == _switches.end())

    // If true value, swap
    if (*_iter)
    {
      Element temp(x);
      x = y;
      y = temp;
      
#ifdef TRACE
      clog << "Swapped elements" << endl;
#endif // TRACE

     } // if (*_iter)
#ifdef TRACE
    else
      clog << "Did not swap elements" << endl;
#endif // TRACE

    // return with flag of *_iter and update _iter
    return *_iter++;

  } // bool boolean_switch::apply(Element& x, Element& y)
  
  template <class Element> 
  inline bool boolean_switch::applyTranspose(Element& x, Element& y)
  {
    // If at end of vector, extend it
    if (_riter == _switches.rend())
    {
      if (_switches.size() == 0)
	return false;
      else   
	_riter = _switches.rbegin();
    } // if (_iter == _switches.rend())

    // If true value, swap
    if (*_riter)
    {
      Element temp(x);
      x = y;
      y = temp;
      
#ifdef TRACE
      clog << "Swapped elements" << endl;
#endif // TRACE

     } // if (*_riter)
#ifdef TRACE
    else
      clog << "Did not swap elements" << endl;
#endif // TRACE

    // return with flag of *_iter and update _iter
    return *_riter++;

  } // bool boolean_switch::applyTranspose(Element& x, Element& y)
  
 } // namespace LinBox

#endif // _BOOLEAN_SWITCH_

#ifndef __LINBOX_METHOD_TRAITS_H__
#define __LINBOX_METHOD_TRAITS_H__


#ifndef _DEFAULT_EarlyTerm_THRESHOLD_
#define _DEFAULT_EarlyTerm_THRESHOLD_ 20
#endif


struct WiedemannTraits {

  WiedemannTraits(unsigned long thres = _DEFAULT_EarlyTerm_THRESHOLD_) : _ett(thres) {}
  unsigned long Early_Term_Threshold( ) { return _ett; }

  private:
  unsigned long _ett;

};




struct MethodTrait {

  typedef WiedemannTraits Wiedemann;

};

#endif

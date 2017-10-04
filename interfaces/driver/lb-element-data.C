#include "lb-element-data.h"
#include "lb-domain-function.inl"

EltAbstract* constructElt(const DomainKey &key){
	EltAbstract *e;
	CreateEltFunctor Fct(key);
	DomainFunction::call(e, key, Fct);
	return e;
}
 

#include "lb-element-data.h"

EltAbstract* constructElt(const DomainKey &key){
	EltAbstract *e;
	CreateEltFunctor Fct(key);
	DomainFunction::call(e, key, Fct);
	return e;
}
 

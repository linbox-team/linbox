/* abstract base types */
class objectbase
{public:
  // basic object interface

  /// copy constructor - makes an instance initialized from b.
  virtual void cstor(const objectbase& b) =0;

  // Derived classes should provide a no-arg constructor which allocates the memory
  // in such a way that subsequent assignment to the object is valid.

  /// assignment operator - modifies an instance by copying from b.
  virtual objectbase& operator=(const objectbase& b) =0;
};

class domainbase : public objectbase
{public:
  typedef objectbase element;
  //typedef void* element;
  //class element;
  virtual element& mul(element& r, const element& a, const element& b) const=0;
};


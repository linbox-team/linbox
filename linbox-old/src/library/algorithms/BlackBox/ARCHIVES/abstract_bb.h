// ==========================================================================
// $Source$
// Copyright(c)'99 by LINBOX Team
// Authors: Erich Kaltofen, Thierry Gautier
// $Id$
// Description:
// - Abstract BB for the LinBox project. 
// One key point is the allocation of objects through the base classes
// BB_Base and BB_Base::Entry. Currently, the memory manager for the 
// concrete objects use killclone and clone provided in the set of 
// virtual methods. BB_Abstract and BB_Abstract::Entry are envelop for
// BB_Base* and BB_Base::Entry* and use these two methods to provide 
// a "copy by value" semantic.
// ==========================================================================
#ifndef _LINBOX_ABSTRACT_BB_H_
#define _LINBOX_ABSTRACT_BB_H_


// --
// - Base class from which derived concrete BB classes.
// - Memory management routines are given through virtual
// - clone and killclone method. 
// --
class BB_Base {
public:
  // - Base class from which derives concret class for the entries of apply
  class Entry{
  public:
    virtual ~Entry() {}

    // - make a copy of the object
    virtual Entry* clone() const =0;

    // - delete a copy make by clone: (delete this ???)
    virtual void killclone() =0;

    // - print an entry
    virtual ostream& write(ostream& o) const { return o; }
  };

  // - Dstor
  virtual ~BB_Base() {}

  // - Make a copy of the BB_Base object 
  virtual BB_Base* clone() const =0;

  // - Delete a copy make by clone: (delete this ???)
  virtual void killclone() =0;

  // - Apply method: should take pointers to entry, may be pointer ==0 (un-allocated !)
  virtual void apply( Entry*& res, const Entry* const& src ) =0;
};


// --
// - BB_Abstract: abstract BB 
// - This class manages BB_Base object through a pointer using a "copy by value"
// - semantic (assignement operation copies the value). There are no sharing of
// - data.
// - Construction of a new BB_Abstract needs a pointer to BB_Base that could
// - be deallocated using killclone.
// --
class BB_Abstract {
public:
  // - Inner type for both prefered in and out in template version
  class Entry {
  public:
    // - Default Cstor:
    Entry( BB_Base::Entry* entry =0) : _entry(entry) {};

    // - Entry recopy constructor
    Entry( const Entry& E ) : _entry(0) 
    { if (E._entry !=0) _entry = E._entry->clone(); }

    // - Assignement
    Entry& operator=( const Entry& E ) 
    { 
      if (_entry !=0) _entry->killclone();
      if (E._entry !=0) _entry = E._entry->clone(); 
    }

    // - DCstor:
    ~Entry() { if (_entry !=0) _entry->killclone(); _entry=0; }

    friend ostream& operator<< (ostream& o, const Entry& e )
    { if (e._entry !=0) e._entry->print(o); 
      return o;
    }
  protected:
    BB_Base::Entry* _entry;
    friend class BB_Abstract;
  };

  // - BB_Abstract constructor: if bb!=0, it should be deleted by killclone.
  BB_Abstract( BB_Base* bb =0) : _bb(bb) {};

  // - BB_Abstract recopy constructor
  BB_Abstract( const BB_Abstract& BA ) : _bb(0) 
  {
    if (BA._bb !=0) _bb = BA._bb->clone();
  }

  // - Assignement
  BB_Abstract& operator=( const BB_Abstract& BA ) 
  { 
    if (_bb !=0) { _bb->killclone(); _bb =0; }
    if (BA._bb !=0) _bb = BA._bb->clone(); 
  }

  // - BB_Abstract destructor
  ~BB_Abstract( ) {
    if (_bb !=0) _bb->killclone();
  }

  // - Ok !
  void apply( Entry& res, const Entry& src )
  {
    // if null pointer then crash: ok
    _bb->apply( res._entry, src._entry ); 
  }

protected:
  BB_Base* _bb;
};

#endif

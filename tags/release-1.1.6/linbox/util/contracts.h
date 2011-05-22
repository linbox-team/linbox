// Object class is the base class for all
// objects in the system. All classes inheriting from this class need 
// to define a method IsValid. This method should perform a
// consistency check on the state of the object. Note that 
// this method needs to be defined only when a debug build is made
#include <cassert>
#include <cstddef>
/*
class Object
{
public:
#ifdef _DEBUG
    virtual bool IsValid() const = 0;
#endif
    
};
*/
#if _DEBUG == 2

// The debug mode also defines the following macros. Failure of any of these macros leads to
// program termination. The user is notified of the error condition with the right file name
// and line number. The actual failing operation is also printed using the stringizing operator #

#define ASSERT(bool_expression) assert(bool_expression)
#define IS_VALID(obj) ASSERT((obj) != NULL && (obj)->IsValid())
#define REQUIRE(bool_expression) ASSERT(bool_expression)
#define ENSURE(bool_expression) ASSERT(bool_expression)
#define STATE(expression) expression

#elif _DEBUG == 1

#define ASSERT(bool_expression) assert(bool_expression)
#define IS_VALID(obj) ASSERT((obj) != NULL && (obj)->IsValid())
#define REQUIRE(bool_expression) if( ! bool_expression ) throw InvalidOperationException(__FILE__,__LINE__);
#define ENSURE(bool_expression) ASSERT(bool_expression)
#define STATE(expression) expression

#elif _DEBUG == 0

// When built in release mode, the _DEBUG flag would not be defined, thus there will be no overhead
// in the final release from these checks.

#define ASSERT(ignore) ((void) 0)
#define IS_VALID(ignore) ((void) 0) 
#define REQUIRE(ignore) ((void) 0)
#define ENSURE(ignore) ((void) 0)
#define STATE(ignore) ((void) 0)

#endif

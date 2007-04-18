#ifndef _KAAPI_ALGORITHM_H_
#define _KAAPI_ALGORITHM_H_

#include "communicate.h"
#include <athapascan-1>

// forward declarations   
namespace LinBox {

template < class _A, class _B> struct IntegerModularDet;
template < class _A, class _B> struct Residue;

}

namespace a1 {

/**************************************************************************************************************************
 * transform algorithm, similar to std::transform with unary operator
 * a generic algorithm that supposes that the unary opearator is communicable is avaailable in kaapi
 * - the best way to use them is to make your operator communicable
 * - there is always the possibility to specialize the algorithm for special operators, as we do for the Residue algorithm
 */

/*******************************************************************************************************
 * specialization for Residue and IntegerModularDet
 */
template <class InputIterator, class OutputIterator, class Blackbox, class Method, class Domain>
struct internal_transform<InputIterator, OutputIterator,  LinBox::Residue<LinBox::IntegerModularDet<Blackbox,Method>,Domain> >
{
    typedef typename std::iterator_traits<InputIterator>::difference_type diff_t;
    typedef LinBox::IntegerModularDet<Blackbox,Method> Function;
    typedef LinBox::Residue< Function, Domain > UnaryOperator;

    void operator()( const_remote<InputIterator> ibeg, const_remote<InputIterator> iend, remote<OutputIterator> obeg, remote<OutputIterator> oend,  std::string str) {
        const diff_t size = iend - ibeg;
        if( size < 10/*floor*/ ) { 
            fetch(ibeg,iend);
            fetch(obeg,oend);
            // create the black box from the string
            typename Blackbox::Field field;
            Blackbox bb(field);
            std::istringstream iss(str);
            bb.read(iss);

            Method m;

            // init IntegerModularDet
            Function f(bb,m);
            UnaryOperator op(f);

            // apply the operator
            std::transform(ibeg,iend,obeg,op);
        }   
        else {
            const diff_t size_2 = size /2;
            Fork<internal_transform< InputIterator, OutputIterator, UnaryOperator> > () ( ibeg, ibeg + size_2, obeg, obeg + size_2, str );
            Fork<internal_transform< InputIterator, OutputIterator, UnaryOperator> > () ( ibeg + size_2, iend, obeg + size_2, oend, str );
        }   
    }
};

template <class InputIterator, class OutputIterator, class Blackbox, class Method, class Domain, class Attribut>
struct transform_task < InputIterator,  OutputIterator,  LinBox::Residue<LinBox::IntegerModularDet<Blackbox,Method>,Domain>, Attribut >
{
    typedef LinBox::IntegerModularDet<Blackbox,Method> Function;
    typedef LinBox::Residue< Function, Domain > UnaryOperator;

    void operator()( InputIterator begin, InputIterator end, OutputIterator to_fill , LinBox::Residue<LinBox::IntegerModularDet<Blackbox,Method>,Domain> op, Attribut att )
    {
        typedef  typename std::iterator_traits<OutputIterator>::value_type OutputType;
        typedef  typename std::iterator_traits<InputIterator>::value_type InputType;

        // serialize the blackbox in a stream
        std::stringstream os;
        op.f.A.write(os, LinBox::FORMAT_GUILLAUME);

        /* init remote iterator */
        const_remote<InputIterator> ibeg, iend;
        init(ibeg, iend, begin,end);

        remote<OutputIterator> obeg,oend;
        init(obeg,oend, to_fill, to_fill + (end-begin) );

        /* add a frame to local thread */
        DFG::Thread* thread = dynamic_cast<DFG::Thread*>(Core::Thread::get_current());
        DFG::Frame frame;
        thread->push(&frame);

        /* fill the frame */
        {
            a1::Fork< internal_transform< InputIterator, OutputIterator, UnaryOperator > > ( att ) ( ibeg, iend, obeg, oend, os.str() );
        }

        /* execute the frame and remove it */
        thread->execute(&frame);
        thread->pop();
    }
};
    
} //namespace
#endif


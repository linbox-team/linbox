#ifndef DENSE_SUB_MATRIX_H
#define DENSE_SUB_MATRIX_H
#include <linbox/util/debug.h>
#include <linbox/blackbox/dense-matrix1.h>
#include <linbox/blackbox/archetype.h>

namespace LinBox
{
  template<class Field>
    class DenseSubMatrix : public BlackboxArchetype<DenseMatrix<Field>::Vector>
  {
  public:
 
    typedef DenseMatrix<Field>::Element Element;
    typedef DenseMatrix<Field>::Vector Vector;
    typedef DenseMatrix<Field>::pointer pointer;
    class RawIterator;
    class ConstRawIterator;
    typedef DenseMatrix<Field>::RowIterator RowIterator;
    typedef DenseMatrix<Field>::ConstRowIterator ConstRowIterator;
    typedef DenseMatrix<Field>::Row Row;
    typedef DenseMatrix<Field>::ConstRow ConstRow;
    typedef DenseMatrix<Field>::ColOfRowsIterator ColOfRowsIterator;
    typedef DenseMatrix<Field>::ConstColOfRowsIterator ConstColOfRowsIterator;
    typedef DenseMatrix<Field>::ColIterator ColIterator;
    typedef DenseMatrix<Field>::ConstColIterator ConstColIterator;
    typedef DenseMatrix<Field>::Col Col;
    typedef DenseMatrix<Field>::ConstCol ConstCol;
    typedef DenseMatrix<Field>::RowOfColsIterator RowOfColsIterator;
    typedef DenseMatrix<Field>::ConstRowOfColsIterator ConstRowOfColsIterator;

    DenseSubMatrix() {}
    DenseSubMatrix(DenseMatrix<Field>* _M, size_t _beg_row, size_t _end_row, size_t _beg_col, size_t _end_col);    
    DenseSubMatrix(const DenseSubMatrix<Field>& SM,size_t _beg_row, size_t _end_row,size_t _beg_col, size_t _end_col);      
    DenseSubMatrix(const DenseSubMatrix& SM);
    DenseSubMatrix& operator=(const DenseSubMatrix<Field>& SM);
    
    BlackboxArchetype<Vector>* clone() const
      { return new DenseSubMatrix<Field>(*this);}
	
    class RawIterator;   
    class ConstRawIterator;
   
    RawIterator rawBegin();     
    RawIterator rawEnd();
        
    ConstRawIterator rawBegin() const;       
    ConstRawIterator rawEnd() const;  
    
    ColOfRowsIterator colOfRowsBegin();
    ColOfRowsIterator colOfRowsEnd();

    ConstColOfRowsIterator colOfRowsBegin() const;
    ConstColOfRowsIterator colOfRowsEnd() const;
 
    RowOfColsIterator rowOfColsBegin();
    RowOfColsIterator rowOfColsEnd();

    ConstRowOfColsIterator rowOfColsBegin() const;
    ConstRowOfColsIterator rowOfColsEnd() const;
       
    void setEntry(size_t i, size_t j, const Element& a_ij);   
    Element& getEntry(size_t i, size_t j, Element& a_ij);
   
    size_t rowdim(void) const;  
    size_t coldim(void) const;
    
    std::ostream &write(std::ostream &os = std::cout) const;

    template<class Vect1, class Vect2>
      Vect1& apply (Vect1& y, const Vect2& x) const;
    
    Vector& apply (Vector &y, const Vector &x) const
      {  return apply<Vector,Vector>(y,x); }

    template<class Iterator1, class Iterator2>
      Iterator1& apply (Iterator1 in, const Iterator2& outbegin, const Iterator2& outend) const;

    template<class Vect1>
      Vect1& applyIn(Vect1& y) const
      {
	std::vector<Element> x(y.begin(),y.end());
	apply(y,x);
	return y;
      }
    
    Vector& applyIn (Vector& y) const
      { return applyIn<Vector>(y);}

    template<class Vect1, class Vect2>
      Vect1& applyTranspose (Vect1& y, const Vect2& x) const;

    Vector& applyTranspose (Vector& y, const Vector& x) const
      { return applyTranspose<Vector,Vector>(y,x); }

    template<class Iterator1, class Iterator2>
      Iterator1& applyTranspose (Iterator1 in, const Iterator2& outbegin, const Iterator2& outend) const;
    
    template<class Vect>
      Vect& applyTransposeIn (Vect& y) const
      {
	std::vector<Element> x(y.begin(),y.end());
	applyTranspose(y,x);
	return y;
      }
    
    
    Vector& applyTransposeIn (Vector& y) const
      { return applyTransposeIn<Vector>(y);}
    
    
    const Field& field() const
     { return M->field();}

    Row operator[] (int i);
    ConstRow operator[] (int i) const;

  protected:
    DenseMatrix<Field>* M;
    size_t beg_row;
    size_t end_row;
    size_t beg_col;
    size_t end_col;
  };
}

#include <linbox/blackbox/densesubmatrix.C>
#endif

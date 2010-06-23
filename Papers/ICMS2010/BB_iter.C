template<class BB, 
class Vect = TypeChooser<BB::InVector, 
                         BB::OutVector >::Vector >
class BB_Sequential_Iterator {     
public:
  typedef typename Vect::value_type  value_type;
private:
  BB * _mat; // matrix A
  Vect v0, v, u;

  value_type _value;

public:
[...]

  void operator ++() {
      // v <- Au 
    _mat->Apply(v, u); 
      // _value <- v0^T . v
    _domain.dotproduct(_value, v0, v);  
      // exchange, u <- A^i . u0
    swap(v, u);
  }

  value_type operator *() const {
    return _value;
  }
};

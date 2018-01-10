#ifndef __LINBOX_alt_blackbox_block_container_H
#define __LINBOX_alt_blackbox_block_container_H

namespace LinBox {

template <class Field,class Blackbox,class Value>
class AltBlackboxBlockContainer {
public:

	typedef typename MatrixDomain<Field>::OwnMatrix Block;

	AltBlackboxBlockContainer () {}

	/*
	AltBlackboxBlockContainer(const Blackbox *M,
	                       const Field &F,
	                       const Block &V)
		: U_(F,V.coldim(),M->rowdim()), V_(V),
		  M_(M), F_(F), b_(V.coldim())
	{

	}
	*/

	AltBlackboxBlockContainer(const Blackbox *M,
	                       const Field &F,
	                       const Block &U,
	                       const Block &V)
		: b_(V.coldim()), F_(&F),
		  U_(U), V_(V),
		  LastOdd_(F,V.rowdim(),V.coldim()),
		  LastEven_(F,V.rowdim(),V.coldim()),
		  val_(F,V.coldim(),V.coldim()),
		  M_(M), MD_(F) {}

	class const_iterator {
	protected:
		AltBlackboxBlockContainer &c_;
	public:
		const_iterator () {}
		const_iterator(AltBlackboxBlockContainer<Field,Blackbox,Value> &c) :
			c_(c) {}

		const_iterator& operator++() {c_.incr(); return *this;}
		const Value& operator*() {return c_.getValue();}
	};

	const_iterator begin() {
		LastOdd_=V_;
		LastEven_=V_;
		isOdd_=true;
		return const_iterator(*this);
	}

	void incr() {
		Block result(*(const_cast<Field*>(F_)),b_,b_);
		if (isOdd_) {
			M_->applyLeft(LastEven_,LastOdd_);
			MD_.mul(result,U_,LastEven_);
		} else {
			M_->applyLeft(LastOdd_,LastEven_);
			MD_.mul(result,U_,LastOdd_);
		}
		isOdd_=!isOdd_;
		for (size_t i=0;i<b_;++i) {
			for (size_t j=0;j<b_;++j) {
				typename Field::Element d;
				result.getEntry(d,i,j);
				val_.setEntry(i,j,d);
			}
		}
	}

	const Value& getValue() const {
		return val_;
	}

	size_t rowdim() const {return b_;}
	size_t coldim() const {return b_;}

	const Blackbox* getBB() const {return M_;}

protected:

	size_t b_;

	const Field *F_;

	Block U_, V_, LastOdd_, LastEven_;

	Value val_;

	bool isOdd_;

	const Blackbox *M_;

	MatrixDomain<Field> MD_;
};

}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
#ifndef __LIE_MATRIX_H__
#define __LIE_MATRIX_H__

#include <vector>
namespace LinBox{ 

template <class Ring, class Operator>
class LieMatrix {

public:
	typedef Ring Field;
	typedef typename Field::Element Element;
	typedef std::vector<Operator*> Ops;

	void scalarMulIn (const Element& s) {
		Operator* fop = op. front();
		typename Operator::RawIterator raw_p;
		for (raw_p = fop -> rawBegin(); raw_p != fop -> rawEnd(); ++ raw_p)
			f. mulin (*raw_p, s);
	}

	template <class OVect>
	LieMatrix(const OVect& ops, const Field& _f) : op(ops. begin(), ops. end()), f(_f) {
		dim = ops. front() -> rowdim();
		tmp. resize(dim);
	}

	int rowdim () const {return dim;}
	int coldim () const {return dim;}

	template <class Out, class In>
	Out& apply (Out& out, const In& in) const {

		int size = op. size();
		typename std::vector<Operator*>::const_iterator op_p;
		if (size & 1) 
			op. front() -> apply (out, in);
		else
			op. front() -> apply (tmp, in);

		-- size;	
		for (op_p = op. begin() + 1; op_p != op. end(); ++ op_p) {
			if (size & 1) 
				(*op_p) -> apply (out, tmp);
			else
				(*op_p) -> apply (tmp, out);
			-- size;
		}

		return out;
	}

	template <class Out, class In>
	In& applyTranspose (Out& out, const In& in) const {

		int size = op. size();

		typename std::vector<Operator*>::const_reverse_iterator op_p;
		if (size & 1) 
			op. back() -> applyTranspose (out, in);
		else
			op. back() -> applyTranspose (tmp, in);

		-- size;	
		for (op_p = op. rbegin() + 1; op_p != op. rend(); ++ op_p) {
			if (size & 1) 
				(*op_p) -> applyTranspose (out, tmp);
			else
				(*op_p) -> applyTranspose (tmp, out);
			-- size;
		}
			
		return out;
	}

	~LieMatrix ( ) {
		typename std::vector<Operator*>::iterator op_p;
		for (op_p = op. begin(); op_p != op. end(); ++ op_p)
			delete *op_p;
	}

	const Field& field() const {
		return f;
	}

	const Ops& getOps () const {
		return op;
	}

	void write (std::ostream& out) {
		typename std::vector<Operator*>::iterator op_p;
		int i = 1;
		for (op_p = op. begin(); op_p != op. end(); ++ op_p) {
			std::cout << i <<"th Matrix:\n";
			(*op_p) -> write (out);
			++ i;
		}
	}

protected:
	std::vector<Operator*> op;
	mutable std::vector<Element> tmp;
	int dim;
	Field f;
};
}

#include <linbox/algorithms/matrix-mod.h>

namespace LinBox {

	template <class Blackbox, class Field>
	struct MatrixModTrait;

	template <class Ring, class Ops, class Field>
	struct MatrixModTrait<LieMatrix<Ring, Ops>, Field> {
		typedef LieMatrix<Field, typename MatrixModTrait<Ops, Field>::value_type> value_type;
	};

	namespace MatrixMod {
		template <class FMatrix, class IMatrix, class Field>
		void mod (FMatrix* & Ap, const IMatrix& A, Field F);

		template <class FOp, class Field, class Ring, class ROp>
		void mod (LieMatrix<Field, FOp>*& Ap, const LieMatrix<Ring, ROp>& A, Field F) {
			typedef typename LieMatrix<Ring,  ROp>::Ops ROps;
			typedef typename LieMatrix<Field, FOp>::Ops FOps;
			const ROps& rop = A. getOps();
		
			FOps fop (rop.size());
			typename FOps::iterator fop_p; 
			typename ROps::const_iterator rop_p;
			for (fop_p = fop. begin(), rop_p = rop. begin(); 
				 fop_p != fop. end(); ++ fop_p, ++ rop_p)
				mod (*fop_p, **rop_p, F);

			Ap = new LieMatrix<Field, FOp> (fop, F);
			return;
		}
	}
}

#endif

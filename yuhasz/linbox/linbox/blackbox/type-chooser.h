namespace LinBox
{

template <class T, class S, class D> class TypeChooser { 
	public:
		typedef D TYPE;
		TypeChooser() { };

};

template <class T, class D> class TypeChooser<T , T, D> { 
	public:
		typedef T TYPE;
		TypeChooser() { }
};

}


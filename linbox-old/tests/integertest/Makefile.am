TESTS=integertry

# integertry works with integer.h 

%:%.cc
	$(CXX) $(OPTFLAGS) $< -o $@ -I$(prefix)/include \
	-R /usr/local/algebra/gmp/lib -R /usr/local/gnu/gcc/lib \
	$(GMP_CFLAGS) $(GMP_LIBS) 

all: integertry

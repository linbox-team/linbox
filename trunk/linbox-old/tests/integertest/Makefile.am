TESTS=integertry

# integertry works with integer.h 

%:%.cc
	$(CXX) $(OPTFLAGS) $< -o $@ -I$(prefix)/include $(GMP_CFLAGS) $(GMP_LIBS) 

all: integertry

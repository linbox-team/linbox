CXXFLAGS+=${OPTFLAGS}

test: test.o densegfrank.o
	g++ ${CXXFLAGS} $^ -o $@

clean:
	\rm -rf test test.o densegfrank.o

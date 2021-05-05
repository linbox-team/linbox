# this Makefile is a temporary hack.  
# codegen/ is at top level so that ../examples/Makefile finds libtool ok.

# just to show gen itself is not (yet) LinBox dependent.
gen: gen.C
	clang++ gen.C -o gen

# borrow a makefile
%:%.C
	make -f ../examples/Makefile $@

run-sample: sample
	charpoly-d ../examples/data/lowmat.sms
	minpoly-s ../examples/data/lowmat.sms
	det-d-a-2 ../examples/data/lowmat.sms 
	rank-s-4097 ../examples/data/lowmat.sms 1000019
	trace-s-2 ../examples/data/lowmat.sms 3
	valence-d-a-4097 ../examples/data/lowmat.sms 
	isPositiveDefinite-d ../examples/data/lowmat.sms 
	smithForm-s ../examples/data/lowmat.sms

sample: charpoly-d minpoly-s det-d-a-2 isPositiveDefinite-d rank-s-4097 smithForm-s trace-s-2 valence-d-a-4097 
	touch sample

# build sources by code gen

#charpoly
charpoly-d.C: gen
	gen charpoly d

charpoly-s.C: gen
	gen charpoly s

charpoly-d-a-2.C: gen
	gen charpoly d a 2

charpoly-s-4097.C: gen
	gen charpoly s 4097

#minpoly
minpoly-d.C: gen
	gen minpoly d

minpoly-s.C: gen
	gen minpoly s

minpoly-d-a-2.C: gen
	gen minpoly d a 2

minpoly-s-4097.C: gen
	gen minpoly s 4097

#det
det-d.C: gen
	gen det d

det-s.C: gen
	gen det s

det-d-a-2.C: gen
	gen det d a 2

det-s-4097.C: gen
	gen det s 4097

#isPositiveDefinite
isPositiveDefinite-d.C: gen
	gen isPositiveDefinite d

isPositiveDefinite-s.C: gen
	gen isPositiveDefinite s

#rank
rank-d.C: gen
	gen rank d

rank-s.C: gen
	gen rank s

rank-d-a-2.C: gen
	gen rank d a 2

rank-s-4097.C: gen
	gen rank s 4097

#smith-form
smithForm-d.C: gen
	gen smithForm d

smithForm-s.C: gen
	gen smithForm s

#trace
trace-d.C: gen
	gen trace d

trace-s.C: gen
	gen trace s

trace-d-4097.C: gen
	gen trace d a 4097

trace-s-2.C: gen
	gen trace s 2

#valence
valence-d.C: gen
	gen valence d

valence-s.C: gen
	gen valence s

valence-d-a-4097.C: gen
	gen valence d a 4097

valence-s-2.C: gen
	gen valence s 2

run: run-rank run-det run-charpoly

run-rank: rank
	rank-d vsmallmat
	rank-s vsmallmat
	rank-d-a-2 vsmallmat
	rank-s-4097 vsmallmat

rank: rank-d rank-s rank-d-a-2 rank-s-4097
	touch rank

run-det: det
	det-d vsmallmat
	det-s vsmallmat
	det-d-a-2 vsmallmat
	det-s-4097 vsmallmat

det: det-d det-s det-d-a-2 det-s-4097
	touch det

clean:
	rm -f gen rank* det* charpoly* minpoly* isP* smith* trace* valence* 
	rm -f sample

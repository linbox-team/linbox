TESTDIRS=Random Gauss Sequence
include Makefile.rules

tests:
	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)

clean:
	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)

distclean:
	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)



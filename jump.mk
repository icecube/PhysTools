# executes the actual makefile from the build directory
# based on example by Paul D. Smith

.SUFFIXES:

OBJDIR := build

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f '$(CURDIR)'/Makefile BASEDIR='$(CURDIR)' $(MAKECMDGOALS)

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

Makefile : ;
%.mk :: ;

% :: $(OBJDIR) ; @:

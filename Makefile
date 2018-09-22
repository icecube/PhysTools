# magic to do the build from the build directory, as though it were in the source directory
ifneq ($(notdir $(CURDIR)),build)
include jump.mk
else
VPATH = $(BASEDIR)/src:$(BASEDIR)/PhysTools
include config.mk


NAME = PhysTools
STAT_PRODUCT = lib$(NAME).a
DYN_PRODUCT = lib$(NAME)$(DYN_SUFFIX)

INCLUDE = -I$(BASEDIR) -I$(INCDIR)
CXXFLAGS += $(INCLUDE) -fPIC

MAKEDEPEND = $(CXX) $(CXXFLAGS) -MM -MF $*.d $<

#libraries that users will need to link with
LIBRARIES=boost_iostreams

# core source files
SRCS = gnuplot.cpp gnuplot-normalize.cpp pipe.cpp plottable_histograms.cpp hdf5_serialization.cpp axis.cpp bin_types.cpp memory_dump.cpp tableio.cpp optimization/lbfgsb/linpack.cpp optimization/lbfgsb/lbfgsb.cpp optimization/ParameterSet.cpp
HEADERS = $(shell ls $(BASEDIR)/PhysTools/*.h $(BASEDIR)/PhysTools/*.tcpp)
DETAIL_HEADERS = $(shell ls $(BASEDIR)/PhysTools/detail/*.h)
LBFGSB_HEADERS = $(shell ls $(BASEDIR)/PhysTools/optimization/lbfgsb/*.h)
LIKELIHOOD_HEADERS = $(shell ls $(BASEDIR)/PhysTools/likelihood/*.h)
OPTIMIZATION_HEADERS = $(shell ls $(BASEDIR)/PhysTools/optimization/*.h)

OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.P)

#
# RULES
#

.PHONY: all clean install test

all: $(STAT_PRODUCT) $(DYN_PRODUCT)

$(DYN_PRODUCT) : $(OBJS)
	@echo Linking $(DYN_PRODUCT)
	@$(CXX) $(LDFLAGS) $(DYN_OPT) -shared $(OBJS) -o $(DYN_PRODUCT)
$(STAT_PRODUCT) : $(OBJS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJS)


#handle dependiecies automatically using Tom Tromey's method, as described by Paul D. Smith

%.o : %.cpp
	@echo Computing dependecies of $(<F)
	@$(MAKEDEPEND); \
            cp $*.d $*.P; \
            sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
                -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
            rm -f $*.d
	@echo Compiling $(<F)
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

  -include $(SRCS:.cpp=.P)

clean: 
	rm -f $(OBJS)
	rm -f $(DEPS)
	rm -f $(STAT_PRODUCT)
	rm -f $(DYN_PRODUCT)

install : 
	@mkdir -p $(PREFIX)
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(STAT_PRODUCT) $(PREFIX)/lib/$(STAT_PRODUCT)
	@cp $(DYN_PRODUCT) $(PREFIX)/lib/$(DYN_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/$(NAME)
	@mkdir -p $(PREFIX)/include/$(NAME)
	@cp $(HEADERS) $(PREFIX)/include/$(NAME)/
	@mkdir -p $(PREFIX)/include/$(NAME)/detail
	@cp $(DETAIL_HEADERS) $(PREFIX)/include/$(NAME)/detail/
	@mkdir -p $(PREFIX)/include/$(NAME)/likelihood
	@cp $(LIKELIHOOD_HEADERS) $(PREFIX)/include/$(NAME)/likelihood/
	@mkdir -p $(PREFIX)/include/$(NAME)/optimization
	@cp $(OPTIMIZATION_HEADERS) $(PREFIX)/include/$(NAME)/optimization/
	@mkdir -p $(PREFIX)/include/$(NAME)/optimization/lbfgsb
	@cp $(LBFGSB_HEADERS) $(PREFIX)/include/$(NAME)/optimization/lbfgsb/
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp phystools.pc $(PREFIX)/lib/pkgconfig
	@../check_install.sh phystools "$(PREFIX)"

uninstall :
	@echo Removing headers from $(PREFIX)/include/PhysTools
	@rm -rf $(PREFIX)/include/PhysTools
	@echo Removing libraries from $(PREFIX)/lib
	@rm -f $(PREFIX)/lib/$(DYN_PRODUCT)
	@rm -f $(PREFIX)/lib/$(STAT_PRODUCT)
	@echo Removing config information from $(PREFIX)/lib/pkgconfig
	@rm -f $(PREFIX)/lib/pkgconfig/phystools.pc

docs : doxygen

doxygen : ../docs/html/index.html

../docs/html/index.html : $(BASEDIR)/doxyfile $(HEADERS) $(DETAIL_HEADERS) $(LBFGSB_HEADERS) $(LIKELIHOOD_HEADERS)
	@cd .. && doxygen doxyfile
	@touch ../docs/html/index.html

test : all
	@cd ../test ; ./run_tests
check : all
	@cd ../test ; ./run_tests

# close build directory check
endif

# NO_UNIT_TESTS=1

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

EXECUTABLE = LBMflow
all: ${EXECUTABLE}

OBJ = \
	warnings.o \
	control.o \
	main.o

.PHONY: clean allclean
clean:
	-rm -f *.o *~ *.d

allclean: clean
	-rm -f ${EXECUTABLE}

CXX = g++

SRC ?= ./src
INC =

ERRFLAGS =
FEAFLAGS =
FEAFLAGS += -D_GIT_KEY=$(GITKEY) 

ifdef NO_UNIT_TESTS
	## for production: eliminate the unit tests
	FEAFLAGS += -D NO_UNIT_TESTS
endif

FEAFLAGS += -D DEVEL

OPTFLAGS = -O0
OPTFLAGS += -g -pedantic -Wall -Wno-format-security -Wno-format


EXTFLAGS ?=

LDFLAGS += -lm
##	LDFLAGS += -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

CXXFLAGS = -std=c++11 ${FLAGS}


DEP=$(OBJ:%.o=%.d)

${EXECUTABLE}: ${OBJ}
	${CXX} ${CXXFLAGS} $^ ${LDFLAGS}  -o $@

%.o:%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.d:%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}

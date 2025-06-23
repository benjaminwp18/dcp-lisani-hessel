# Copyright 2011 Jose-Luis Lisani <joseluis.lisani@uib.es>
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.


# C source code
CSRC	= library/libImageFormats.c 
# C++ source code
CXXSRC	= library/parser.cpp

# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= dehazeDCP

default	: $(BIN)

# use DEBUG mode with `make DEBUG=1`
ifdef DEBUG
  # C optimization flags
  COPT = -g
else
  # C optimization flags
  COPT	= -O3 -ftree-vectorize -funroll-loops
endif


# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) 
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra

# link flags
LDFLAGS	= -lpng -lm -ltiff -ljpeg 



# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
CXXFLAGS	+= -fopenmp
LDFLAGS += -lgomp 
endif



# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS) 
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS) 
    
# link all the object code
$(BIN): % : %.o $(OBJ) 
	$(CXX) -o $@ $^ $(LDFLAGS) 

# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
distclean	: clean
	$(RM) $(BIN) *.o *~


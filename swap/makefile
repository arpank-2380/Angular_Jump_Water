COMPILER = gfortran
OUTPUT := ./swap.x
#FLAGS = -static -ffpe-trap=invalid -fbacktrace -fdump-core 
FLAGS = -g -O3 #-Wno-maybe-uninitialized
SOURCE = io_module tools swap
EXECUTABLE = main

OBJS = source/$(SOURCE:%=%.o)
ALL_OBJS = $(SOURCE:%=%.o) $(EXECUTABLE:%=%.o)

.PHONY : all
all: $(ALL_OBJS) $(EXECUTABLE)

.PHONY : clean
clean:
	@echo "Cleaning old objective files and .mod files"
	@\rm -rf *.o *.mod  $(OUTPUT)
%.o: ./%.f90
	@echo "Compiling: " $*
	$(COMPILER) -c $(FLAGS) $< -o $@

%: %.o $(ALL_OBJS)
	@echo "   Linking program:   " $a
	$(COMPILER) $(FLAGS) $(ALL_OBJS) -o $(OUTPUT)


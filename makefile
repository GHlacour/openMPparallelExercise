# Makefile

# ### Default settings
ifeq ($(MACHINE), )
  CC = icc
  # LIBDIR = -L/cm/shared/apps/intel/compilers_and_libraries/linux/lib/intel64
  CPPFLAGS = -qopenmp
  CC_FLAGS = -B -O3 -lm -no-prec-div -qopenmp
  LAPACK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -liomp5
endif

LIBS=$(LAPACKL) $(LAPACK) $(LIBDIR)

all: echo aggregate

echo:
	@echo Compiling on the machine $(MACHINE) with $(CC)

aggregate: aggregate.o randomlib.o
	$(CC) randomlib.o aggregate.o $(LIBS) -o aggregate $(CC_FLAGS)

aggregate.o: aggregate.c aggregate.h lapack.h

randomlib.o: randomlib.h randomlib.c

clean:
	\rm *.o

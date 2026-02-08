ifeq ($(ENABLE_MPI),true)
CC   = mpiicx
DEFINES  = -D_MPI
else
CC = icx
endif

GCC  = gcc
LINKER = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -qopenmp
endif

ifeq ($(strip $(DEBUG)),true)
CFLAGS   = -O0 -g -std=c99 $(OPENMP)
else
CFLAGS   = -O3 -xHost -qopt-zmm-usage=high -std=c99 $(OPENMP) -Wno-unused-command-line-argument
endif

VERSION  = --version
LFLAGS   = $(OPENMP)
DEFINES  += -D_GNU_SOURCE
INCLUDES =
LIBS     =

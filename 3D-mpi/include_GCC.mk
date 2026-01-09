ifeq ($(ENABLE_MPI),true)
CC   = mpicc
DEFINES  = -D_MPI
else
CC = gcc
endif

GCC  = gcc
LINKER = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -fopenmp
endif

ifeq ($(strip $(DEBUG)),true)
CFLAGS   = -O0 -g -std=c99 $(OPENMP)
else
CFLAGS   = -Ofast -ffreestanding -std=c99 $(OPENMP)
endif

VERSION  = --version
LFLAGS   = $(OPENMP)
DEFINES  += -D_GNU_SOURCE
INCLUDES =
LIBS     = -lm

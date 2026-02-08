CC   = clang
GCC  = cc
LINKER = $(CC)

ifeq ($(strip $(ENABLE_OPENMP)),true)
OPENMP   = -fopenmp
#OPENMP   = -Xpreprocessor -fopenmp #required on Macos with homebrew libomp
LIBS     = # -lomp
endif
ifeq ($(strip $(DEBUG)),true)
CFLAGS   = -O0 -g -std=c17
else
CFLAGS   = -O3 -std=c17 $(OPENMP)
endif

VERSION  = --version
LFLAGS   = $(OPENMP) -lm
DEFINES  = -D_GNU_SOURCE
INCLUDES =

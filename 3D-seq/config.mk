# Supported: GCC, CLANG, ICX
TAG ?= ICX
ENABLE_OPENMP ?= false
# Supported: rb, mg
SOLVER ?= mg
# Run in debug settings
DEBUG ?= false

#Feature options
OPTIONS +=  -DARRAY_ALIGNMENT=64
OPTIONS +=  -DVERBOSE
#OPTIONS +=  -DDEBUG

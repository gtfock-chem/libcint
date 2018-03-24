alignlen = 64

# For normal cluster
# CXX = icpc
# CC  = icc
# FC  = ifort

# For Cori (use Intel compilers by default)
CXX = CC
CC  = cc -std=gnu99
FC  = ftn

AR  = xiar rcs

ERD_LIB    = /global/homes/h/huangh/scratch/OptErd_Makefile/lib
SIMINT_DIR = /global/homes/h/huangh/scratch/gtfock-simint/build-avx512/install

INC=-I. -I${SIMINT_DIR}/include

OPTFLAGS = -m64 -qno-offload

#SRC = $(wildcard *.c)
SRC = cint_basisset.c \
      erd_integral.c  \
      oed_integral.c  \
      cint_simint.c
CFLAGS = -O3 -Wall -w2 -qopenmp -g
CFLAGS += -Wunknown-pragmas -Wunused-variable
CFLAGS += ${OPTFLAGS}
CFLAGS += -D__ALIGNLEN__=${alignlen}

LIBCINT = libcint.a
OBJS := $(addsuffix .o, $(basename $(SRC)))

all: ${LIBCINT} test_cint_simint test_cint_opterd

test_cint_simint: test_cint_simint.o libcint.a
	$(CC) -qopenmp -o test_cint_simint test_cint_simint.o libcint.a -L${ERD_LIB} -lerd -loed -L${SIMINT_DIR}/lib64 -lsimint

test_cint_opterd: test_cint_opterd.o libcint.a
	$(CC) -qopenmp -o test_cint_opterd test_cint_opterd.o libcint.a -L${ERD_LIB} -lerd -loed

runtest:
	./test_cint_simint sto-3g-cart.gbs water.xyz
	./test_cint_opterd sto-3g-sph.gbs water.xyz

${LIBCINT}: ${OBJS}
	${AR} $@ $^

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@ 

clean:
	rm -f *.o *.s *.d *~ *.a

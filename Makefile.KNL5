alignlen = 64

# Options for KNL5
CC  = icc -std=gnu99
CXX = icpc
FC  = ifort

# options for Cori (uses Intel by default)
#CC  = cc -std=gnu99
#CXX = CC
#FC  = ftn

AR  = xiar rcs

INC=-I. -I/home/huangh/gtfock-simint/build-avx512/install/include

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
	$(CC) -qopenmp -o test_cint_simint test_cint_simint.o libcint.a -L/home/huangh/OptErd_Makefile/lib -lerd -loed -L/home/huangh/gtfock-simint/build-avx512/install/lib64 -lsimint -lifcore

test_cint_opterd: test_cint_opterd.o libcint.a
	$(CC) -qopenmp -o test_cint_opterd test_cint_opterd.o libcint.a -L/home/huangh/OptErd_Makefile/lib -lerd -loed -lifcore

runtest:
	./test_cint_simint sto-3g-cart.gbs water.xyz
	./test_cint_opterd sto-3g-sph.gbs water.xyz

${LIBCINT}: ${OBJS}
	${AR} $@ $^

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@ 

clean:
	rm -f *.o *.s *.d *~ *.a

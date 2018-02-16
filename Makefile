alignlen = 64

CC  = icc -std=gnu99
CXX = icpc
FC  = ifort

# options for Cori (uses Intel by default)
CC  = cc -std=gnu99
CXX = CC
FC  = ftn

AR  = xiar rcs

INC=-I.

OPTFLAGS = -m64 -qno-offload


SRC = $(wildcard *.c)
CFLAGS = -O3 -Wall -w2 -qopenmp
CFLAGS += -Wunknown-pragmas -Wunused-variable
CFLAGS += ${OPTFLAGS}
CFLAGS += -D__ALIGNLEN__=${alignlen}

LIBCINT = libcint.a
OBJS := $(addsuffix .o, $(basename $(SRC)))

all: ${LIBCINT}

${LIBCINT}: ${OBJS}
	${AR} $@ $^

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@ 

clean:
	rm -f *.o *.s *.d *~ *.a

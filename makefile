F77 = gfortran
AR = ar

.PHONY: all clean

CANONICAL_SRCS = alpha.F convergence.F factor.F gaunt.F spherical.F sym_to_rect.F
# add after fixing mesh.F

CANONICAL_OBJS += $(CANONICAL_SRCS:.F=.o)

CANONICAL_LIB = libcanonical.a

.F.o: Makefile
	$(F77) -c $<

$(CANONICAL_LIB): $(CANONICAL_OBJS)
	$(AR) cr $@ $(CANONICAL_OBJS)

all: $(CANONICAL_LIB)

clean:
	rm -f *.o $(CANONICAL_LIB)

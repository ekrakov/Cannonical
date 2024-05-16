F77 = gfortran
AR = ar

.PHONY: all clean

CANONICAL_SRCS = VECGEN.F SMTRX_STRUCTURE.F CROSS.F FACTOR.F CONV_GLF.F LIMAO.F YLMRK.F GAUNT.F SYM_TO_REC_KTRNSF.F MESH.F 


CANONICAL_OBJS += $(CANONICAL_SRCS:.F=.o)

CANONICAL_LIB = ../lib/libcanonical.a

all: $(CANONICAL_LIB) STR_mod 

.F.o: Makefile
	$(F77) -c $<

$(CANONICAL_LIB): $(CANONICAL_OBJS)
	$(AR) cr $@ $(CANONICAL_OBJS)

STR_mod: STR_mod.F
	$(F77) -o $@ $<

clean:
	rm -f *.o $(CANONICAL_LIB) STR_mod

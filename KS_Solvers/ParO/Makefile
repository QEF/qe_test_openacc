# Makefile for ParO

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../ELPA/src $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

PARO = \
paro_gamma.o paro_gamma_new.o \
paro_k.o paro_k_new.o \
pcg_gamma.o \
pcg_k.o \
bpcg_gamma.o \
bpcg_k.o \
paro_gamma_new_gpu.o  \
paro_gamma_new_acc.o  \
paro_k_new_gpu.o \
paro_k_new_acc.o \
bpcg_gamma_gpu.o \
bpcg_gamma_acc.o \
bpcg_k_gpu.o \
bpcg_k_acc.o

all : libparo.a

libparo.a:  $(PARO)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend

# Makefile for DAVID

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../ELPA/src $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

DAVID = \
cegterg.o \
regterg.o \
cegterg_gpu.o \
cegterg_acc.o \
regterg_gpu.o \
regterg_acc.o 

ifdef PGI_POWER_WORKAROUND
# Circumvents bug in the Power implementation of PGI compilers
F90FLAGS:= -O1 $(filter-out -fast,$(F90FLAGS))
endif

all : libdavid.a


libdavid.a:  $(DAVID)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend

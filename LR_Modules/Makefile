# Makefile for LR_Modules

include ../make.inc

# location of needed modules and included files (if any)
MODFLAGS=$(BASEMOD_FLAGS) $(MOD_FLAG)$(TOPDIR)/PW/src

LR_MODULES = \
adddvscf.o \
addusdbec.o \
addusdbec_nc.o \
apply_dpot.o \
apply_dpot_mod.o \
addusldos.o \
adddvepsi_us.o \
addusddenseq.o \
cft_wave.o \
cgsolve_all.o \
ccgsolve_all.o \
cg_psi.o \
ccg_psi.o \
ch_psi_all.o \
ch_psi_all_complex.o \
commutator_Hx_psi.o \
commutator_Vhubx_psi.o \
Coul_cut_2D_ph.o \
check_q_points_sym.o \
dfpt_tetra_mod.o \
dynmat_sub.o \
h_prec.o \
incdrhoscf.o \
incdrhoscf_nc.o \
dv_of_drho.o \
dgradcorr.o \
dnonloccorr.o \
dv_vdW_DF.o \
dv_rVV10.o \
newdq.o \
compute_intq.o \
set_intq_nc.o \
lr_sym_mod.o \
orthogonalize.o \
setup_offset_beta.o \
setup_nscf.o \
set_dbecsum_nc.o \
set_int3_nc.o \
set_kplusq_nc.o \
smallgq.o \
lr_sm1_psi.o \
lr_addusddens.o \
localdos.o \
setup_dgc.o \
compute_vsgga.o \
setup_dmuxc.o \
setup_nbnd_occ.o \
setup_alpha_pv.o \
set_small_group_of_q.o \
sgam_lr.o \
star_q.o \
setqmod.o \
lr_orthoUwfc.o \
qdipol_cryst.o \
mode_group.o \
mix_pot.o \
lrcom.o \
lr_dot.o \
lanczos_pseudohermitian.o \
lanczos_nonhermitian.o 

TLDEPS=mods pwlibs

all : liblrmod.a

liblrmod.a: $(LR_MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend

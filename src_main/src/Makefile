### ### Makefile for FrontFlow Red 
###
### $Date: 2004/07/12 14:14:14 $
### $Revision: 1.2 $
###
 
### include via platform 
include Makefile.in
 
 
### do not edit the following lines
MODDIR = ../mods_main_$(OS)
INC_DIR = -I../..
GARBAGE = *.o *~ *.mod $(MODDIR)/*.mod
SHELL = /bin/sh
RM = rm -f

### Fortran compiler related macros
FFLAGS = $(FOPTIONS) $(INC_DIR) $(MPI_INC)
#FLIBS = -lls
#FLIBS = -llssrv  # (if needed)

### C compiler related macros
CFLAGS = $(COPTIONS) $(INC_DIR)


### Fortran source files
FSRCS1 = module_total.f
FSRCS2 = bc_total.f \
        cal_total.f  conv_term.f \
        dc_total.f diff_total.f \
        grad_total.f metric_fflow.f rad_fvm.f \
        particle.f \
        frontflow.f  vel_admin.f \
        nml_total.f print_total.f \
        read_total.f reset_vars.f set_deltat.f \
        set_trnsprp.f solve_total.f \
        src_total.f utl_total.f hpcmw.f \
        chem_unstdy.f \
        move.f potential.f 
MAIN    = ../frontflow_main.f
OBJMAIN    = ./frontflow_main.o
#OBJMAIN = $(MAIN:.f=.o)

### C source files
###################################
####  No longer used *c   #########
###################################
#CSRCS = gfall.c gfrd_c.c gfutil_c.c gfwrt_c.c


### products
BIN_DIR = ..
INSTALL_DIR = ../../..
FFR_SINGLE = $(BIN_DIR)/fflowS_$(OS)
FFR_PARALLEL = $(BIN_DIR)/fflowHPC_$(OS)

### FFR static library
FFR_SINGLE_LIB = $(BIN_DIR)/libserial_$(OS).a
FFR_PARALLEL_LIB = $(BIN_DIR)/libhpc_$(OS).a


### Object files
USRC_DIR = ../../src_user
USRCS = $(USRC_DIR)/USER_BC.f $(USRC_DIR)/USER_chem_RATE.f\
        $(USRC_DIR)/USER_ini.f $(USRC_DIR)/USER_ini_E2P.f\
        $(USRC_DIR)/USER_ini_particle.f $(USRC_DIR)/USER_moving_mesh.f\
        $(USRC_DIR)/USER_output.f $(USRC_DIR)/USER_rad_GasPara.f\
        $(USRC_DIR)/USER_rad_pro.f $(USRC_DIR)/USER_src_T.f\
        $(USRC_DIR)/USER_src_cavitation.f $(USRC_DIR)/USER_src_e2p.f\
        $(USRC_DIR)/USER_src_fire.f $(USRC_DIR)/USER_src_r.f\
        $(USRC_DIR)/USER_src_uvw.f $(USRC_DIR)/USER_suf_chem_rate.f\
        $(USRC_DIR)/USER_property.f\
	$(USRC_DIR)/USER_PARTICLE_TEMP.f\
        $(USRC_DIR)/USER_spray_evaporation.f\
        $(USRC_DIR)/USER_particle_HeatExch.f\
        $(USRC_DIR)/USER_particle_mov_f1.f\
        $(USRC_DIR)/USER_rough_wall.f\
        $(USRC_DIR)/USER_UserDefine_region.f


UOBJS = $(USRCS:.f=.o)

COBJS = $(CSRCS:.c=.o)
FOBJS1 = $(FSRCS1:.f=.o)
FOBJS2 = $(FSRCS2:.f=.o)
SRCSINGLE = $(MAIN) module_dumyhpc.f $(FSRCS1) dumyhpc.f $(FSRCS2) $(CSRCS) $(USRCS)
OBJSINGLE = $(OBJMAIN) module_dumyhpc.o $(FOBJS1) dumyhpc.o $(FOBJS2) $(COBJS) $(UOBJS) 
OBJPARALLEL = $(OBJMAIN) module_hpc.o $(FOBJS1) hpc.o $(FOBJS2) $(COBJS) $(UOBJS)
#module_dumyhpc.o: module_dumyhpc.f
#dumypc.o: dumyhpc.f

### suffix rules
.f.o:
#	echo '//     compile $< ...'
#	echo $(F90) $(FFLAGS) -c $<
	$(F90) $(FFLAGS) -c $<
.c.o:
#	echo '//     compile $< ...'
	$(CC) $(CFLAGS) -c $<


### targets
#..SILENT:
#all: single single_lib clean_single_obj parallel parallel_lib
all: 
	make clean_parallel_obj
	make single
	make single_lib
	make clean_single_obj
	make parallel
	make parallel_lib
single: $(OBJSINGLE)
#	echo '//     link $< ...'
#	echo $(F90LINKER) $(FOBJS) $(COBJS) $(UOBJS) -o $(FFR_SINGLE) $(FLIBS)
	$(F90LINKER) $(FFLAGS) $(OBJSINGLE) -o $(FFR_SINGLE) $(FLIBS)
	cp -f $(FFR_SINGLE) $(INSTALL_DIR)/fflowS
	cp -f $(FFR_SINGLE) ../../bin_FFR/fflowS

single_lib: $(OBJSINGLE) install
	$(AR) $(FFR_SINGLE_LIB) $(OBJSINGLE)
	$(RANLIB) $(FFR_SINGLE_LIB)

parallel: $(OBJPARALLEL)
#	echo '//     link $< ...'
#	echo '//     link $< ...'
	$(MPI_F90) $(FFLAGS) $(MPI_INC) $(OBJPARALLEL) -o $(FFR_PARALLEL) $(FLIBS) $(MPI_LIBS)
	cp -f $(FFR_PARALLEL) $(INSTALL_DIR)/fflowHPC
	cp -f $(FFR_PARALLEL) ../../bin_FFR/fflowHPC

parallel_lib: $(OBJPARALLEL) install
	$(AR) $(FFR_PARALLEL_LIB) $(OBJPARALLEL)
	$(RANLIB) $(FFR_PARALLEL_LIB)

$(UOBJS): $(USRCS)
	$(F90) $(FFLAGS) -c $*.f -o $@
#	$(F90) $(FFLAGS) -c $< -o $@

$(OBJMAIN): $(MAIN)
	$(F90) $(FFLAGS) -c $? -o $@

install:
	-rm -fr $(MODDIR)
	-mkdir  $(MODDIR)
	-cp -f *.mod $(MODDIR)

clean_single_obj: 
	$(RM) $(GARBAGE) $(OBJSINGLE)

clean_parallel_obj: 
	$(RM) $(GARBAGE) $(OBJPARALLEL)

clean:
	$(RM) $(GARBAGE) $(OBJMAIN) $(UOBJS) $(FFR_SINGLE) $(FFR_SINGLE_LIB) $(FFR_PARALLEL) $(FFR_PARALLEL_LIB)

echo:
	echo '//     F90 = $(F90)'
	echo '//     FFLAGS = $(FFLAGS)'
	echo '//     F90LINKER = $(F90LINKER)'

	echo 
	echo '//     CC = $(CC)'
	echo '//     CFLAGS = $(CFLAGS)'
	echo
	echo '1', $(FFLAGS)
	echo '2', $(FOBJS)
	echo '3', $(COBJS)
	echo '4', $(UOBJS)
	echo '5', $(FLIBS)
	echo '6', -o $(FFR_SINGLE)
# eof


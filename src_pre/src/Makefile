### ### Makefile for FrontFlow Red
###
### $Date: 2004/08/17 00:25:44 $
### $Revision: 1.8 $
###

### include via platform
include Makefile.in


### do not edit the following lines
MODDIR = ../mods_pre_$(OS)
INC_DIR = -I../..
GARBAGE = *.o *~ *.mod $(MODDIR)/*.mod
SHELL = /bin/sh
RM = rm -f

### Preprocessor
#FPP=-DSTARCD -DGAMBIT -DSCRYU

### Fortran compiler related macros
FFLAGS = $(FOPTIONS) $(FPP) $(INC_DIR)
DEBUGFLG1 = :conform_check=1:div_check=3
DEBUGFLG2 = :subscript_check=1:trap_uninitialized=1
DEBUGFLG3 = :varargs_interface_check=1:verbose_runtime=ON:fullwarn=1
FDEBUG = -g -w -c -64 -O2 -static -mips4 -DEBUG$(DEBUGFLG1)$(DEBUGFLG2)$(DEBUGFLG3) $(FPP) $(INC_DIR)

#FLIBS = -lls
#FLIBS = -llssrv  # (if needed)

### C compiler related macros
CFLAGS = $(COPTIONS) $(INC_DIR)


### Fortran source files
FSRCS1 = pre_module_total.f
FSRCS2 = pre_array_prefflow.f pre_metric_total.f pre_list_total.f pre_bc_total.f\
         pre_gfutil.f pre_nml_total.f pre_print_total.f pre_read_total.f\
         pre_src_total.f pre_part.f pre_part_comm.f pre_part_cyc.f\
         pre_part_data.f pre_part_output.f pre_part_paraset.f\
         pre_part_wall.f pre_read_hpc.f pre_list_hpc.f pre_read_griddata.f\
         pre_rad_function.f pre_rad_interface.f pre_rad_mczone.f pre_rad_sxf.f
MAIN    = ../pre_frontflow.f
OBJMAIN = ./pre_frontflow.o

### C source files 
###########################
#### No longer used *c ####
###########################
#CSRCS = gfall.c gfrd_c.c gfutil_c.c gfwrt_c.c


### products
BIN_DIR = ..
INSTALL_DIR = ../../..
FFR_PRE = $(BIN_DIR)/prefflow_$(OS)

### FFR static library
FFR_PRE_LIB = $(BIN_DIR)/libprefflow_$(OS).a


### Object files
COBJS = $(CSRCS:.c=.o)
FOBJS1 = $(FSRCS1:.f=.o)
FOBJS2 = $(FSRCS2:.f=.o)
OBJPRE = $(FOBJS1) $(OBJMAIN) $(FOBJS2) $(COBJS)

### Lib files
LIBPRE = -L../../src_metis_4.1_fflow -lmetis

### suffix rules
.f.o:
	$(F90) $(FFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<


### targets
#..SILENT:
all: pre pre_lib
pre: $(OBJPRE)
	$(F90LINKER) $(FFLAGS) $(OBJPRE) -o $(FFR_PRE) $(FLIBS) $(LIBPRE)
	cp -f $(FFR_PRE) $(INSTALL_DIR)/prefflow
	cp -f $(FFR_PRE) ../../bin_FFR/prefflow

pre_lib: $(OBJPRE) install
	$(AR) $(FFR_PRE_LIB) $(OBJPRE)
	$(RANLIB) $(FFR_PRE_LIB)

debug: $(OBJPRE)
	$(F90LINKER) $(FDEBUG) $(OBJPRE) -o $(FFR_PRE) $(FLIBS) $(LIBPRE)
	cp -f $(FFR_PRE) $(INSTALL_DIR)/prefflow

$(OBJMAIN): $(MAIN)
	$(F90)  $(FFLAGS) -c $? -o $@

install:
	-rm -fr $(MODDIR)
	-mkdir  $(MODDIR)
	-cp -f *.mod $(MODDIR)

clean:
	$(RM) $(GARBAGE) $(OBJMAIN) $(FFR_PRE) $(FFR_PRE_LIB)

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
	echo '6', -o $(FFR_PRE)
# eof


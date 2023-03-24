### ### Makefile for FFR2VIZP
###
### $Date: 2006/11/02 11:11:11 $
### $Revision: 1.11 $
###

### include via platform
include Makefile.in


### do not edit the following lines
MODDIR = ../mods_vizp_$(OS)
INC_DIR = -I../.. -L../../
GARBAGE = *.o *~ *.mod $(MODDIR)/*.mod
SHELL = /bin/sh
RM = rm -f

### Preprocessor
#FPP=-DSTARCD -DGAMBIT -DSCRYU -DAVS -DFLUENT -DFVASCII -DENSIGHT
#FPP=-cpp
### Fortran compiler related macros
FFLAGS = $(FOPTIONS) $(FPP) $(INC_DIR)
#FLIBS = -lls
#FLIBS = -llssrv  # (if needed)

### C compiler related macros
CFLAGS = $(COPTIONS) $(INC_DIR)


### Fortran source files
FSRCS1 =
FSRCS2 = 
MAIN    = ffr2fvpc.f
OBJMAIN    = ffr2fvpc.o
#OBJMAIN = $(MAIN:.f=.o)

### C source files
CSRCS = 


### products
BIN_DIR = ..
INSTALL_DIR = ../../..
FFR_VIZP = $(BIN_DIR)/ffr2vizp_$(OS)

### FFR static library
FFR_VIZP_LIB = $(BIN_DIR)/libvizp_$(OS).a


### Object files
COBJS = $(CSRCS:.c=.o)
FOBJS1 = $(FSRCS1:.f=.o)
FOBJS2 = $(FSRCS2:.f=.o)
OBJVIZP  = $(OBJMAIN) $(FOBJS1) $(FOBJS2) $(COBJS)

### suffix rules
.f.o:
	$(F90) $(FFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<

### targets
#..SILENT:
all: vizp vizp_lib
vizp: $(OBJVIZP)
	$(F90LINKER) $(FFLAGS) $(OBJVIZP) -o $(FFR_VIZP) $(FLIBS)
	cp -f $(FFR_VIZP) $(INSTALL_DIR)/ffr2vizp
	cp -f $(FFR_VIZP) ../../bin_FFR/ffr2vizp
vizp_lib: $(OBJVIZP) install
	$(AR) $(FFR_VIZP_LIB) $(OBJVIZP)
	$(RANLIB) $(FFR_VIZP_LIB)

$(OBJMAIN): $(MAIN)
	$(F90) $(FFLAGS) -c $? -o $@

install:
	-rm -fr $(MODDIR)
	-mkdir  $(MODDIR)
	-cp -f *.mod $(MODDIR)

clean:
	$(RM) $(GARBAGE) $(OBJMAIN) $(FFR_VIZP)

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
	echo '6', -o $(FFR_VIZP)
# eof

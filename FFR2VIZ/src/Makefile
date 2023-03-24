### ### Makefile for FFR2VIZ
###
### $Date: 2004/09/13 02:49:01 $
### $Revision: 1.11 $
###

### include via platform
include Makefile.in


### do not edit the following lines
MODDIR = ../mods_viz_$(OS)
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
FSRCS1 = dataClass.f 
FSRCS2 = read_griddata.f read_partdata.f\
         setupGridGrobalInfo.f readAndStoreRslFFR.f\
	 writeFVHeader.f ftn_encode_header.f append_resultfv.f\
	 writeAVSHeader.f append_resultavs.f\
	 writeSCRYUHeader.f append_resultscry.f\
	 writeFLUENTHeader.f append_resultfluent.f\
	 append_resultstar.f writeENSIGHTHeader.f\
	 append_resultensight.f\
         read_fort1.f setNumVarsToZero.f addVal.f\
	 calculateYpUpAv.f tovtk.f
MAIN    = ../ffr2viz.f
OBJMAIN    = ./ffr2viz.o
#OBJMAIN = $(MAIN:.f=.o)

### C source files
CSRCS = 


### products
BIN_DIR = ..
INSTALL_DIR = ../../..
FFR_VIZ = $(BIN_DIR)/ffr2viz_$(OS)
FFR_MOVIE = $(BIN_DIR)/ffrmovie_$(OS)

### FFR static library
FFR_VIZ_LIB = $(BIN_DIR)/libviz_$(OS).a
FFR_MOVIE_LIB = $(BIN_DIR)/libmovie_$(OS).a


### Object files
COBJS = $(CSRCS:.c=.o)
FOBJS1 = $(FSRCS1:.f=.o)
FOBJS2 = $(FSRCS2:.f=.o) tovtkbin.o
OBJVIZ   = $(OBJMAIN) $(FOBJS1) $(FOBJS2) main.o $(COBJS)
OBJMOVIE = $(OBJMAIN) $(FOBJS1) $(FOBJS2) main_m.o $(COBJS)

### suffix rules
.f.o:
	$(F90) $(FFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<


### targets
#..SILENT:
all: viz viz_lib movie movie_lib
viz: $(OBJVIZ)
	$(F90LINKER) $(FFLAGS) $(OBJVIZ) -o $(FFR_VIZ) $(FLIBS)
	cp -f $(FFR_VIZ) $(INSTALL_DIR)/ffr2viz
	cp -f $(FFR_VIZ) ../../bin_FFR/ffr2viz
viz_lib: $(OBJVIZ) install
	$(AR) $(FFR_VIZ_LIB) $(OBJVIZ)
	$(RANLIB) $(FFR_VIZ_LIB)

movie: $(OBJMOVIE)
	$(F90LINKER) $(FFLAGS) $(OBJMOVIE) -o $(FFR_MOVIE) $(FLIBS)
	cp -f $(FFR_MOVIE) $(INSTALL_DIR)/ffrmovie
	cp -f $(FFR_MOVIE) ../../bin_FFR/ffrmovie

movie_lib: $(OBJMOVIE) install
	$(AR) $(FFR_MOVIE_LIB) $(OBJMOVIE)
	$(RANLIB) $(FFR_MOVIE_LIB)

$(OBJMAIN): $(MAIN)
	$(F90) $(FFLAGS) -c $? -o $@

install:
	-rm -fr $(MODDIR)
	-mkdir  $(MODDIR)
	-cp -f *.mod $(MODDIR)

clean:
	$(RM) $(GARBAGE) $(OBJMAIN) $(FFR_VIZ) $(FFR_MOVIE)

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
	echo '6', -o $(FFR_VIZ)
# eof

include ../Makefile.in.msvc

BINDIR = ..

INCLUDES = -I../Lib 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)

LIBS = -lmetis
#LIBS = -lmetis -lm

LD = $(CC) $(LDOPTIONS) -64 -L. -L..
METISLIB = ../libmetis.a

KMETISOBJS = kmetis.o io.o 

.c.o:
	$(CC) $(CFLAGS) -c $*.c

default: $(BINDIR)/kmetis$(VERNUM) 


$(BINDIR)/kmetis$(VERNUM): $(KMETISOBJS) $(METISLIB)
	$(LD)  -o $@ $(KMETISOBJS) $(LIBS)
	chmod 755 $@
clean:
	rm -f *.o

realclean:
	rm -f *.o ;\
        rm -f $(BINDIR)/kmetis$(VERNUM) 

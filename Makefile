PROGRAM  = flat2 
PROGRAM2 = adpgs1
PROGRAM3 = faopgs2
PROGRAM4 = plotp
PROGRAM5 = surfgs1

# Set CC compilation flags

CC = gcc
CFLAGS = -O    
LDLIBS = -lm

# Set FC compilation

FC = gfortran

OBJS2 = adp.o

OBJS3 = aom1.o aow1.o

OBJS5 = surf.o

HEADERS3 = incl.c inclx.c

HEADERS5 = surf.h

$(PROGRAM) : 
	$(FC) -std=legacy flat2.for -o $(PROGRAM)  

$(PROGRAM2) : $(OBJS2)
	$(CC) $(OBJS2) $(LDLIBS) -o $(PROGRAM2)

$(PROGRAM3) : $(OBJS3)
	$(CC) $(OBJS3) $(LDLIBS) -o $(PROGRAM3) 

$(OBJS3) : $(HEADERS3) 

$(PROGRAM4):
	$(FC) plotp01.f -o $(PROGRAM4)

$(PROGRAM5) : $(OBJS5)
	$(CC) $(OBJS5) $(LDLIBS) -o $(PROGRAM5) 

all: $(PROGRAM) $(PROGRAM2) $(PROGRAM3) $(PROGRAM4) $(PROGRAM5)

.PHONY: clean distclean
clean : 
	rm -f core *.o *.*% *.str *.me *_.dat *.fem *.err *.adp *.inf
distclean:
	rm -f core *.o *.*% *.str *.me *_.dat *.fem *.err *.adp *.inf
	rm -f adpgs1 flat2 faopgs2 plotp surfgs1 

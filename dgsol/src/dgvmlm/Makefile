#
# DGVMLM directory
#

L_ARCH   = $(DG_ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O
CFLAGS   =  $(OPTFLAGS) 
FFLAGS   =  $(OPTFLAGS)

AR       = ar -rcv

FILES = dgopt.o dginitial.o dgfun.o dgerr.o dgsel.o

dgvmlm: $(FILES) 
	$(AR) $(LIB_NAME) $(FILES) 
	rm -f *.o

clean:
	- rm -f *.o
	- rm -f *.a

.c.o:
	$(CC) $(CFLAGS) -c $*.c

.f.o:
	$(FC) $(FFLAGS) -c $*.f



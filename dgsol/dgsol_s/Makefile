#
# DGSOL_S directory.
#

L_ARCH   = $(DG_ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O
CFLAGS   =  $(OPTFLAGS) -D$(DG_ARCH)
FFLAGS   =  $(OPTFLAGS)

AR       = ar -rcv

# DGSOL libraries.

DGVMLM     = ../src/dgvmlm/$(LIB_NAME)
VMLM       = ../src/vmlm/$(LIB_NAME)
CSRCH      = ../src/csrch/$(LIB_NAME)
BLAS       = ../src/blas/$(LIB_NAME)
UTILS      = ../src/utils/$(LIB_NAME)
DGSOL_LIBS = $(DGVMLM) $(VMLM) $(CSRCH) $(BLAS) $(UTILS)

FILES = dgsol.o 

exec: $(FILES)
	$(FC) $(OPTFLAGS) -o dgsol $(FILES) $(DGSOL_LIBS) 

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f

clean:
	- rm -f *.o
	- rm -f *.sol
	- rm -f *.sum
	- rm -f dgsol

#
#   DGSOL directory
#

dgsol_s:
	cd dgsol_s; make

dgsol_p:
	cd dgsol_p; make 

libs: dgvmlm blas csrch utils vmlm

dgvmlm:
	cd src/dgvmlm; make 

blas:
	cd src/blas; make

csrch:
	cd src/csrch; make 

utils:
	cd src/utils; make

vmlm:
	cd src/vmlm; make

clean:
	cd dgsol_s;    make clean
	cd dgsol_p;    make clean
	cd src/dgvmlm; make clean
	cd src/blas;   make clean
	cd src/csrch;  make clean
	cd src/utils;  make clean
	cd src/vmlm;   make clean



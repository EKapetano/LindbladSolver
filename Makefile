#Compiler and Linker - Options
#Used Compiler:
FC = gfortran
#Compiler options:
FCOPTS = -std=f2003 -pedantic -Wall -O3
#Used Linker:
LN = $(FC)
#Linker options:
LNOPTS = 
#Paths to the folders containing the libraries blas, lapack and lapack95:
PATHTOBLAS = ../../../lib/BLAS/
PATHTOLAPACK = ../../../lib/lapack-3.5.0/
PATHTOL95 = ../../../lib/LAPACK95/
PATHTOLIBS = /usr/lib/
#Path to the folder containing the lapack95 modules:
PATHTOL95MODULES = ../../../lib/LAPACK95/lapack95_modules/


OBJS = accuracy.o inout.o calc.o debug.o lindblad.o
OBJS2 = accuracy.o inout.o calc.o res.o

#Normal Program
lindblad: $(OBJS)
	$(FC) $(FCOPTS) -fopenmp -o lindblad $(OBJS) -L$(PATHTOLIBS) -lexpokit -llapack95 -llapack -lblas


lindblad.o: lindblad.f90
	$(FC) $(FCOPTS) -fopenmp -c lindblad.f90 -I$(PATHTOL95MODULES) 

accuracy.o:
	$(FC) $(FCOPTS) -c accuracy.f90

inout.o: inout.f90 accuracy.o
	$(FC) $(FCOPTS) -fopenmp -c inout.f90

calc.o: calc.f90 accuracy.o
	$(FC) $(FCOPTS) -fopenmp -c calc.f90 -L$(PATHTOLIBS) -I$(PATHTOL95MODULES) -llapack95 -llapack -lblas

debug.o: debug.f90 accuracy.o calc.o
	$(FC) $(FCOPTS) -c debug.f90

#Program which runs the code for multiple dephasing strengths
res: $(OBJS2)
	$(FC) $(FCOPTS) -fopenmp -o res $(OBJS2) -L$(PATHTOLIBS) -I$(PATHTOL95MODULES) -lexpokit -llapack95 -llapack -lblas

res.o: accuracy.o res.f90
	$(FC) $(FCOPTS) -fopenmp -c res.f90 -L$(PATHTOLIBS) -I$(PATHTOL95MODULES) -llapack95 -llapack -lblas


#Cleaning-Options
.PHONY: clean realclean

clean:
	rm -f *.mod *.o *.dat

realclean: clean
	rm -f lindblad res

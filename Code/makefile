
CC := g++
FC := gfortran
NC := nvcc

bdir := bin

all : compile clean

compile :
	$(FC) TDFDAWFS2DSG.f90 -o $(bdir)/TDFDAWFS2DSG-F90
	$(FC) TDFDEWFS2DSG.f90 -o $(bdir)/TDFDEWFS2DSG-F90
	$(CC) TDFDEWFS2DSG.cpp -o $(bdir)/TDFDEWFS2DSG-CPP
	$(NC) TDFDEWFS2DSG.cu -o $(bdir)/TDFDEWFS2DSG-CU

clean :
	-rm -f *.mod

clear :
	-rm -f $(bdir)/*

# vim:ft=make:noet

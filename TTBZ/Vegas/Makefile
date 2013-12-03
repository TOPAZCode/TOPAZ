comp = ifort -O2 -D_WriteTmpHisto=1


vegas: vegas.o
	$(comp) -c vegas.f
	mpicc -D_WriteTmpHisto=1 -c -lpthread  pvegas_mpi.c


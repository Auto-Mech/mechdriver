FC = x86_64-conda_cos6-linux-gnu-gfortran
FFLAGS = -O3 -fPIC
FOR_SOURCES = pot_aux.f dummy_corr.f mol_corr.f
OBJECTS = $(FOR_SOURCES:.f=.o)
libcorrpot.so : ${OBJECTS}
	${FC} -shared -o $@ ${OBJECTS}
clean:
	rm -f *.o

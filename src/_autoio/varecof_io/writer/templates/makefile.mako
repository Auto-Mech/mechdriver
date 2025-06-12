FC = ${fc}
FFLAGS = -O3 -fPIC

FOR_SOURCES = pot_aux.f dummy_corr.f ${corr_potential_names} 

<%text>
OBJECTS = $(FOR_SOURCES:.f=.o)

libcorrpot.so : ${OBJECTS} 
	${FC} -shared -o $@ ${OBJECTS}  
</%text>

clean:
	rm -f *.o

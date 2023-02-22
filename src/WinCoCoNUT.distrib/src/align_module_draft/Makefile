#LD=/usr/local/bin/g++
#GNULIB=/usr/local/lib
#GNULIB64=/usr/local/lib/sparcv9
# Linux
#LD=g++ -O3 -Wall
#GNULIB=/usr/lib
#GNULIB64=/usr/bin/lib/sparcv9

#LIBGNU=-L$(GNULIB) -R$(GNULIB)

include ../Makedef

OFILES=alichainer.o fragment.o align_chain.o memmap.o contigs.o alignment.o 
CFILES=alichainer.cpp fragment.cpp align_chain.cpp memmap.cpp contigs.cpp alignment.cpp 
EXEFILE=alichainer 

$(EXEFILE):	$(OFILES)
	$(LD)  $(CFLAGS) $(LIBGNU) $(OFILES) -o $(EXEFILE)
alichainer.o:	alichainer.cpp fragment.hpp align_chain.hpp memmap.hpp contigs.hpp alignment.cpp alignment.hpp
	$(LD)  $(CFLAGS) $(LIBGNU) -c alichainer.cpp
alignment.o:	alignment.cpp alignment.hpp
	$(LD)  $(CFLAGS) $(LIBGNU) -c alignment.cpp
fragment.o:	fragment.cpp fragment.hpp 
	$(LD)  $(CFLAGS) $(LIBGNU) -c fragment.cpp
align_chain.o:	align_chain.cpp align_chain.hpp fragment.hpp memmap.hpp
	$(LD)  $(CFLAGS) $(LIBGNU) -c align_chain.cpp
memmap.o: memmap.cpp memmap.hpp
	$(LD)  $(CFLAGS) $(LIBGNU) -c memmap.cpp
contigs.o: contigs.cpp contigs.hpp
	$(LD)  $(CFLAGS) $(LIBGNU) -c contigs.cpp
	$(LD)  $(CFLAGS) $(LIBGNU) invert_matches.cpp -o invert_matches

64bits:
	$(LD) $(LIBGNU) $(COPT) $(CFILES) -o $(EXEFILE)	

.PHONY:	clean_all clean
clean_all:
	rm -f $(OFILES) $(EXEFILE)
clean:
	rm -f $(OFILES) $(EXEFILE)


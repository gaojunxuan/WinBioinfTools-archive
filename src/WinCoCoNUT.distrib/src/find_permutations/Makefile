#LD=/usr/local/bin/g++

# for debug
#LD=g++ -Wall -O3

#LD=g++
#GNULIB=/usr/local/lib
#GNULIB64=/usr/local/lib/sparcv9
#LIBGNU=-L$(GNULIB) -R$(GNULIB)

include ../Makedef

OFILES=main.o   Interval.o    Point.o  Transformer.o kdTree.o ContigSorting.o
CFILES=main.cpp   Interval.cpp    Point.cpp  Transformer.cpp kdTree.cpp ContigSorting.cpp
EXEFILE=chainer2permutation.x

$(EXEFILE):	$(OFILES)
	$(LD) $(CFLAGS) $(LIBGNU) $(OFILES) -o $(EXEFILE)
main.o:	main.cpp Transformer.cpp Interval.cpp  Point.cpp  Transformer_fixed.cpp   kdTree.cpp Interval.h  Point.h  Transformer.h kdTree.h ContigSorting.cpp
	$(LD) $(CFLAGS) $(LIBGNU) -c main.cpp
Interval.o:	Interval.cpp Interval.h
	$(LD) $(CFLAGS) $(LIBGNU) -c Interval.cpp 
Point.o:	Point.cpp Point.h
	$(LD) $(CFLAGS) $(LIBGNU) -c Point.cpp
kdTree.o:	kdTree.cpp kdTree.h
	$(LD) $(CFLAGS) $(LIBGNU) -c kdTree.cpp
Transformer.o:	Transformer.cpp Transformer.h 
	$(LD) $(CFLAGS) $(LIBGNU) -c Transformer.cpp
ContigSorting.o: ContigSorting.cpp ContigSorting.h
	$(LD) $(CFLAGS) $(LIBGNU) -c ContigSorting.cpp

64bits:
	$(LD) $(LIBGNU) $(COPT) $(CFILES) -o $(EXEFILE)	

.PHONY:	clean_all clean
clean_all:
	rm -f $(OFILES) $(EXEFILE) 
clean:
	rm -f $(OFILES) $(EXEFILE) 

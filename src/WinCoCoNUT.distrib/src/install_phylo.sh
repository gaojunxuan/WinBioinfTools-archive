export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mohamed/Myprojects/CoCoNUT/bin/phylobin
cd revDist
make clean
make
cd ../../bin/phylobin
ln -s ../../src/revDist/librevDist.so.1.0 librevDist.so
cd ../../src/revDist
make -f Makefile.start
mv revDist ../../bin/phylobin
cd ..

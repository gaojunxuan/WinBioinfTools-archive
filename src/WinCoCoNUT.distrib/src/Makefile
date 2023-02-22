all:
	$(MAKE) -C align_cdna_mod
	$(MAKE) -C align_module_draft
	$(MAKE) -C cdna_postprocessing
	$(MAKE) -C chainer 
	$(MAKE) -C filter_alignment
	$(MAKE) -C find_permutations
	$(MAKE) -C clustalw1.83

64bit:
	$(MAKE) -C align_cdna_mod 64bits
	$(MAKE) -C align_module_draft 64bits
	$(MAKE) -C cdna_postprocessing 64bits
	$(MAKE) -C chainer 64bits
	$(MAKE) -C filter_alignment 64bits
	$(MAKE) -C find_permutations 64bits
	$(MAKE) -C clustalw1.83

clean:
	$(MAKE) -C align_cdna_mod clean
	$(MAKE) -C align_module_draft clean
	$(MAKE) -C cdna_postprocessing clean
	$(MAKE) -C chainer clean
	$(MAKE) -C filter_alignment clean
	$(MAKE) -C find_permutations clean
	$(MAKE) -C clustalw1.83 clean

install:
	mkdir -p ../bin/align
	mkdir -p ../bin/chainer
	cp align_cdna_mod/estchainer../bin/align/
	cp align_cdna_mod/invert_matches_cdna ../bin/align/
	cp clustalw1.83/clustalw ../bin/align/
	cp align_module_draft/alichainer ../bin/align/
	cp align_module_draft/invert_matches ../bin/align/
	cp cdna_postprocessing/estchainer2cluster ../bin/chainer/
	cp chainer/chainer ../bin/chainer/
	cp chainer/match2chainer ../bin/chainer/
	cp chainer/rep2chainer  ../bin/chainer/
	cp chainer/multimat2chainer ../bin/chainer/
	cp chainer/ramaco2chainer ../bin/chainer/
	cp chainer/rep2chainer4memspe. ../bin/chainer/
	cp filter_alignment/filteralign.x ../bin/align/
	cp find_permutations/chainer2permutation.x  ../bin/chainer/
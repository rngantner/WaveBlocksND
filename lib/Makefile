all:
	$(MAKE) python -C Arnoldi_PyCpp/src
	$(MAKE) python -C Pade_PyCpp/src
	cp Arnoldi_PyCpp/src/carnoldi.so ../src/WaveBlocksND/carnoldi.so
	cp Pade_PyCpp/src/cpade.so ../src/WaveBlocksND/cpade.so

clean:
	$(MAKE) clean -C Arnoldi_PyCpp/src
	$(MAKE) clean -C Pade_PyCpp/src
	$(RM) ../src/WaveBlockND/*.so

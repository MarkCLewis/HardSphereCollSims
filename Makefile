MPICC	= /home/mlewis/openmpi-1.2.6/build/bin/mpiCC
CCFLAGS	= -O -openmp

%: %.cpp
	$(MPICC) -o $@ $(CCFLAGS) $<

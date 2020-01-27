FLAGS=-Wall -pedantic -std=c++14 -Ofast -fopenmp -I/home/mlewis/workspace/RingsVersion3_OMP/
#FLAGS=-Ofast -LNO -march=auto -mso -mp -I/users/mlewis/workspace/RingsVersion3_OMP/

RingsSim: RingsSim.cpp
	g++ $(FLAGS) -o RingsSim RingsSim.cpp

MPIRingsSim: MPIRingsSim.cpp
	mpiCC $(FLAGS) -o MPIRingsSim MPIRingsSim.cpp
#	/users/mlewis/openmpi/bin/mpiCC $(FLAGS) -o MPIRingsSim MPIRingsSim.cpp

RecoverSim: RecoverSim.cpp *.h
	g++ $(FLAGS) -Wall -pedantic -o RecoverSim RecoverSim.cpp

MPIRecoverSim: MPIRecoverSim.cpp ../*.h
	mpiCC $(FLAGS) -o MPIRecoverSim MPIRecoverSim.cpp

BasicMain: BasicMain.cpp *.h
	g++ $(FLAGS) -Wall -pedantic -o BasicMain BasicMain.cpp

TestMain: TestMain.cpp *.h
	g++ $(FLAGS) -Wall -pedantic -o TestMain TestMain.cpp

SpinTestMain: SpinTestMain.cpp *.h
	g++ $(FLAGS) -Wall -pedantic -o SpinTestMain SpinTestMain.cpp

GracenRingsSim: GracenRingsSim.cpp
	g++ $(FLAGS) -o GracenRingsSim GracenRingsSim.cpp

MPIGracenRingsSim: MPIGracenRingsSim.cpp
	g++ $(FLAGS) -o MPIGracenRingsSim MPIGracenRingsSim.cpp

clean:

all:
	make RingsSim
	make BasicMain

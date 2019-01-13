// ProcessorCommunication.h
// This is a class that is used to wrap most of the MPI functionality.  One
// of these should be declared in the program and then references of it should
// be stored in the classes that are going to need to communicate with other
// processes.  This mainly includes the boundary conditions for moving particles
// around and then later on gravity routines.

#ifdef PARALLEL

#ifndef PROCESSOR_COMMUNICATION
#define PROCESSOR_COMMUNICATION

#include <mpi.h>

class ProcessorCommunication;

ProcessorCommunication *firstMade=NULL;

class ProcessorCommunication {
	public:
		ProcessorCommunication(int argc,char **argv) {
			if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
				fprintf(stderr, "MPI initialization error\n");
				exit(EXIT_FAILURE);
			}
			MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
			MPI_Comm_rank(MPI_COMM_WORLD, &thisProc);
			mytag=1;

			// Set up neighbors.
			neighbors[0]=-1;
			neighbors[1]=-1;
			neighbors[2]=-1;
			neighbors[3]=(thisProc+numProcs-1)%numProcs;
			neighbors[4]=-1;
			neighbors[5]=(thisProc+1)%numProcs;
			neighbors[6]=-1;
			neighbors[7]=-1;
			neighbors[8]=-1;

			if(firstMade==NULL) {
				firstMade=this;
			}
		}

		~ProcessorCommunication() {
			MPI_Finalize();
		}

		int getProcessNum() {
			return thisProc;
		}

		int getNumProcesses() {
			return numProcs;
		}

		int getNeighbor(int neighborNum) {
			return neighbors[neighborNum];
		}

		// Returns error code.
		int sendTo(int num,std::vector<double> &buffer) {
			if (MPI_Send(&(buffer[0]),buffer.size(),MPI_DOUBLE,num,mytag,MPI_COMM_WORLD) == MPI_SUCCESS) {
				return 0;
			} else {
				return -1;
			}
		}

		int sendToNeighbor(int neighborNum,std::vector<double> &buffer) {
			if (MPI_Send(&(buffer[0]),buffer.size(),MPI_DOUBLE,neighbors[neighborNum],mytag,MPI_COMM_WORLD) == MPI_SUCCESS) {
				return 0;
			} else {
				return -1;
			}
		}

		// Returns error code.
		int isendTo(int num,std::vector<double> &buffer,MPI_Request *req) {
			if (MPI_Isend(&(buffer[0]),buffer.size(),MPI_DOUBLE,num,mytag,MPI_COMM_WORLD,req) == MPI_SUCCESS) {
				return 0;
			} else {
				return -1;
			}
		}

		int isendToNeighbor(int neighborNum,std::vector<double> &buffer,MPI_Request *req) {
			if (MPI_Isend(&(buffer[0]),buffer.size(),MPI_DOUBLE,neighbors[neighborNum],mytag,MPI_COMM_WORLD,req) == MPI_SUCCESS) {
				return 0;
			} else {
				return -1;
			}
		}

		// Returns error code.
		int broadcast() {
			// !!! Not written yet.
			return 0;
		}

		// Returns error code or node read from.
		int readFrom(std::vector<double> &buffer,int &numRead) {
			MPI_Status status;
			if(MPI_Recv(&(buffer[0]),buffer.size(),MPI_DOUBLE,MPI_ANY_SOURCE,mytag,MPI_COMM_WORLD,&status) == MPI_SUCCESS) {
//				MPI_Get_count(&status,MPI_INT,&numRead);
				return status.MPI_SOURCE;
			} else {
				return -1;
			}
		}

		// Returns error code.
		int readFrom(int num,std::vector<double> &buffer,int &numRead) {
			MPI_Status status;
			if(MPI_Recv(&(buffer[0]),buffer.size(),MPI_DOUBLE,num,mytag,MPI_COMM_WORLD,&status) == MPI_SUCCESS) {
//				MPI_Get_count(&status,MPI_INT,&numRead);
				return 0;
			} else {
				return -1;
			}
		}

		int readFromNeighbor(int neighborNum,std::vector<double> &buffer,int &numRead) {
			MPI_Status status;
			if(MPI_Recv(&(buffer[0]),buffer.size(),MPI_DOUBLE,neighbors[neighborNum],mytag,MPI_COMM_WORLD,&status) == MPI_SUCCESS) {
//				MPI_Get_count(&status,MPI_INT,&numRead);
				return 0;
			} else {
				return -1;
			}
		}

		int ireadFromNeighbor(int neighborNum,std::vector<double> &buffer,int &numRead,MPI_Request *req) {
			if(MPI_Irecv(&(buffer[0]),buffer.size(),MPI_DOUBLE,neighbors[neighborNum],mytag,MPI_COMM_WORLD,req) == MPI_SUCCESS) {
				return 0;
			} else {
				return -1;
			}
		}

	private:
		int numProcs;
		int thisProc;

		// Process numbers of adjacent cells.
		int neighbors[9];

		int mytag;
};

void printProcessNumber() {
	printf("Process %d\n",firstMade->getProcessNum());
}

#endif

#endif

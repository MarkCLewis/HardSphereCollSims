// ArrayOfLists.h
// This class is intended to help with the finding of collisions.  The lists
// should be doubly linked lists so that just having an index into a list is
// sufficient for deleting.  This will be important becuase we want to be able
// to "randomly" delete items.  When two particles undergo a collision, all
// future potential collisions they might have been involved in need to be
// removed.

// For the first element in the list, the prev will be the negative of
// which list it is in minus 1.
//
// This class has been synchronized to support OpenMP parallelization.

#ifndef ARRAY_LIST
#define ARRAY_LIST

#include<vector>
#include<stdio.h>

template<class Data>
class ArrayOfLists {
	public:
		ArrayOfLists() {
			firstFree=-1;
			maxSize=-1;
			omp_init_lock(&lock);
			wait.resize(omp_get_max_threads());
			for (int i = 0; i < (int)wait.size(); i++)
			{
				omp_init_lock(&wait[i]);
				omp_set_lock(&wait[i]);
			}
			waiting.resize(omp_get_max_threads());
			for (int i = 0; i < (int)waiting.size(); i++)
				waiting[i]=false;
		}

		~ArrayOfLists() {
			omp_destroy_lock(&lock);
			for (int i = 0; i < (int)wait.size(); i++)
			{
				omp_unset_lock(&wait[i]);
				omp_destroy_lock(&wait[i]);
			}
		}

		void init(double maxValue) {
			if(maxSize<100) {
				head.resize(100000);
			} else {
				head.resize(maxSize/2);
			}
			spacing=head.size()/maxValue;
			for(unsigned int i=0; i<head.size(); ++i) head[i]=-1;
			pool.resize(0);
			firstFree=-1;
			pullingFrom=0;
			maxSize=0;
			numActive=0;
		}

		// Returns the node in the pool it was added at for quick deletes
		// later on.
		int add(Data &d) {
			omp_set_lock(&lock);
			int newNode;
			if(firstFree>-1) {
				newNode=firstFree;
				firstFree=pool[firstFree].next;
			} else {
				newNode=pool.size();
				pool.resize(pool.size()+1);
			}
			pool[newNode].data=d;
			int binNum=(int)(d.getValue()*spacing);
			if(binNum<0 || binNum>=(int)head.size()) {
				printf("Tried to add element to bin %d.  Max is %lu\n",binNum,head.size());
				fflush(stdout);
				binNum=head.size()-1;
			}
			if(binNum<pullingFrom) {
				//printf("Adding behind pulling cell.\n");
				//fflush(stdout);
				pullingFrom=binNum;
			}
			int rover=head[binNum],prev=-1;
			while(rover!=-1 && pool[rover].data.getValue()<d.getValue()) {
				prev=rover;
				rover=pool[rover].next;
			}
			if(prev==-1) {
				pool[newNode].next=head[binNum];
				pool[newNode].prev=-binNum-1;
				if(head[binNum]>-1) pool[head[binNum]].prev=newNode;
				head[binNum]=newNode;
			} else {
				pool[newNode].next=rover;
				pool[newNode].prev=prev;
				if(rover>-1) pool[rover].prev=newNode;
				pool[prev].next=newNode;
			}
			if((int)pool.size()>maxSize) maxSize=pool.size();
			omp_unset_lock(&lock);
			return newNode;
		}

		// Get the next element returns false if there wasn't one.
		template <class GridType>
		bool getNext(Data &d, GridType &grid) {
			omp_set_lock(&lock);
			while(pullingFrom<(int)head.size() && head[pullingFrom]<0) {
				++pullingFrom;
			}
			if(pullingFrom>=(int)head.size())
			{
				omp_unset_lock(&lock);
				return false;
			}
			int bin=pullingFrom;
			int ret=head[bin];
			d=pool[ret].data;

			while (!checkIfSafe(d,grid))
			{
			   if (pool[ret].next>-1)
			   {
				ret=pool[ret].next;
				d=pool[ret].data;
			   }
			   else
			   {
				do bin++;
				while (bin<(int)head.size() && head[bin]==-1);

				if (bin>=(int)head.size())
				{
				   if (numActive<=0)
				   {
				   	omp_unset_lock(&lock);
				   	return false;
				   }
				   else
				   {
					//Set this thread's wait value to true, so that another
					//thread will release its lock in finishCollision().
					//We also have to unset the synchronization lock to make
					//the finishCollision() method accessible.
					waiting[omp_get_thread_num()]=true;
					omp_unset_lock(&lock);
					//This will block until this thread's lock is released.
					//This is behaviorally similar to wait() in Java.
					omp_set_lock(&wait[omp_get_thread_num()]);
				   	return getNext(d,grid);
				   }
				}
				else
				{
				   ret=head[bin];
				   d=pool[ret].data;
				}
			   }
			}
			numActive++;
			omp_unset_lock(&lock);
			return true;
		}

		// Returns false if that was not a valid node to delete.
		bool removeAt(int which) {
			omp_set_lock(&lock);
			if(which<0 || which>=(int)pool.size() || pool[which].prev==-2000000000)
			{
				omp_unset_lock(&lock);
			   	return false;
			}
			if(pool[which].prev<0) {
				// Head of list so deal with that.
				head[-pool[which].prev-1]=pool[which].next;
			} else {
				// Not the head of a list.
				pool[pool[which].prev].next=pool[which].next;
			}
			if(pool[which].next>-1) {
				pool[pool[which].next].prev=pool[which].prev;
			}
			pool[which].prev=-2000000000;	// This line should allow me to
				// check quickly if an element is already on the free list.
			pool[which].next=firstFree;
			firstFree=which;
			omp_unset_lock(&lock);
			return true;
		}

		Data &getData(int which) {
			omp_set_lock(&lock);
			if(which<0 || which>=(int)pool.size()) {
				printf("Warning: attempt to get invalid data from pool in ArrayOfLists.\n");
				fflush(stdout);
				omp_unset_lock(&lock);
				return pool[0].data;
			}
			omp_unset_lock(&lock);
			return pool[which].data;
		}

		template <class GridType>
		void finishCollision(Data &d, GridType &grid) {
			if(d.p2<0) return;
			omp_set_lock(&lock);
			grid.setInUse(d.p1,false);
			if (d.p1!=d.p2)
				grid.setInUse(d.p2,false);
			numActive--;

			//Release the locks for all threads that are currently waiting.
			//This is behaviorally similar to notifyAll() in Java.
			for (int i = 0; i < (int)waiting.size(); i++)
			{
				if (waiting[i])
				{
					waiting[i]=false;
					omp_unset_lock(&wait[i]);
				}
			}
			omp_unset_lock(&lock);
		}

	private:
		template <class GridType>
		bool checkIfSafe(Data &d, GridType &grid) {
			if(d.p2<0) return true;
			if (grid.isSafe(d.p1) && (d.p1==d.p2 || grid.isSafe(d.p2)))
			{
				grid.setInUse(d.p1,true);
				if (d.p1!=d.p2)
					grid.setInUse(d.p2,true);
				return true;
			}
			return false;
		}

		class DLLNode {
			public:
				int prev,next;
				Data data;
		};

		std::vector<int> head;
		std::vector<DLLNode> pool;
		int firstFree;
		int maxSize;
		double spacing;
		int pullingFrom;
		omp_lock_t lock;
		std::vector<bool> waiting;
		std::vector<omp_lock_t> wait;
		int numActive;
};

#endif

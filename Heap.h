// Heap.h
// This is a heap based data structure for storing the collisions that
// particles will undergo if they maintain their current paths.

// I wrote this thinking that I would use it with the PotentialCollisions
// but it doesn't work well because I have to delete all elements later
// on that involve one of those particles.  That is worst case O(n log n)
// and average case of O(n) so it doesn't really help me much.  I think that
// having multiple lists will be better.

template <class Data>
class Heap {
	public:
		void add(Data &d) {
			if(vect.size()<2) {
				vect.resize(2);
				vect[1]=d;
			}
			int i=vect.size();
			vect.resize(vect.size()+1);
			while(i>0 && d<vect[i/2]) {
				vect[i]=vect[i/2];
				i/=2;
			}
			vect[i]=d;
		}

		Data &getFirst() {
			return vect[0];
		}

		bool removeFirst() {
			if(vect.size()<=0) return false;
			percolateDown(0);
			return true;
		}

		bool removeAt(int which) {
			if(vect.size()<which) return false;
			percolateDown(which);
			return true;
		}
	private:
		// This version of percolate down will also overwrite what had been
		// at location which.  It takes the last thing in the heap and decides
		// where to put it.  It finishes with a resize.
		void percolateDown(int which) {
			int i=which,end=vect.size()-1;
			boolean flag=true;
			while(i<end && flag) {
				if(vect[i*2]<vect[i*2+1]) {
					if(vect[i*2]<vect[end]) {
						vect[i]=vect[i*2];
						i*=2;
					} else {
						flag=false;
					}
				} else {
					if(vect[i*2+1]<vect[end]) {
						vect[i]=vect[i*2+1];
						i=i*2+1;
					} else {
						flag=false;
					}
				}
			}
			vect[i]=vect[end];
			vect.resize(end);
		}

		vector<Data> vect;
}

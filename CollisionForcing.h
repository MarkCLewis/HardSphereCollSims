// CollisionForcing.h
// This is the class to process forcing on particles by collisions.  The
// applyForce function takes a population and goes through and calculates all
// of the collisions that occur on that population during the timestep,
// adjusting the velocities appropriately.  The beginning of the applyForce
// call builds the hash for the particles.  I'm thinking of using a "smart
// hash" here that does variables sizes so that it can be more efficient at
// handling the ringlets when the form.  It also has built into it the logic
// for keeping track of the collisions in some type of queue structure.  I think
// that here I might want to write both a PriorityQueue and an array of linked
// lists.  Perhaps maybe even an array of heaps.  If I write it in a nice
// generic way then it shouldn't matter too much which one I do.  This is also
// something that doesn't have to be shipped over MPI so I can use vector.

#ifndef COLL_FORCING
#define COLL_FORCING

#ifdef PARALLEL
extern void printProcessNumber();
#endif

#include <vector>
#include "ArrayOfLists.h"
#include "AdhesionForces.h"
#include "CollisionEvents.h"
#include "BinIndex.h"

class PotentialCollision {
	public:
		PotentialCollision():p1{-1},p2{-1},time(-1.0),nextP1(-1),nextP2(-1),prevP1(-1),prevP2(-1) {
		}

		ParticleIndex p1,p2;
		double time;
		int nextP1,nextP2;		// This assumes that I can reference
			// the PotentialCollision objects by number indexes.  This does
			// set certain limitations, but I think it is better than pointers
			// in this case.
		int prevP1,prevP2;		// These need to be doubly linked so that
			// I can remove elements at random.
		bool operator< (PotentialCollision &pc) {
			return time<pc.time;
		}
		double getValue() { return time; }
};

template<class HashStructure,class AdhesionForce=NoAdhesionForce,class OtherCollisionEvents=NoEvents>
class CollisionForcing {
	public:
		CollisionForcing() {
			omp_init_lock(&lock);
			omp_init_lock(&addLock);
		}

		~CollisionForcing() {
			omp_destroy_lock(&lock);
			omp_destroy_lock(&addLock);
		}

		template<class Population>
		void applyForce(Population &pop) {
			// Init all the variables that are needed.
			printf("Collision Forcing - init %d\n",pop.getNumBodies());
			fflush(stdout);
			if(pop.getNumBodies()<=1) return;
			if(firstPotential.size()<(unsigned int)pop.getNumBodies()) {
				firstPotential.resize(pop.getNumBodies());
			}
			for(unsigned int i=0; i<firstPotential.size(); ++i)
				firstPotential[i]=-1;
			collisions.init(pop.getTimeStep());

			// Go through and do the binning by particle positions.
			printf("Build hash\n");
			fflush(stdout);
			collisionHash.build(pop);

			af.doAdhesionForce(pop,*this);

#ifndef COLLISIONLESS
			// Now find the collisions that would result given the original
			// trajectories during the timestep and add them to the "list".
#ifdef PARALLEL
			printProcessNumber();
#endif
			printf("Find Initial\n");
			fflush(stdout);

			findInitialCollisions(pop);
			if(OtherCollisionEvents::hasEvents) {
				otherEvents.findEventsForAll(pop,*this);
			}

			// Now run through the collisions.
#ifdef PARALLEL
			printProcessNumber();
#endif
			printf("Do processing\n");
			fflush(stdout);
			//int collisionCount=0;

			#pragma omp parallel
			{
				PotentialCollision pc;
				while(collisions.getNext(pc,collisionHash))
				{
					if(OtherCollisionEvents::hasEvents && pc.p2<0) {
						otherEvents.handleEvent(pop,pc.p1,pc.p2.i,pc.time);
					} else {
						pop.processCollision(pc.p1,pc.p2,pc.time);
					}
					omp_set_lock(&lock);
					// Remove all collisions between p1 and p2.  I can do this
					// sloppy in theory and walk the list removing all elements
					// that are on it then just set the first collisions for those
					// particles to -1;
					int poolNum=firstPotential[pc.p1.i];
					while(poolNum>-1) {
						collisions.removeAt(poolNum);
						PotentialCollision &pcRef=collisions.getData(poolNum);
						removeFromParticleLists(pcRef);
						if(pcRef.p1==pc.p1) {
							poolNum=pcRef.nextP1;
						} else {
							poolNum=pcRef.nextP2;
						}
					}	
					firstPotential[pc.p1.i]=-1;
					if(!OtherCollisionEvents::hasEvents || pc.p2>=0) {
						poolNum=firstPotential[pc.p2.i];
						while(poolNum>-1) {
							collisions.removeAt(poolNum);
							PotentialCollision &pcRef=collisions.getData(poolNum);
							removeFromParticleLists(pcRef);
							if(pcRef.p2==pc.p2) {
								poolNum=pcRef.nextP2;
							} else {
								poolNum=pcRef.nextP1;
							}
						}
						firstPotential[pc.p2.i]=-1;
					}
					omp_unset_lock(&lock);
					if(OtherCollisionEvents::hasEvents && pc.p2<0) {
						otherEvents.findEventsForOne(pop,pc.p1,*this);
					} else {
						findLaterCollisions(pop,pc);
					}
					collisions.finishCollision(pc,collisionHash);
					//collisionCount++;
				}
			}


#ifdef PARALLEL
			printProcessNumber();
#endif
			//printf("%d collisions\n",collisionCount);
			//fflush(stdout);
#endif
		}

		double getMinSpacing() { return collisionHash.getMinSpacing(); }

		AdhesionForce &getAdhesionForce() { return af; }

		template<class Population,class ActionType>
		void runThroughCollisionPairs(Population &pop,ActionType &action) {
			int maxX = collisionHash.getMaxX();
			#pragma omp parallel for schedule(dynamic)
			for(int ii=0; ii<maxX; ++ii) {
				BinIndex bi = {ii};
				for(int jj=0; jj<collisionHash.getMaxY(); ++jj) {
					BinIndex bj = {jj};
					for(ParticleIndex p1=collisionHash.getFirst(bi,bj); p1.i>-1; p1=collisionHash.getNext(p1)) {
						int offsetx[5]={0,1,-1,0,1};
						int offsety[5]={0,0,1,1,1};
						for(int k=0; k<5; ++k) {
							BinIndex tx=bi+offsetx[k],ty=bj+offsety[k];
							if(tx.valid(collisionHash.getMaxX()) && ty.valid(collisionHash.getMaxY())) {
								ParticleIndex p2;
								if(k==0) p2=collisionHash.getNext(p1);
								else p2=collisionHash.getFirst(tx,ty);
								while(p2>-1) {
									action(pop,p1,p2);
									p2=collisionHash.getNext(p2);
								}
							}
						}
					}
				}
			}
		}

		template<class Population>
		void operator()(Population &pop,ParticleIndex p1,ParticleIndex p2) {
			double collTime=pop.collisionTime(p1,p2);
			if(collTime>=0.0 && collTime<pop.getTimeStep())
				addPotentialWithLinks(p1,p2,collTime);
		}

		void addPotentialWithLinks(ParticleIndex p1,ParticleIndex p2,double time,bool later=false) {
			omp_set_lock(&lock);
			PotentialCollision pc;
			pc.p1=p1;
			pc.p2=p2;
			pc.time=time;
			pc.nextP1=firstPotential[p1.i];
			pc.prevP1=-1;
			if(p2>=0) {
				pc.nextP2=firstPotential[p2.i];
				pc.prevP2=-1;
			}
			int poolNum=collisions.add(pc);
			firstPotential[p1.i]=poolNum;
			if(pc.nextP1>-1) {
				PotentialCollision &pcRef=collisions.getData(pc.nextP1);
				if(pcRef.p1==p1) {
					pcRef.prevP1=poolNum;
				} else {
					pcRef.prevP2=poolNum;
				}
			}
			if(p2>=0) {
				firstPotential[p2.i]=poolNum;
				if(pc.nextP2>-1) {
					PotentialCollision &pcRef=collisions.getData(pc.nextP2);
					if(pcRef.p2==p2) {
						pcRef.prevP2=poolNum;
					} else {
						pcRef.prevP1=poolNum;
					}
				}
			}
			omp_unset_lock(&lock);
		}

	private:
		template<class Population>
		void findInitialCollisions(Population &pop) {
			runThroughCollisionPairs(pop,*this);
		}

		template<class Population>
		void findLaterCollisions(Population &pop,PotentialCollision &pc) {
			ParticleIndex p1=pc.p1;
			BinIndex i=collisionHash.getBinX(p1);
			BinIndex j=collisionHash.getBinY(p1);
			int offsetx[9]={-1,0,1,-1,0,1,-1,0,1};
			int offsety[9]={-1,-1,-1,0,0,0,1,1,1};
			for(int k=0; k<9; k++) {
				BinIndex tx=i+offsetx[k],ty=j+offsety[k];
				if(tx.valid(collisionHash.getMaxX()) && ty.valid(collisionHash.getMaxY())) {
					ParticleIndex p2=collisionHash.getFirst(tx,ty);
					while(p2>-1) {
						if(p2!=p1 && p2!=pc.p2) {
							double collTime=pop.collisionTime(p1,p2);
							if(collTime>=0.0 && collTime<pop.getTimeStep())
								addPotentialWithLinks(p1,p2,collTime,true);
						}
						p2=collisionHash.getNext(p2);
					}
				}
			}

			p1=pc.p2;
			i=collisionHash.getBinX(p1);
			j=collisionHash.getBinY(p1);
			for(int k=0; k<9; k++) {
				BinIndex tx=i+offsetx[k],ty=j+offsety[k];
				if(tx.valid(collisionHash.getMaxX()) && ty.valid(collisionHash.getMaxY())) {
					ParticleIndex p2;
					p2=collisionHash.getFirst(tx,ty);
					while(p2>-1) {
						if(p2!=p1 && p2!=pc.p1) {
							double collTime=pop.collisionTime(p1,p2);
							if(collTime>=0.0 && collTime<pop.getTimeStep())
								addPotentialWithLinks(p1,p2,collTime,true);
						}
						p2=collisionHash.getNext(p2);
					}
				}
			}
		}

		void removeFromParticleLists(PotentialCollision &d) {
			// Adjust the lists for the first particle.
			if(d.prevP1>-1) {
				PotentialCollision &prev=collisions.getData(d.prevP1);
				if(prev.p1==d.p1) {
					prev.nextP1=d.nextP1;
				} else {
					prev.nextP2=d.nextP1;
				}
			} else {
				firstPotential[d.p1.i]=d.nextP1;
			}
			if(d.nextP1>-1) {
				PotentialCollision &next=collisions.getData(d.nextP1);
				if(next.p1==d.p1) {
					next.prevP1=d.prevP1;
				} else {
					next.prevP2=d.prevP1;
				}
			}

			// Adjust the lists for the second particle.
			if(d.p2>=0) {
				if(d.prevP2>-1) {
					PotentialCollision &prev=collisions.getData(d.prevP2);
					if(prev.p2==d.p2) {
						prev.nextP2=d.nextP2;
					} else {
						prev.nextP1=d.nextP2;
					}
				} else {
					firstPotential[d.p2.i]=d.nextP2;
				}
				if(d.nextP2>-1) {
					PotentialCollision &next=collisions.getData(d.nextP2);
					if(next.p2==d.p2) {
						next.prevP2=d.prevP2;
					} else {
						next.prevP1=d.prevP2;
					}
				}
			}
		}

		ArrayOfLists<PotentialCollision> collisions;
		std::vector<int> firstPotential;
		HashStructure collisionHash;
		AdhesionForce af;
		OtherCollisionEvents otherEvents;
		omp_lock_t lock;
		omp_lock_t addLock;
};

#endif

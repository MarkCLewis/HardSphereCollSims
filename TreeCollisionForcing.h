// TreeCollisionForcing.h
// This is the class to process forcing on particles by collisions found on
// a tree.

#include <vector>
#include "ArrayOfLists.h"
#include "AdhesionForces.h"
#include "CollisionEvents.h"

class PotentialCollision {
	public:
		PotentialCollision():p1{-1},p2{-1},time(-1.0),nextP1(-1),nextP2(-1),prevP1(-1),prevP2(-1) {
		}

		ParticleIndex p1,p2; // If p2 is negative it means we have another event type.
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

template<class TreeType,class AdhesionForce=NoAdhesionForce,class OtherCollisionEvents=NoEvents>
class TreeCollisionForcing {
	public:
		TreeCollisionForcing(TreeType &t):tree(t),otherEvents(-1) {
			omp_init_lock(&lock);
		}

		~TreeCollisionForcing() {
			omp_destroy_lock(&lock);
		}

		template<class Population>
		void applyForce(Population &pop) {
			// Init all the variables that are needed.
			printf("Collision Forcing - init %d\n",pop.getNumBodies());
			const unsigned int nb=pop.getNumBodies();
			if(nb<=1) return;
			firstPotential.resize(nb);
			for(unsigned int i=0; i<nb; ++i)
				firstPotential[i]=-1;
			collisions.init(pop.getTimeStep());

			// Go through and do the binning by particle positions.
//			printf("Build tree\n");
			tree.build(pop);

			af.doAdhesionForce(pop,*this);

			// Now find the collisions that would result given the original
			// trajectories during the timestep and add them to the "list".
			printf("Find Initial\n");
			curTime=0.0;
			tree.findAllCollisions(pop,*this);
			if(OtherCollisionEvents::hasEvents) {
				otherEvents.findEventsForAll(pop,*this);
			}

			// Now run through the collisions.
			printf("Do processing\n");
//			int collisionCount=0;

			#pragma omp parallel
			{
				PotentialCollision pc;
				while(collisions.getNext(pc,tree)) {
					if(OtherCollisionEvents::hasEvents && pc.p2<0) {
						otherEvents.handleEvent(pop,pc.p1,pc.p2.i,pc.time);
					} else {
						pop.processCollision(pc.p1,pc.p2,pc.time);
					}
					// Remove all collisions between p1 and p2.  I can do this
					// sloppy in theory and walk the list removing all elements
					// that are on it then just set the first collisions for those
					// particles to -1;
					omp_set_lock(&lock);

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

					curTime=pc.time;
					omp_unset_lock(&lock);
					
					if(OtherCollisionEvents::hasEvents && pc.p2<0) {
						otherEvents.findEventsForOne(pop,pc.p1,*this);
					} else {
						tree.findCollisions(pop,*this,pc.p1,pc.p2);
						tree.findCollisions(pop,*this,pc.p2,pc.p1);
					}

					collisions.finishCollision(pc,tree);
//					collisionCount++;
				}
			}
//			printf("%d collisions\n",collisionCount);
		}

		double getMinSpacing() { return tree.getMinSpacing(); }

		AdhesionForce &getAdhesionForce() { return af; }

		void addPotentialCollision(ParticleIndex p1,ParticleIndex p2,double t) {
			addPotentialWithLinks(p1,p2,t);
			//if(t>=curTime) addPotentialWithLinks(p1,p2,t);
			//else printf("Rejected potential %d %d %e\n",p1,p2,t);
		}

		template<class Population,class ActionType>
		void runThroughCollisionPairs(Population &pop,ActionType &action) {
			tree.runThroughCollisionPairs(pop,action);
		}

		template<class Population>
		void operator()(Population &pop,ParticleIndex p1,ParticleIndex p2) {
			double t=pop.collisionTime(p1,p2);
			if(t>=0.0 && t<pop.getTimeStep()) {
				addPotentialCollision(p1,p2,t);
			}
		}

		void addPotentialWithLinks(ParticleIndex p1,ParticleIndex p2,double time) {
			PotentialCollision pc;
			pc.p1=p1;
			pc.p2=p2;
			pc.time=time;
			omp_set_lock(&lock);
			pc.nextP1=firstPotential[p1.i];
			pc.prevP1=-1;
			if(!OtherCollisionEvents::hasEvents || p2>=0) {
				pc.nextP2=firstPotential[p2.i];
				pc.prevP2=-1;
			}
			const int poolNum=collisions.add(pc);
			firstPotential[p1.i]=poolNum;
			if(pc.nextP1>-1) {
				PotentialCollision &pcRef=collisions.getData(pc.nextP1);
				if(pcRef.p1==p1) {
					pcRef.prevP1=poolNum;
				} else {
					pcRef.prevP2=poolNum;
				}
			}
			if(!OtherCollisionEvents::hasEvents || p2>=0) {
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
		// I'm not using this function right now, but I'm keeping it around
		// just in case.  The reason I'm not using it is that I can do a
		// sloppy cleanup on that lists.
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
			if(!OtherCollisionEvents::hasEvents || d.p2>=0) {
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
		AdhesionForce af;
		TreeType &tree;
		OtherCollisionEvents otherEvents;
		double curTime;
		omp_lock_t lock;
};

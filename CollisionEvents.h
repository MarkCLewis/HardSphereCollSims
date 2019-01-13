/*
 * This holds the basic other types of collision events and handlers.  This can be
 * used to implement bouncing boundaries or the like.
 */

#ifndef COLLISION_EVENTS_H
#define COLLISION_EVENTS_H

class NoEvents {
	public:
		static const int hasEvents=false;

		NoEvents(int t=-1) {}

		template<class Population,class CollisionForce>
		void findEventsForAll(Population &pop,CollisionForce &collForce) {
		}

		template<class Population,class CollisionForce>
		void findEventsForOne(Population &pop,ParticleIndex pi,CollisionForce &collForce) {
		}

		template<class Population>
		void handleEvent(Population &pop,ParticleIndex pi,int type,double time) {
		}
	private:
};

template<class Event1,class Event2>
class DoubleEvents {
	public:
		static const int hasEvents=true;

		DoubleEvents(int type):typeNum(type),e1(type),e2(type-1) {}

		template<class Population,class CollisionForce>
		void findEventsForAll(Population &pop,CollisionForce &collForce) {
			e1.findEventsForAll(pop,collForce);
			e2.findEventsForAll(pop,collForce);
		}

		template<class Population,class CollisionForce>
		void findEventsForOne(Population &pop,int p,CollisionForce &collForce) {
			e1.findEventsForOne(pop,p,collForce);
			e2.findEventsForOne(pop,p,collForce);
		}

		template<class Population>
		void handleEvent(Population &pop,ParticleIndex pi,int type,double time) {
			if(type==typeNum) e1.handleEvent(pop,pi,type,time);
			else e2.handleEvent(pop,pi,type,time);
		}
	private:
		int typeNum;
		Event1 e1;
		Event2 e2;
};

template<class EpsPerp,class EpsPar>
class ZeroPlaneEvents {
	public:
		static const int hasEvents=true;

		ZeroPlaneEvents(int type):typeNum(type) {}

		template<class Population,class CollisionForce>
		void findEventsForAll(Population &pop,CollisionForce &collForce) {
			int nb=pop.getNumBodies();
			double step=pop.getTimeStep();
			#pragma omp parallel for
			for(int i=0; i<nb; ++i) {
				ParticleIndex pi = {i};
				if(pop.getz(pi)<pop.getRadius(pi)) {
					collForce.addPotentialWithLinks(pi,typeNum,0.0);
				} else {
					double vz=pop.getvz(pi);
					if(vz<0) {
#ifndef CONST_ACCEL
						double time=-(pop.getz(pi)-pop.getRadius(pi))/vz;
#else
						double time=2*step;
						double root=vz*vz-2.0*CONST_ACCEL*(pop.getz(pi)-pop.getRadius(pi));
						if(root>=0) {
							time=(-vz-sqrt(root))/CONST_ACCEL;
						}
#endif
						if(time<step) {
							collForce.addPotentialWithLinks(i,typeNum,time);
						}
					}
				}
			}
		}

		template<class Population,class CollisionForce>
		void findEventsForOne(Population &pop,ParticleIndex pi,CollisionForce &collForce) {
			double minStep=std::min(pop.getRadius(pi)*0.001/fabs(pop.getvz(pi)),pop.getTimeStep()*0.5);
			if(pop.getz(pi)<pop.getRadius(pi)) {
				if(pop.getTime(pi)+minStep<pop.getTimeStep())
					collForce.addPotentialWithLinks(pi,typeNum,pop.getTime(pi)+minStep);
			} else {
				double vz=pop.getvz(pi);
				if(vz<0) {
#ifndef CONST_ACCEL
					double time=-(pop.getz(pi)-pop.getRadius(pi))/vz;
#else
					double time=2*step;
					double root=vz*vz-2.0*CONST_ACCEL*(pop.getz(pi)-pop.getRadius(pi));
					if(root>=0) {
						time=(-vz-sqrt(root))/CONST_ACCEL;
					}
#endif
					if(time<step) {
						if(time<minStep) time=minStep;
						collForce.addPotentialWithLinks(pi,typeNum,time+pop.getTime(pi));
					}
				}
			}
		}

		template<class Population>
		void handleEvent(Population &pop,ParticleIndex pi,int type,double time) {
			pop.advanceParticleTo(pi,time);
			if(pop.getz(pi)<pop.getRadius(pi)) {
				pop.setz(pi,pop.getRadius(pi)*1.0000001);
			}
			if(pop.getvz(pi)<0) {
				double vPar=sqrt(pop.getvx(pi)*pop.getvx(pi)+pop.getvy(pi)*pop.getvy(pi));
				pop.setvx(pi,pop.getvx(pi)*EpsPar::epsilon(vPar));
				pop.setvy(pi,pop.getvy(pi)*EpsPar::epsilon(vPar));
				pop.setvz(pi,-pop.getvz(pi)*EpsPerp::epsilon(pop.getvz(pi)));
			}
		}
	private:
		int typeNum;
};

template<class EpsPerp,class EpsPar>
class CylinderEvents {
	public:
		static const int hasEvents=true;

		CylinderEvents(int type):typeNum(type),radius(1.0),height(1.0) {
#ifdef CYLINDER_RADIUS
			radius=CYLINDER_RADIUS;
#endif
#ifdef CYLINDER_HEIGHT
			height=CYLINDER_HEIGHT;
#endif
		}

		template<class Population,class CollisionForce>
		void findEventsForAll(Population &pop,CollisionForce &collForce) {
			int nb=pop.getNumBodies();
			double step=pop.getTimeStep();
			#pragma omp parallel for
			for(int ii=0; ii<nb; ++ii) {
				ParticleIndex pi = {ii};
				double x=pop.getx(pi);
				double y=pop.gety(pi);
				double z=pop.getz(pi);
				if(z-pop.getRadius(pi)<height) {
					double dist=sqrt(x*x+y*y);
					if(dist-pop.getRadius(pi)<radius) {
						if(dist+pop.getRadius(pi)>radius) {
							collForce.addPotentialWithLinks(pi,typeNum,0.0);
						} else {
							double vx=pop.getvx(pi);
							double vy=pop.getvy(pi);
							double vz=pop.getvz(pi);
							double a=vx*vx+vy*vy;
							double b=x*vx+y*vy;
							double c=x*x+y*y-(radius-pop.getRadius(pi))*(radius-pop.getRadius(pi));
							if(a==0.0) {
								if(b!=0.0) {
									double time=-c/b;
									if(time>=0 && time<step) {
										collForce.addPotentialWithLinks(pi,typeNum,time);
									}
								}
							} else {
								double root=b*b-4*a*c;
								if(root>=0) {
									double sroot=sqrt(root);
									double time=(-b-sroot)/(2*a);
									if(time<0) time=(-b+sroot)/(2*a);
									if(time>=0 && time<step) {
										collForce.addPotentialWithLinks(pi,typeNum,time);
									}
								}
							}
						}
					}
				}
			}
		}

		template<class Population,class CollisionForce>
		void findEventsForOne(Population &pop,ParticleIndex pi,CollisionForce &collForce) {
			double vx=pop.getvx(pi);
			double vy=pop.getvy(pi);
			double vz=pop.getvz(pi);
			if(vx==0.0 && vy==0.0 && vz==0.0) {
				return;
			}
			double minStep=std::min(pop.getRadius(pi)*0.001/(fabs(vx)+fabs(vy)+fabs(vz)),pop.getTimeStep()*0.5);
			double x=pop.getx(pi);
			double y=pop.gety(pi);
			double z=pop.getz(pi);
			if(z-pop.getRadius(pi)<height) {
				double dist=sqrt(x*x+y*y);
				if(dist-pop.getRadius(pi)<radius) {
					double step=pop.getTimeStep();
					if(dist+pop.getRadius(pi)>radius) {
						if(pop.getTime(pi)+minStep<pop.getTimeStep()) {
							collForce.addPotentialWithLinks(pi,typeNum,pop.getTime(pi)+minStep);
						}
					} else {
						double vx=pop.getvx(pi);
						double vy=pop.getvy(pi);
						double vz=pop.getvz(pi);
						double a=vx*vx+vy*vy;
						double b=2.0*(x*vx+y*vy);
						double c=x*x+y*y-(radius-pop.getRadius(pi))*(radius-pop.getRadius(pi));
						if(a==0.0) {
							if(b!=0.0) {
								double time=-c/b;
								if(time>=pop.getTime(pi) && time<step) {
									if(time<minStep) time=minStep;
									collForce.addPotentialWithLinks(pi,typeNum,time+pop.getTime(pi));
								}
							}
						} else {
							double root=b*b-4*a*c;
							if(root>=0) {
								double sroot=sqrt(root);
								double time=(-b-sroot)/(2*a);
								if(time<0) time=(-b+sroot)/(2*a);
								if(time>=pop.getTime(pi) && time<step) {
									if(time<minStep) time=minStep;
									collForce.addPotentialWithLinks(pi,typeNum,time+pop.getTime(pi));
								}
							}
						}
					}
				}
			}
		}

		template<class Population>
		void handleEvent(Population &pop,ParticleIndex pi,int type,double time) {
			pop.advanceParticleTo(pi,time);
			double x=pop.getx(pi);
			double y=pop.gety(pi);
			double z=pop.getz(pi);
			double dist=sqrt(x*x+y*y);
			if(radius-dist<pop.getRadius(pi)) {
				double frac=(radius-pop.getRadius(pi))*0.99999999/dist;
				pop.setx(pi,x*frac);
				pop.sety(pi,y*frac);
				x=pop.getx(pi);
				y=pop.gety(pi);
			}
			double vx=pop.getvx(pi);
			double vy=pop.getvy(pi);
			double vz=pop.getvz(pi);
			if(vx*x+vy*y>0) {
				double rx=x/dist;
				double ry=y/dist;
				double vn=vx*rx+vy*ry;
				double vnx=rx*vn;
				double vny=ry*vn;
				double vtx=vx-vnx;
				double vty=vy-vny;
				double vtz=vz;
				double vt=sqrt(vtx*vtx+vty*vty+vtz*vtz);
				pop.setvx(pi,vtx*EpsPar::epsilon(vt)-vnx*EpsPerp::epsilon(vn));
				pop.setvy(pi,vty*EpsPar::epsilon(vt)-vny*EpsPerp::epsilon(vn));
				pop.setvz(pi,vtz*EpsPar::epsilon(vt));
			}
		}
	private:
		int typeNum;
		double radius;
		double height;
};
#endif

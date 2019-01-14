// MoonForcing.h
// This class is used to provide the forcing of a moon's gravity as it pulls
// on the particles when they have a close approach.  It has code to allow
// the moon to be on an elliptical orbit.  If an eccentricity greater than 0
// is used, then the boundary conditions that are applied must be able to
// handle significant gradients in the forced eccentricity and the output
// methods need to dump things that make sense.
//
// It can also handle an inclined moon in which case similar restrictions
// apply for handling inclination gradients.  The BC that is probably best
// suited is one that spans a full orbital period of the moon in synodic
// space so that the particles at opposite edges both went by the moon at the
// same part of its orbit.

#ifndef MOON_FORCING
#define MOON_FORCING

class MoonForcing {
	public:
		MoonForcing(double m):mass(m),ecc(0.0),phi(0.0),inc(0.0),zeta(0.0),X(0.0),Y(0.0),num(0) {}

		MoonForcing(double m,double e_i,double phi_i):mass(m),ecc(e_i),phi(phi_i),inc(0.0),zeta(0.0),X(0.0),Y(0.0),num(0) {}

		MoonForcing(double m,double e_i,double phi_i,double i_i,double zeta_i):mass(m),ecc(e_i),phi(phi_i),inc(i_i),zeta(zeta_i),X(0.0),Y(0.0),num(0) {}

		MoonForcing(int n,double X_i,double Y_i,double m):mass(m),ecc(0.0),phi(0.0),inc(0.0),zeta(0.0),X(X_i),Y(Y_i),num(n) {}

		MoonForcing(int n,double X_i,double Y_i,double m,double e_i,double phi_i):mass(m),ecc(e_i),phi(phi_i),inc(0.0),zeta(0.0),X(X_i),Y(Y_i),num(n) {}

		MoonForcing(int n,double X_i,double Y_i,double m,double e_i,double phi_i,double i_i,double zeta_i):mass(m),ecc(e_i),phi(phi_i),inc(i_i),zeta(zeta_i),X(X_i),Y(Y_i),num(n) {}

		template <class Population>
		void applyForce(Population &pop) {
			if(X==0.0 && Y==0.0) {
				doOriginMoon(pop);
			} else {
				printf("Moon X=%e, Y=%e\n",X,Y);
				doExternalMoon(pop);
			}
		}
		
		template<class Population>
		void doOriginMoon(Population &pop) {
			if(pop.getNumBodies()<1) return;
			if(ecc==0.0 && inc==0.0) {
				#pragma omp parallel for schedule(static)
				for(int ii=0; ii<pop.getNumBodies(); ii++) {
					ParticleIndex pi = {ii};
					double sepx=-pop.getx(pi);
					double sepy=-pop.gety(pi);
					double dx=(1+sepx)*cos(sepy)-1.0;
					double dy=(1+sepx)*sin(sepy);
//					double dx=cos(pop.gety(i))-1.0-pop.getx(i);
//					double dy=-sin(pop.gety(i));
					double dist=sqrt(dx*dx+dy*dy+pop.getz(pi)*pop.getz(pi));
					double mag=pop.getTimeStep()*mass/(dist*dist*dist);
					pop.setvx(pi,pop.getvx(pi)+dx*mag);
					pop.setvy(pi,pop.getvy(pi)+dy*mag);
					pop.setvz(pi,pop.getvz(pi)+(pop.getvz(pi)-pop.getz(pi))*mag);
					pop.adjustAfterForce(pi);
				}
			} else {
				double px=-ecc*cos(phi);
				double py=GCCoords::BETA*ecc*sin(phi);
				double pz=inc*cos(zeta);
				#pragma omp parallel for schedule(static)
				for(int ii=0; ii<pop.getNumBodies(); ii++) {
					ParticleIndex pi = {ii};
					double sepx=px-pop.getx(pi);
					double sepy=py-pop.gety(pi);
					double dx=(1+sepx)*cos(sepy)-1.0;
					double dy=(1+sepx)*sin(sepy);
//					double dx=cosy*(1+px)+siny*py-(1+pop.getx(i));
//					double dy=cosy*py-siny*(1+px);
					double dz=pz-pop.getz(pi);
					double dist=sqrt(dx*dx+dy*dy+dz*dz);
					double mag=pop.getTimeStep()*mass/(dist*dist*dist);
					pop.setvx(pi,pop.getvx(pi)+dx*mag);
					pop.setvy(pi,pop.getvy(pi)+dy*mag);
					pop.setvz(pi,pop.getvz(pi)+dz*mag);
					pop.adjustAfterForce(pi);
				}
				phi+=pop.getTimeStep();
				zeta+=pop.getTimeStep();
			}
		}

		template<class Population>
		void doExternalMoon(Population &pop) {
			if(pop.getNumBodies()<1) return;
//			double esum=0.0;
			double px=X-ecc*cos(phi);
			double pz=inc*cos(zeta);
			#pragma omp parallel for schedule(static)
			for(int ii=0; ii<pop.getNumBodies(); ii++) {
				ParticleIndex pi = {ii};
//				double sumy=0.0;
				for(int n=0; n<num; ++n) {
					double py=(Y+n*2*3.14159/num)+GCCoords::BETA*ecc*sin(phi);
					double sepx=px-pop.getx(pi);
					double sepy=py-pop.gety(pi);
					while(sepy<-3.1415927) sepy+=2*3.1415927;
					while(sepy>3.1415927) sepy-=2*3.1415927;
					double dx=(1+sepx)*cos(sepy)-1.0;
					double dy=(1+sepx)*sin(sepy);
//					double dx=cosy*(1+px)+siny*sepy-(1+pop.getx(i));
//					double dy=cosy*sepy-siny*(1+px);
					double dz=pz-pop.getz(pi);
					double dist=sqrt(dx*dx+dy*dy+dz*dz);
					double mag=pop.getTimeStep()*mass/(dist*dist*dist);
//					if(i==0) printf("Sepy=%e %d of %d %e %e %e %e\n",sepy,n,num,dx,dy,dz,mag);
					pop.setvx(pi,pop.getvx(pi)+dx*mag);
					pop.setvy(pi,pop.getvy(pi)+dy*mag);
					pop.setvz(pi,pop.getvz(pi)+dz*mag);
//					sumy+=dy*mag;
				}
				pop.adjustAfterForce(pi);
//				esum+=pop.gete(pi);
//				if(i==0) printf("sumy=%e\n",sumy);
			}
			phi+=pop.getTimeStep();
			zeta+=pop.getTimeStep();
			Y-=2.0*GCCoords::A0*X*pop.getTimeStep();
//			printf("Average e is %e\n",esum/pop.getNumBodies());
		}

		bool doWrapAround(double miny,double maxy) {
			double size=2*(maxy-miny);
			for(int n=0; n<num; ++n) {
				double py=(Y+n*2*3.14159/num);
				while(py-miny<-3.1415927) py+=2*3.1415927;
				while(py-miny>3.1415927) py-=2*3.1415927;
				if(py>miny-1.0*size && py<maxy+1.0*size) return false;
			}
			return true;
		}

	private:
		double mass;
		double ecc,phi;
		double inc,zeta;
		double X,Y;
		int num;
};

template<class GCType>
class ModulatedExternalMoonForcing {
	public:
		ModulatedExternalMoonForcing(double X_i,double Y_i,double m,double len):mass(m),gc(X_i,Y_i,0.0,0.0,0.0,0.0),L(len) {}

		ModulatedExternalMoonForcing(double X_i,double Y_i,double m,double e_i,double phi_i,double len):mass(m),gc(X_i,Y_i,e_i,phi_i,0.0,0.0),L(len) {}

		ModulatedExternalMoonForcing(double X_i,double Y_i,double m,double e_i,double phi_i,double i_i,double zeta_i,double len):mass(m),gc(X_i,Y_i,e_i,phi_i,i_i,zeta_i),L(len) {}

		template <class Population>
		void applyForce(Population &pop) {
			if(pop.getNumBodies()<1) return;
			if(gc.Y<-L || gc.Y>L) return;
//			double esum=0.0;
			double moonMag=0.5*(1+cos(gc.Y*3.14159/L));
			#pragma omp parallel for schedule(static)
			for(int ii=0; ii<pop.getNumBodies(); ii++) {
				ParticleIndex pi = {ii};
//				double sumy=0.0;
				for(int n=-2; n<3; ++n) {
					GCType mgc(gc);
					mgc.advance(n*pop.getTimeStep());
					CartCoords cc(mgc);
					double px=cc.x;
					double py=cc.y;
					double pz=cc.z;
					double sepx=px-pop.getx(pi);
					double sepy=py-pop.gety(pi);
					while(sepy<-3.1415927) sepy+=2*3.1415927;
					while(sepy>3.1415927) sepy-=2*3.1415927;
					double dx=(1+sepx)*cos(sepy)-1.0;
					double dy=(1+sepx)*sin(sepy);
//					double dx=cosy*(1+px)+siny*sepy-(1+pop.getx(i));
//					double dy=cosy*sepy-siny*(1+px);
					double dz=pz-pop.getz(pi);
					double dist=sqrt(dx*dx+dy*dy+dz*dz);
					double mag=moonMag*pop.getTimeStep()*mass/(dist*dist*dist);
//					if(i==0) printf("Sepy=%e %d of %d %e %e %e %e\n",sepy,n,num,dx,dy,dz,mag);
					pop.setvx(pi,pop.getvx(pi)+dx*mag);
					pop.setvy(pi,pop.getvy(pi)+dy*mag);
					pop.setvz(pi,pop.getvz(pi)+dz*mag);
//					sumy+=dy*mag;
				}
				pop.adjustAfterForce(pi);
//				esum+=pop.gete(pi);
//				if(i==0) printf("sumy=%e\n",sumy);
			}
			gc.advance(pop.getTimeStep());
			while(gc.Y<-3.14159) gc.Y+=2*3.14159;
			while(gc.Y>3.14159) gc.Y-=2*3.14159;
//			printf("Average e is %e\n",esum/pop.getNumBodies());
		}

		bool doWrapAround(double miny,double maxy) {
			return true;
		}

	private:
		double mass;
		GCType gc;
		double L;
};

#endif


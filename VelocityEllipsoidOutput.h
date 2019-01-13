// VelocityEllipsoidOutput.h
// This class is for outputing the components of the velocity ellipsoid.  It
// is mainly used for testing of code to make sure that my results agree with
// Wisdom and Tremaine.  This version work only with the single processor model.

#ifndef VELOCITY_ELLIPSOID_OUTPUT
#define VELOCITY_ELLIPSOID_OUTPUT

#include <stdio.h>

class VelocityEllipsoidOutput {
	public:
		VelocityEllipsoidOutput() {
			fout=fopen("VelocityEllipsoid.txt","wt");
			time=0.0;
		}

		~VelocityEllipsoidOutput() {
			fclose(fout);
		}

		template<class Population>
		void output(Population &pop) {
		    int i;
			double sig1,sig2,sig3,ang;
		    double vr,vt,vz;
		    double p_rr=0.0,p_tt=0.0,p_zz=0.0,p_rt=0.0;
		    double cnt=0.0;

		    for(i=0; i<pop.getNumBodies(); i++) {
				vr=pop.gete(i)*sin(pop.getPhi(i));
				vt=(2.0-2.0*GCCoords::A0)*pop.gete(i)*cos(pop.getPhi(i));
				vz=-pop.geti(i)*sin(pop.getZeta(i));
				p_rr+=vr*vr;
				p_tt+=vt*vt;
				p_zz+=vz*vz;
				p_rt+=vr*vt;
				cnt+=1.0;
		    }
		    ang=0.5*atan2(2.0*p_rt,p_rr-p_tt);
		    sig1=(p_rr*pow(cos(ang),2.0)+p_rt*sin(2.0*ang)+p_tt*pow(sin(ang),2.0))/cnt;
		    sig2=(p_rr*pow(sin(ang),2.0)-p_rt*sin(2.0*ang)+p_tt*pow(cos(ang),2.0))/cnt;
		    sig3=p_zz/cnt;
			if((sig1<0.0) || (sig2<0.0) || (sig3<0.0)) {
				fprintf(fout,"Error calculating sigmas\n");
				fflush(fout);
				return;
			}
		    sig1=sqrt(sig1)*pop.VelocityToCMperS();
		    sig2=sqrt(sig2)*pop.VelocityToCMperS();
		    sig3=sqrt(sig3)*pop.VelocityToCMperS();
			time+=1.0;
			fprintf(fout,"%e %e %e %e %e %e\n",time*pop.getTimeStep(),sig1,sig2,sig3,sig2/sig1,sig3/sig1);
			fflush(fout);
		}
	private:
		FILE *fout;
		double time;
};

#endif
	

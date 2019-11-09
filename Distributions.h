// Distributions.h
// This file contains some classes that can produce different distributions
// of particles to be used with the populations.

#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

#include <sys/types.h>
#include <unistd.h>

//#include "GCPopulation.h"

using std::vector;

class RadiusDistrib {
	public:
		static const int MONO=0;
		static const int MODAL=1;
		static const int POWER=2;

		RadiusDistrib(double size) {
			type=MONO;
			data.resize(1);
			data[0]=size;
		}

		RadiusDistrib(double min,double max,double slope) {
			type=POWER;
			data.resize(3);
			data[2]=1.0/(1.0-slope);
			data[1]=pow(min,1-slope);
			data[0]=pow(max,1-slope)-data[1];
		}

		RadiusDistrib(int numBins) {
			type=MODAL;
			data.resize(numBins*2);
		}

		RadiusDistrib(RadiusDistrib &rd):type(rd.type),data(rd.data) {}

		void setBin(int which,double size,double fraction) {
			data[2*which]=size;
			data[2*which+1]=fraction;
		}

		double getRadius() {
			switch(type) {
				case MONO:
					return data[0];
					break;
				case MODAL: {
					double rand=drand48();
					for(unsigned int i=0; i<data.size()/2; ++i) {
						if(rand<data[2*i+1]) return data[2*i];
						rand-=data[2*i+1];
					}
					return data[data.size()-2];
				}
				case POWER:
					return pow(drand48()*data[0]+data[1],data[2]);
			}
			return 0.0;
		}
	private:
		int type;
		std::vector<double> data;
};

class CubeDistrib {
	public:
		CubeDistrib(int n,double minX,double maxX,double minY,double maxY,double minZ,double maxZ,double maxV,RadiusDistrib &rad):tot(n),minx(minX),miny(minY),maxx(maxX),maxy(maxY),minz(minZ),maxz(maxZ),maxv(maxV),rd(rad) {
			cnt=0;
			useff=false;
		}
		
		CubeDistrib(double ff,double minX,double maxX,double minY,double maxY,double minZ,double maxZ,double maxV,RadiusDistrib &rad):totff(ff),minx(minX),miny(minY),maxx(maxX),maxy(maxY),minz(minZ),maxz(maxZ),maxv(maxV),rd(rad) {
			cnt=0;
			vol=0;
			useff=true;
		}
		
		bool usesGC() {
			return false;
		}

		bool moreParticles() {
			if(useff) {
				return vol/((maxx-minx)*(maxy-miny)*(maxz-minz))<totff;
			} else {
				return cnt<tot;
			}
		}

		void setNextParticle(BasicCartCoords &cart,double &rad) {
			rad=rd.getRadius();
			cart.p[0]=minx+drand48()*(maxx-minx);
			cart.p[1]=miny+drand48()*(maxy-miny);
			cart.p[2]=minz+drand48()*(maxz-minz);
			cart.p[3]=maxv*(drand48()-0.5);
			cart.p[4]=maxv*(drand48()-0.5);
			cart.p[5]=maxv*(drand48()-0.5);
			cnt++;
			vol+=1.3333333*3.14159*rad*rad*rad;
		}

		void setNextParticle(GCCoords &gc,double &rad) {
		}

	private:
		bool useff;
		int tot,cnt;
		double totff,vol;
		double minx,miny,maxx,maxy,minz,maxz;
		double maxv;
		RadiusDistrib rd;
};

class TextFileDistrib {
	public:
		TextFileDistrib(char *fn) {
			fin=fopen(fn,"rt");
			read=fscanf(fin,"%le",&nextx);
		}
		
		bool usesGC() {
			return false;
		}

		bool moreParticles() {
			return read>0;
		}

		void setNextParticle(BasicCartCoords &cart,double &rad) {
			double y,z,vx,vy,vz;
			fscanf(fin,"%le %le %le %le %le %le",&y,&z,&vx,&vy,&vz,&rad);
			cart.p[0]=nextx;
			cart.p[1]=y;
			cart.p[2]=z;
			cart.p[3]=vx;
			cart.p[4]=vy;
			cart.p[5]=vz;
			read=fscanf(fin,"%le",&nextx);
		}

		void setNextParticle(GCCoords &gc,double &rad) {
		}
	private:
		FILE *fin;
		double nextx;
		int read;
};

template<class GCType>
class RandomSquareEIOrbits {
	public:
		RandomSquareEIOrbits(int n,double minX,double maxX,double minY,double maxY,double maxE,double maxI,RadiusDistrib &rad):tot(n),minx(minX),miny(minY),maxx(maxX),maxy(maxY),maxe(maxE),maxi(maxI),rd(rad) {
			cnt=0;
		}

		bool usesGC() {
			return true;
		}

		bool moreParticles() {
			return cnt<tot;
		}

		void setNextParticle(CartCoords &cart,double &rad) {
		}

		void setNextParticle(GCType &gc,double &rad) {
			rad=rd.getRadius();
			gc.X=minx+drand48()*(maxx-minx);
			gc.Y=miny+drand48()*(maxy-miny);
			gc.e=maxe*drand48();
			gc.i=maxi*drand48();
			gc.phi=6.28*drand48();
			gc.zeta=6.28*drand48();
			cnt++;
		}

	private:
		int tot,cnt;
		double minx,miny,maxx,maxy;
		double maxe,maxi;
		RadiusDistrib rd;
};

template <class PositionCheck,class GCType>
class RandomSquareEIOrbitsWithCheck {
	public:
		RandomSquareEIOrbitsWithCheck(int n,double minX,double maxX,double minY,double maxY,double maxE,double maxI,RadiusDistrib &rad,PositionCheck &pc):tot(n),minx(minX),miny(minY),maxx(maxX),maxy(maxY),maxe(maxE),maxi(maxI),rd(rad),posCheck(pc) {
			cnt=0;
			printf("dist=%e %e %e %e\n",minx,maxx,miny,maxy);
			printf("BC  =%e %e %e %e\n",pc.getMinX(),pc.getMaxX(),pc.getMinY(),pc.getMaxY());
		}

		bool usesGC() {
			return true;
		}

		bool moreParticles() {
			return cnt<tot;
		}

		void setNextParticle(CartCoords &cart,double &rad) {
		}

		void setNextParticle(GCType &gc,double &rad) {
			rad=rd.getRadius();
			do {
				gc.X=minx+drand48()*(maxx-minx);
				gc.Y=miny+drand48()*(maxy-miny);
			} while(!posCheck.checkIfIn(gc.X,gc.Y));
			gc.e=maxe*drand48();
			gc.i=maxi*drand48();
			gc.phi=6.28*drand48();
			gc.zeta=6.28*drand48();
			cnt++;
		}

	private:
		int tot,cnt;
		double minx,miny,maxx,maxy;
		double maxe,maxi;
		RadiusDistrib rd;
		PositionCheck &posCheck;
};

template <class PositionCheck,class GCType>
class TauRandomSquareEIOrbitsWithCheck {
	public:
		TauRandomSquareEIOrbitsWithCheck(double t,double minX,double maxX,double minY,double maxY,double maxE,double maxI,RadiusDistrib &rad,PositionCheck &pc):tau(t),minx(minX),miny(minY),maxx(maxX),maxy(maxY),maxe(maxE),maxi(maxI),rd(rad),posCheck(pc) {
			curArea=0.0;
			totArea=pc.getArea();
			printf("area=%e\n",totArea);
		}

		bool usesGC() {
			return true;
		}

		bool moreParticles() {
			return curArea/totArea<tau;
		}

		void setNextParticle(CartCoords &cart,double &rad) {
		}

		void setNextParticle(GCType &gc,double &rad) {
			rad=rd.getRadius();
			do {
				gc.X=minx+drand48()*(maxx-minx);
				gc.Y=miny+drand48()*(maxy-miny);
				curArea+=3.141592654*rad*rad;
			} while(!posCheck.checkIfIn(gc.X,gc.Y));
			gc.e=maxe*drand48();
			gc.i=maxi*drand48();
			gc.phi=6.28*drand48();
			gc.zeta=6.28*drand48();
		}

	private:
		double tau,curArea,totArea;
		double minx,miny,maxx,maxy;
		double maxe,maxi;
		RadiusDistrib rd;
		PositionCheck &posCheck;
};

template <class PositionCheck,class GCType>
class TauRandomGaussianEIOrbitsWithCheck {
	public:
		TauRandomGaussianEIOrbitsWithCheck(double t,double minX,double maxX,double minY,double maxY,double maxE,double maxI,RadiusDistrib &rad,PositionCheck &pc):tau(t),minx(minX),miny(minY),maxx(maxX),maxy(maxY),maxe(maxE),maxi(maxI),rd(rad),posCheck(pc) {
			curArea=0.0;
			totArea=(maxX-minX)*(maxY-minY);//pc.getArea();
			printf("area=%e\n",totArea);
		}

		bool usesGC() {
			return true;
		}

		bool moreParticles() {
			return curArea/totArea<tau;
		}

		void setNextParticle(CartCoords &cart,double &rad) {
		}

		void setNextParticle(GCType &gc,double &rad) {
			double midx=(maxx+minx)/2.0;
			double widthx=(maxx-minx)/6.0;
			rad=rd.getRadius();
			do {
				double r1,r2,s;
				do {
					r1=2.0*drand48()-1.0;
					r2=2.0*drand48()-1.0;
					s=r1*r1+r2*r2;
				} while(s>=1.0);
				gc.X=midx+widthx*r1*sqrt(-2.0*log(s)/s);
				gc.Y=miny+drand48()*(maxy-miny);
				curArea+=3.141592654*rad*rad;
			} while(!posCheck.checkIfIn(gc.X,gc.Y));
			gc.e=maxe*drand48();
			gc.i=maxi*drand48();
			gc.phi=6.28*drand48();
			gc.zeta=6.28*drand48();
		}

	private:
		double tau,curArea,totArea;
		double minx,miny,maxx,maxy;
		double maxe,maxi;
		RadiusDistrib rd;
		PositionCheck &posCheck;
};

/* Builds an eccentric ring.  ratio is what ffraction of the x range is the
Gaussian width.
*/
template <class GCType>
class TauRandomGaussianEccentric {
	public:
		TauRandomGaussianEccentric(double t,double minX,double maxX,double minY,double maxY,double ringE,double maxE,double ringI,double maxI,double widthRatio,RadiusDistrib &rad):tau(t),minx(minX),maxx(maxX),miny(minY),maxy(maxY),ringe(ringE),maxe(maxE),ringi(ringI),maxi(maxI),ratio(widthRatio),rd(rad) {
			curArea=0.0;
			totArea=(maxy-miny)*(maxx-minx);
			printf("area=%e\n",totArea);
		}

		bool usesGC() {
			return true;
		}

		bool moreParticles() {
			return curArea/totArea<tau;
		}

		void setNextParticle(CartCoords &cart,double &rad) {
		}

		void setNextParticle(GCType &gc,double &rad) {
			double midx=(maxx+minx)/2.0;
			double widthx=(maxx-minx)/ratio;
			rad=rd.getRadius();
			do {
				double r1,r2,s;
				do {
					r1=2.0*drand48()-1.0;
					r2=2.0*drand48()-1.0;
					s=r1*r1+r2*r2;
				} while(s>=1.0);
				gc.X=midx+widthx*r1*sqrt(-2.0*log(s)/s);
				gc.Y=miny+drand48()*(maxy-miny);
				curArea+=3.141592654*rad*rad;
			} while(gc.X<minx || gc.X>maxx);
			gc.e=ringe+maxe*drand48();
			gc.i=ringi+maxi*drand48();
			gc.phi=gc.Y+0.01*drand48();
			gc.zeta=gc.Y+0.01*drand48();
//			printf("Putting particle at %e %e %e %e %e %e\n",gc.X,gc.Y,gc.e,gc.i,gc.phi,gc.zeta);
//			printf("area %e/%e=%e\n",curArea,totArea,curArea/totArea);
//			fflush(stdout);
		}

	private:
		double tau,curArea,totArea;
		double minx,maxx;
		double miny,maxy;
		double ringe,ringi;
		double maxe,maxi;
		double ratio;
		RadiusDistrib rd;
};

class FileRecover {
	public:
		FileRecover(int step,double minX,double maxX,double minY,double maxY):minx(minX),miny(minY),maxx(maxX),maxy(maxY) {
			char buf[40];
			sprintf(buf,"CartAndRad.%d.bin",step);
			fin=fopen(buf,"rb");
			readLoc=0;
			fread(&totalParticles,sizeof(int),1,fin);
			findNextParticle();
			printf("%d - %e %e %e %e\n",totalParticles,miny,maxy,xyz.y,radius);
		}

		~FileRecover() {
			fclose(fin);
		}

		bool usesGC() {
			return false;
		}

		bool moreParticles() {
			return readLoc<totalParticles+1;
		}

#ifdef SPIN
		void setNextParticle(CartCoords &cart,double &rad,SpinVector &sp) {
#else
		void setNextParticle(CartCoords &cart,double &rad) {
#endif
			cart.x=xyz.x;
			cart.y=xyz.y;
			cart.z=xyz.z;
			cart.vx=xyz.vx;
			cart.vy=xyz.vy;
			cart.vz=xyz.vz;
			rad=radius;
			findNextParticle();
#ifdef SPIN
			sp=spin;
#endif
		}

		void setNextParticle(GCCoords &gc,double &rad) {
		}

	private:
		void findNextParticle() {
			fseek(fin,sizeof(int)+readLoc*sizeof(CartCoords),SEEK_SET);
			do {
				fread(&xyz,sizeof(CartCoords),1,fin);
				readLoc++;
//				printf("read %d (%e %e) (%e %e)\n",readLoc,xyz.x,xyz.y,miny,maxy);
			} while((xyz.y<miny || xyz.y>=maxy || xyz.x<minx*0.1) && readLoc<totalParticles+1);
			if(xyz.y>=miny && xyz.y<maxy) {
				fseek(fin,sizeof(int)+totalParticles*sizeof(CartCoords)+(readLoc-1)*sizeof(double),SEEK_SET);
				fread(&radius,sizeof(double),1,fin);
#ifdef SPIN
				fseek(fin,sizeof(int)+totalParticles*(sizeof(CartCoords)+sizeof(double))+(readLoc-1)*sizeof(SpinVector),SEEK_SET);
				fread(&spin,sizeof(SpinVector),1,fin);
#endif
			} else readLoc=totalParticles+5;
		}

		FILE *fin;
		double minx,miny,maxx,maxy;
		int totalParticles;
		int readLoc;
		CartCoords xyz;
		double radius;
#ifdef SPIN
		SpinVector spin;
#endif
};

class FileRecoverAll {
	public:
		FileRecoverAll(int step) {
			char buf[40];
			sprintf(buf,"CartAndRad.%d.bin",step);
			fin=fopen(buf,"rb");
			readLoc=0;
			fread(&totalParticles,sizeof(int),1,fin);
			findNextParticle();
#ifdef SPIN
			fseek(fin, 0, SEEK_END);
			fileLen = ftell(fin);
#endif
		}

		~FileRecoverAll() {
			fclose(fin);
		}

		bool usesGC() {
			return false;
		}

		bool moreParticles() {
			return readLoc<totalParticles+1;
		}

#ifdef SPIN
		void setNextParticle(CartCoords &cart,double &rad,SpinVector &sp) {
#else
		void setNextParticle(CartCoords &cart,double &rad) {
#endif
			cart=xyz;
			rad=radius;
#ifdef SPIN
			sp=spin;
#endif
			findNextParticle();
		}

		void setNextParticle(GCCoords &gc,double &rad) {
		}

		void setNextParticle(BasicCartCoords &cart,double &rad) {
			cart.p[0]=xyz.x;
			cart.p[1]=xyz.y;
			cart.p[2]=xyz.z;
			cart.p[3]=xyz.vx;
			cart.p[4]=xyz.vy;
			cart.p[5]=xyz.vz;
			rad=radius;
			findNextParticle();
		}

	private:
		void findNextParticle() {
			fseek(fin,sizeof(int)+readLoc*sizeof(CartCoords),SEEK_SET);
			fread(&xyz,sizeof(CartCoords),1,fin);
			readLoc++;
			fseek(fin,sizeof(int)+totalParticles*sizeof(CartCoords)+(readLoc-1)*sizeof(double),SEEK_SET);
			fread(&radius,sizeof(double),1,fin);
#ifdef SPIN
			size_t pos = sizeof(int)+totalParticles*(sizeof(CartCoords)+sizeof(double))+(readLoc-1)*sizeof(SpinVector);
			if (pos < fileLen) {
				fseek(fin, pos, SEEK_SET);
				fread(&spin,sizeof(SpinVector),1,fin);
			} else {
				spin.x = spin.y = spin.z = 0.0;
			}
#endif
		}

		FILE *fin;
		int totalParticles;
		int readLoc;
		CartCoords xyz;
		double radius;
#ifdef SPIN
		long fileLen;
		SpinVector spin;
#endif
};

class FileRecoverBasic {
	public:
		FileRecoverBasic(int step) {
			char buf[40];
			sprintf(buf,"CartAndRad.%d.bin",step);
			fin=fopen(buf,"rb");
			readLoc=0;
			fread(&totalParticles,sizeof(int),1,fin);
			findNextParticle();
		}

		~FileRecoverBasic() {
			fclose(fin);
		}

		bool usesGC() {
			return false;
		}

		bool moreParticles() {
			return readLoc<totalParticles+1;
		}

		void setNextParticle(BasicCartCoords &cart,double &rad) {
			cart.p[0]=xyz.x;
			cart.p[1]=xyz.y;
			cart.p[2]=xyz.z;
			cart.p[3]=xyz.vx;
			cart.p[4]=xyz.vy;
			cart.p[5]=xyz.vz;
			rad=radius;
			findNextParticle();
		}

		void setNextParticle(GCCoords &gc,double &rad) {
		}

	private:
		void findNextParticle() {
			fseek(fin,sizeof(int)+readLoc*sizeof(CartCoords),SEEK_SET);
			fread(&xyz,sizeof(CartCoords),1,fin);
			readLoc++;
			fseek(fin,sizeof(int)+totalParticles*sizeof(CartCoords)+(readLoc-1)*sizeof(double),SEEK_SET);
			fread(&radius,sizeof(double),1,fin);
		}

		FILE *fin;
		int totalParticles;
		int readLoc;
		CartCoords xyz;
		double radius;
};

// An inherently unbalanced distribution, used to test load balancing
// Added by mmaly
template <class PositionCheck>
class UnbalancedTauRandomSquareEIOrbitsWithCheck {
        public:
                UnbalancedTauRandomSquareEIOrbitsWithCheck(double t,double minX,double maxX,double minY,double maxY,double maxE,double maxI,RadiusDistrib &rad,PositionCheck &pc):tau(t),minx(minX),miny(minY),maxx(maxX),maxy(maxY),maxe(maxE),maxi(maxI),rd(rad),posCheck(pc) {
                        curArea=0.0;
                        totArea=(maxX-minX)*(maxY-minY);
                        printf("area=%e\n",totArea);
                }

                bool usesGC() {
                        return true;
                }

                bool moreParticles() {
                        return curArea/totArea<tau;
                }

                void setNextParticle(CartCoords &cart,double &rad) {
                }

                void setNextParticle(GCCoords &gc,double &rad) {
                        rad=rd.getRadius();
                        do {
                                gc.X=minx+drand48()*(maxx-minx);
				double input = 2*(drand48() - 0.5);
				//gc.Y=miny+(4*input*input*input + 0.5)*(maxy-miny);
				gc.Y = (input*input*input + input*0.1/1.1)*(maxy-miny);
                                curArea+=3.141592654*rad*rad;
                        } while(!posCheck.checkIfIn(gc.X,gc.Y));
                        
			gc.e=maxe*drand48();
                        gc.i=maxi*drand48();
                        gc.phi=6.28*drand48();
                        gc.zeta=6.28*drand48();
                }

        private:
                double tau,curArea,totArea;
                double minx,miny,maxx,maxy;
                double maxe,maxi;
                RadiusDistrib rd;
                PositionCheck &posCheck;
};


#ifdef PARALLEL

/**
 * This class was created because using the normal file loading with a large
 * parallel simulation seems to swamp the system because we have many processes
 * reading a large file all at once.  This version has only the master read the
 * file and send out the information.  Because of the way distributions normally
 * work, all of the reading and spreading around of that information must be
 * done during the construction of these objects.  The handing out of the
 * particles is done after each processor knows what it has.
 *
 * When calling this, make sure to pass in the total bounds, not just the
 * local bounds.  The procReadAhead argument is used to compensate for changing
 * the number of processors.  It can be 1 if the number of processors stays the
 * same (technically you can try 0 in that case to get smaller buffers).  If
 * the number of processors is changed it must be at least 1 and should be
 * higher if the number of processors is increased significantly.
 */
template<class GCType>
class ParallelFileRecover {
	public:
		ParallelFileRecover(int step,int procReadAhead,double minX,double maxX,double minY,double maxY,ProcessorCommunication &pc):minx(minX),miny(minY),maxx(maxX),maxy(maxY),cur(0) {
			if(pc.getProcessNum()==0) {
				int writeProc=1;
				int totalParticles;
				char buf[40];
				sprintf(buf,"CartAndRad.%d.bin",step);
				FILE *fin=fopen(buf,"rb");
				fread(&totalParticles,sizeof(int),1,fin);
				printf("Reading %d particles.\n",totalParticles);
				CartCoords coord;
				vector<CartCoords> buffer(2*totalParticles/pc.getNumProcesses());
				vector<double> cartMessage(6*totalParticles/pc.getNumProcesses());
				vector<double> radiusMessage(totalParticles/pc.getNumProcesses());
				vector<int> partNum(totalParticles/pc.getNumProcesses());
				vector<int> proc(totalParticles/pc.getNumProcesses());
				buffer.resize(0);
				partNum.resize(0);
				proc.resize(0);
				int firstFree=-1;
				for(int i=0; i<totalParticles; ++i) {
					fread(&coord,sizeof(CartCoords),1,fin);
					int curProc=(int)(pc.getNumProcesses()*(coord.y-miny)/(maxy-miny));
					if(curProc<0) curProc=0;
					else if(curProc>=pc.getNumProcesses()) curProc=pc.getNumProcesses()-1;
					if(curProc==0) {
						xyz.push_back(coord);
						fseek(fin,sizeof(int)+(unsigned)totalParticles*sizeof(CartCoords)+(unsigned)i*sizeof(double),SEEK_SET);
						double tmp;
						fread(&tmp,sizeof(double),1,fin);
						radius.push_back(tmp);
						fseek(fin,sizeof(int)+(unsigned)(i+1)*sizeof(CartCoords),SEEK_SET);
					} else {
						if(firstFree<0) {
							proc.push_back(curProc);
							partNum.push_back(i);
							buffer.push_back(coord);
//							printf("pushed to buffer %d for %d at %d = %e %e\n",i,curProc,buffer.size(),buffer[buffer.size()-1].x,buffer[buffer.size()-1].y);
						} else {
							int j=firstFree;
							firstFree=-proc[j];
							proc[j]=curProc;
							partNum[j]=i;
							buffer[j]=coord;
//							printf("stored at %d %d for %d = %e %e\n",j,i,curProc,buffer[buffer.size()-1].x,buffer[buffer.size()-1].y);
						}
					}
					if(curProc>writeProc+procReadAhead) {
						firstFree=sendFromBuffer(writeProc,partNum,proc,buffer,cartMessage,radiusMessage,firstFree,fin,totalParticles,pc);
						fseek(fin,sizeof(int)+(unsigned)(i+1)*sizeof(CartCoords),SEEK_SET);
						writeProc++;
					}
				}
				for(; writeProc<pc.getNumProcesses(); ++writeProc) {
					firstFree=sendFromBuffer(writeProc,partNum,proc,buffer,cartMessage,radiusMessage,firstFree,fin,totalParticles,pc);
				}
				fclose(fin);
			} else {
				int numRead;
				vector<double> size;
				vector<double> buffer;
				size.resize(1);
				printf("Receiving in %d\n",pc.getProcessNum());
				fflush(stdout);
				pc.readFrom(0,size,numRead);
				printf("size is %d\n",(int)size[0]);
				buffer.resize((int)(6*size[0]));
				printf("Buffer is %d.\n",buffer.size());
				pc.readFrom(0,buffer,numRead);
				xyz.resize((int)(size[0]));
				for(unsigned int i=0; i<xyz.size(); ++i) {
					xyz[i].x=buffer[i*6];
					xyz[i].y=buffer[i*6+1];
					xyz[i].z=buffer[i*6+2];
					xyz[i].vx=buffer[i*6+3];
					xyz[i].vy=buffer[i*6+4];
					xyz[i].vz=buffer[i*6+5];
				}
				radius.resize((int)(size[0]));
				printf("Radius is %d.\n",radius.size());
				pc.readFrom(0,radius,numRead);
				printf("Finished reading.\n");
			}
		}

		bool usesGC() {
			return false;
		}

		bool moreParticles() {
			return cur<xyz.size();
		}

		void setNextParticle(CartCoords &cart,double &rad) {
			cart.x=xyz[cur].x;
			cart.y=xyz[cur].y;
			cart.z=xyz[cur].z;
			cart.vx=xyz[cur].vx;
			cart.vy=xyz[cur].vy;
			cart.vz=xyz[cur].vz;
			rad=radius[cur];
			cur++;
		}

		void setNextParticle(GCType &gc,double &rad) {
		}

	private:
		int sendFromBuffer(int writeProc,vector<int> &partNum,vector<int> &proc,vector<CartCoords> &buffer,vector<double> &cartMessage,vector<double> &radiusMessage,int firstFree,FILE *fin,int totalParticles,ProcessorCommunication &pc) {
			printf("Sending buffer to %d, %d, %d.\n",writeProc,proc.size(),totalParticles);
			cartMessage.resize(0);
			radiusMessage.resize(0);
			for(unsigned int i=0; i<proc.size(); ++i) {
				if(proc[i]==writeProc) {
					cartMessage.push_back(buffer[i].x);
					cartMessage.push_back(buffer[i].y);
					cartMessage.push_back(buffer[i].z);
					cartMessage.push_back(buffer[i].vx);
					cartMessage.push_back(buffer[i].vy);
					cartMessage.push_back(buffer[i].vz);
//					printf("%e==%e\n",buffer[i].x,cartMessage[cartMessage.size()-6]);
					fseek(fin,sizeof(int)+(unsigned)totalParticles*sizeof(CartCoords)+(unsigned)partNum[i]*sizeof(double),SEEK_SET);
					double tmp;
					fread(&tmp,sizeof(double),1,fin);
//					printf("%d radius[%d]=%e at %ud %d\n",i,partNum[i],tmp,sizeof(int)+(unsigned)totalParticles*sizeof(CartCoords)+(unsigned)partNum[i]*sizeof(double),sizeof(unsigned));
					radiusMessage.push_back(tmp);
					proc[i]=-firstFree;
					firstFree=i;
//					printf("Message x=%e y=%e r=%e\n",cartMessage[cartMessage.size()-6],cartMessage[cartMessage.size()-5],radiusMessage[radiusMessage.size()-1]);
				}
			}
			vector<double> size;
			size.push_back(radiusMessage.size());
			printf("Sending messages %d, %d, %d, %d\n",(int)size[0],size.size(),cartMessage.size(),radiusMessage.size());
			pc.sendTo(writeProc,size);
			pc.sendTo(writeProc,cartMessage);
			pc.sendTo(writeProc,radiusMessage);
			printf("Sent %d, %d, %d, %d\n",(int)size[0],size.size(),cartMessage.size(),radiusMessage.size());
			return firstFree;
		}

		double minx,miny,maxx,maxy;
		int cur;
		vector<CartCoords> xyz;
		vector<double> radius;
};

/**
 * This class was created because using the normal file loading with a large
 * parallel simulation seems to swamp the system because we have many processes
 * reading a large file all at once.  This version has only the master read the
 * file and send out the information.  Because of the way distributions normally
 * work, all of the reading and spreading around of that information must be
 * done during the construction of these objects.  The handing out of the
 * particles is done after each processor knows what it has.  This version can
 * also handle big files on a 32-bit system because it doesn't use fseek.
 *
 * When calling this, make sure to pass in the total bounds, not just the
 * local bounds.  The procReadAhead argument is used to compensate for changing
 * the number of processors.  It can be 1 if the number of processors stays the
 * same (technically you can try 0 in that case to get smaller buffers).  If
 * the number of processors is changed it must be at least 1 and should be
 * higher if the number of processors is increased significantly.
 */
class ParallelFileRecover32 {
	public:
		ParallelFileRecover32(int step,int procReadAhead,double minX,double maxX,double minY,double maxY,ProcessorCommunication &pc):minx(minX),miny(minY),maxx(maxX),maxy(maxY),cur(0) {
			if(pc.getProcessNum()==0) {
				int writeProc=1;
				int totalParticles;
				char buf[40];
				sprintf(buf,"CartAndRad.%d.bin",step);
				FILE *fin=fopen(buf,"rb");
				fread(&totalParticles,sizeof(int),1,fin);
				printf("Reading %d particles.\n",totalParticles);
				CartCoords coord;
				vector<int> totProc(totalParticles);
				{
					vector<CartCoords> buffer(2*totalParticles/pc.getNumProcesses());
					vector<double> cartMessage(6*totalParticles/pc.getNumProcesses());
					vector<int> proc(totalParticles/pc.getNumProcesses());
					buffer.resize(0);
					proc.resize(0);
					int firstFree=-1;
					for(int i=0; i<totalParticles; ++i) {
						fread(&coord,sizeof(CartCoords),1,fin);
						int curProc=(int)(pc.getNumProcesses()*(coord.y-miny)/(maxy-miny));
						if(curProc<0) curProc=0;
						else if(curProc>=pc.getNumProcesses()) curProc=pc.getNumProcesses()-1;
						totProc[i]=curProc;
						if(curProc==0) {
							xyz.push_back(coord);
						} else {
							if(firstFree<0) {
								proc.push_back(curProc);
								buffer.push_back(coord);
							} else {
								int j=firstFree;
								firstFree=-proc[j];
								proc[j]=curProc;
								buffer[j]=coord;
							}
						}
						if(curProc>writeProc+procReadAhead) {
							firstFree=sendFromBuffer(writeProc,proc,buffer,cartMessage,firstFree,fin,totalParticles,pc);
							writeProc++;
						}
					}
					for(; writeProc<pc.getNumProcesses(); ++writeProc) {
						firstFree=sendFromBuffer(writeProc,proc,buffer,cartMessage,firstFree,fin,totalParticles,pc);
					}
				}
				sendRadii(totProc,fin,totalParticles,pc);
				fclose(fin);
			} else {
				int numRead;
				vector<double> size;
				vector<double> buffer;
				size.resize(1);
				printf("Receiving in %d\n",pc.getProcessNum());
				fflush(stdout);
				pc.readFrom(0,size,numRead);
				printf("%d size is %d\n",pc.getProcessNum(),(int)size[0]);
				buffer.resize((int)(6*size[0]));
				pc.readFrom(0,buffer,numRead);
				xyz.resize((int)(size[0]));
				for(unsigned int i=0; i<xyz.size(); ++i) {
					xyz[i].x=buffer[i*6];
					xyz[i].y=buffer[i*6+1];
					xyz[i].z=buffer[i*6+2];
					xyz[i].vx=buffer[i*6+3];
					xyz[i].vy=buffer[i*6+4];
					xyz[i].vz=buffer[i*6+5];
				}
				printf("%d x=%e y=%e\n",pc.getProcessNum(),xyz[0].x,xyz[0].y);
				radius.resize((int)(size[0]));
				pc.readFrom(0,radius,numRead);
				printf("%d r=%e\n",pc.getProcessNum(),radius[0]);
				printf("%d Finished reading.\n",pc.getProcessNum());
			}
		}

		bool usesGC() {
			return false;
		}

		bool moreParticles() {
			return cur<xyz.size();
		}

		void setNextParticle(CartCoords &cart,double &rad) {
			cart.x=xyz[cur].x;
			cart.y=xyz[cur].y;
			cart.z=xyz[cur].z;
			cart.vx=xyz[cur].vx;
			cart.vy=xyz[cur].vy;
			cart.vz=xyz[cur].vz;
			rad=radius[cur];
			cur++;
		}

		void setNextParticle(GCCoords &gc,double &rad) {
		}

	private:
		int sendFromBuffer(int writeProc,vector<int> &proc,vector<CartCoords> &buffer,vector<double> &cartMessage,int firstFree,FILE *fin,int totalParticles,ProcessorCommunication &pc) {
			cartMessage.resize(0);
			int numSent=0;
			for(unsigned int i=0; i<proc.size(); ++i) {
				if(proc[i]==writeProc) {
					cartMessage.push_back(buffer[i].x);
					cartMessage.push_back(buffer[i].y);
					cartMessage.push_back(buffer[i].z);
					cartMessage.push_back(buffer[i].vx);
					cartMessage.push_back(buffer[i].vy);
					cartMessage.push_back(buffer[i].vz);
					proc[i]=-firstFree;
					firstFree=i;
					numSent++;
				}
			}
			vector<double> size;
			size.push_back(numSent);
			pc.sendTo(writeProc,size);
			pc.sendTo(writeProc,cartMessage);
			return firstFree;
		}

		void sendRadii(vector<int> &proc,FILE *fin,int totalParticles,ProcessorCommunication &pc) {
			vector<double> radii(totalParticles);
			vector<double> message;
			fread(&radii[0],sizeof(double),totalParticles,fin);
			for(unsigned int j=0; j<totalParticles; ++j) {
				if(proc[j]==0) radius.push_back(radii[j]);
			}
			for(int i=1; i<pc.getNumProcesses(); ++i) {
				message.resize(0);
				double min=1e100,max=-1e100;
				for(unsigned int j=0; j<totalParticles; ++j) {
					if(proc[j]==i) {
						message.push_back(radii[j]);
						if(radii[j]<min) min=radii[j];
						if(radii[j]>max) max=radii[j];
					}
				}
				printf("Send radii to %d, min=%e, max=%e, size=%d\n",i,min,max,message.size());
				pc.sendTo(i,message);
			}
		}

		double minx,miny,maxx,maxy;
		int cur;
		vector<CartCoords> xyz;
		vector<double> radius;
};

#endif

#endif

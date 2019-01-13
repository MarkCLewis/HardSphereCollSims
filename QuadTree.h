#ifndef QUADTREE
#define QUADTREE

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

class Node {
	
	public:	
		Node()
		{
			//////////cout << "Creating a node" << endl;
			m_iPopIndex = -1;
			m_iParentIndex = -1;
			m_vChildren.resize(4);
			m_vChildren[0] = -1;
			m_vChildren[1] = -1;
			m_vChildren[2] = -1;
			m_vChildren[3] = -1;
			divX = 0;
			divY = 0;
		}

		
		int getQuadrant(vector<double> &point) {
			int ret=0;
			if(point[0]>divX) { ret|=1; }
			if(point[1]>divY) { ret|=2; }
			return ret;
		}

		int getPopIndex() const { return m_iPopIndex; }
		int getParentIndex() const { return m_iParentIndex; } 
		double getDivisionX() { return divX; }
		double getDivisionY() { return divY; }

		int getChildValue(int index) { return m_vChildren[index]; }

		bool isEmpty() const 
		{ 
			if(m_iPopIndex == -1)
				return true;
			else
				return false;
		}

		template <class Population>
		void setPopIndex(Population &pop,int i) { m_iPopIndex = i; divX=pop.getx(i); divY=pop.gety(i); }
		void setParentIndex(int i) { m_iParentIndex = i; }

		void setChildValue(int index,int val) {
			m_vChildren[index]=val;
		}

		bool isLeaf()
		{
			return (m_vChildren[0] == -1) && (m_vChildren[1] == -1) && (m_vChildren[2] == -1) && (m_vChildren[3] == -1); 
		}

	private:

		int m_iPopIndex;
		double divX;
		double divY;
		int m_iParentIndex;
		bool m_bEmpty;
		vector<int> m_vChildren;
};
	
class quadtree {

	public:
		quadtree()
		{
			builds = 0;
			radius = 0;
			stop = 0;
		}
		

	// Build a quadtree representation of the population
	// This will call virtually all the functions in this
	// class.
	template<class Population>
	void build(Population &SamplePop)
	{
	        int num = SamplePop.getNumBodies();
		qtree.clear();
	        qtree.resize( (num * 2) );
		m_iRoot = 0;
		FreeNode = 0;
		// radius = 2 * particleRadius + 3*RMS*dt
		if( (builds%100) == 0 )
			radius = (2 * SamplePop.getMaxParticleRadius() ) + (3 * getRMSVelocity(SamplePop) * SamplePop.getTimeStep() ) ;
		builds = builds + 1;
		vector<double> coord(3);
	        for(int i = 0; i < num; i++)
	        {
			coord[0] = SamplePop.getx(i);
			coord[1] = SamplePop.gety(i);
			coord[2] = SamplePop.getz(i);
			Insert(SamplePop, i, coord, m_iRoot);
	        }
		getCoordFile(builds, num);
	}

	template<class TreeForcing,class Population>
	void findAllCollisions(Population &pop, TreeForcing &tf)
	{
		int index;
		int ignore = -1;
		for(int i = 0; i<pop.getNumBodies();i++)
		{
			index = findNode(pop, i, m_iRoot);
			checkSearchRadius(pop,tf, m_iRoot, index, radius, ignore);
		//	cout << "i= " << i << endl;
		}
		//cout << "find all collisions method is over" << endl;
	}
	
		
	template<class TreeForcing,class Population>
	void findCollisions(Population &pop,TreeForcing &tf,int particle, int ignore)
	{
		//cout << "Find collisions called" << endl;
		int index = findNode(pop, particle, m_iRoot);
		checkSearchRadius(pop,tf,m_iRoot, index, radius, ignore);
		//cout << "done with findCollisions, particle= "<< particle << "ignore = " << ignore << endl;
	}

	private:

	// radius = 2 * particleRadius + 3*RMS*dt
	// Done Coding.. need testing
	template<class Population>
	double getRMSVelocity(Population &pop) 
	{
		return 0.3;

/*		double vms, temp;
		int index, parent;
		vms = 0;
		temp = 0;
		double sum = 0;
		int n = 100; 

		if( pop.getNumBodies() < 100) n = pop.getNumBodies();
	
		srand((unsigned)time(NULL));

		for(int i = 0; i < n; i++)
		{
			// Make sure you pick a valid index
			// if a bad index is picked.. keep 
			// picking random numbers till you find
			// a good one.
			do
				{
					index = int(n * (rand()/(RAND_MAX+1.0)));	
				}
			while( qtree[index].getPopIndex() != -1 );

			if(index == 0) index = 1;
			parent = qtree[index].getParentIndex();
			temp = pow((pop.getvx(parent) - pop.getvx(index)),2);
			temp += pow((pop.getvy(parent) - pop.getvy(index)),2);
			temp += pow((pop.getvz(parent) - pop.getvz(index)),2);
			vms = 0.5 * sqrt(temp);
			sum += vms;
			
		}
	
		double average = sum/n;
	
		
		return average; 
*/	}

		// Recursive Insertation of data into the quadtree
		// I might later replace this with an iterative insertation
		// method to get a slight speed increase, but we'll see.
		// For testing purposes this will do.
		template<class Population>
		bool Insert(Population &pop, int &popIndex, vector<double> &point, int parent)
		{
			int Q = qtree[parent].getQuadrant(point);

			if(qtree[parent].getChildValue(Q) != -1)
			{
				Insert(pop, popIndex, point, qtree[parent].getChildValue(Q) );
			}
			else
			{
				int index = getFreeCell();
				if(parent != index)
				{
					qtree[parent].setChildValue(Q,index);
					qtree[index].setPopIndex(pop,popIndex);
					qtree[index].setParentIndex(parent);
					return true;
				}
				// Special case for the root
				if(index == m_iRoot)
				{
					qtree[index].setPopIndex(pop,popIndex);
					qtree[index].setParentIndex(-1);
					return true;
				}
				
				else return false;
			}
			return false;
		}

		// UNDER CONSTRUCTION
		template<class TreeForcing,class Population>
		void checkSearchRadius(Population &pop, TreeForcing &cm, int current, int center, double radius, int ignore)
		{
			//cout << "in the overlap function" << endl;
			if(current >= 0 && (qtree[current].isLeaf() == false) && (current != center) )
			{
				//cout << "in first block " << endl;
				double diffSquared, rSquared;
				int quad = 0;
				rSquared = radius * radius;

				diffSquared = qtree[current].getDivisionX() - qtree[center].getDivisionX();
				diffSquared = diffSquared * diffSquared;

				if(diffSquared <= rSquared) // 0011
					quad = quad ^ 3;
				else // 0001
					quad = quad ^ 1;
				
				diffSquared = qtree[current].getDivisionY() - qtree[center].getDivisionY();
				diffSquared = diffSquared * diffSquared;
				
				if(diffSquared <= rSquared) // 1100
					quad = quad ^ 12;
				else // 0100
					quad = quad ^ 4;
				
				//cout << "quad = " << quad << endl;
				int kid = 0;
				if(quad == 15)
				{
					// recurse down all of children
					kid = qtree[current].getChildValue(0);
					if(kid != -1)
						checkSearchRadius(pop,cm,kid, center, radius, ignore);

					kid = qtree[current].getChildValue(1);
					if(kid != -1)
						checkSearchRadius(pop,cm,kid, center, radius, ignore);

					kid = qtree[current].getChildValue(2);
					if(kid != -1)
						checkSearchRadius(pop,cm,kid, center, radius, ignore);

					kid = qtree[current].getChildValue(3);
					if(kid != -1)
						checkSearchRadius(pop,cm,kid, center, radius, ignore);
				}

				if(quad == 7)
				{
					vector<double> a(2);
					a[0] = qtree[center].getDivisionX();
					a[1] = qtree[center].getDivisionY();
					int Q = qtree[current].getQuadrant(a);
					// Recurse along XN,YN, XO
					kid = qtree[current].getChildValue(Q);
					if(kid != -1)
						checkSearchRadius(pop, cm,kid, center, radius, ignore);

					Q = Q ^ 1;

					kid = qtree[current].getChildValue(Q);
					if(kid != -1)
						checkSearchRadius(pop,cm, kid, center, radius, ignore);
				}
				
				if(quad == 13)
				{
					vector<double> a(2);
					a[0] = qtree[center].getDivisionX();
					a[1] = qtree[center].getDivisionY();
					int Q = qtree[current].getQuadrant(a);
					//Recurse along XN, YN, YO
					kid = qtree[current].getChildValue(Q);
					if(kid != -1)
						checkSearchRadius(pop,cm,kid, center, radius, ignore);
					Q = Q ^ 2;
					kid = qtree[current].getChildValue(Q);
					if(kid != -1)
						checkSearchRadius(pop,cm,kid, center, radius, ignore);
				}

				if(quad == 5)
				{
					vector<double> a(2);
					a[0] = qtree[center].getDivisionX();
					a[1] = qtree[center].getDivisionY();
					////cout << "a = " << a[0] << " , " << a[1] << endl;
					int Q = qtree[current].getQuadrant(a);
					// Recurse down XN, YN
					kid = qtree[current].getChildValue(Q);
					////cout << "Q=" << Q << endl;
					////cout << "kid=" << kid << endl;
					if(kid != -1)
					{
						////cout << "Args: \n" << "kid = " << kid << "\ncenter= " << center;
						////cout << "\nradius = " << radius << "\nignore = "<< ignore << endl;
						checkSearchRadius(pop,cm,kid, center, radius, ignore);
					}
				}

				if(qtree[current].getPopIndex() != ignore)
				{
				//	cout << "first half | ";
				//	cout <<  "ignore: " << ignore << " | current: " << qtree[current].getPopIndex();
				//	cout << " | center: " << qtree[center].getPopIndex() << endl;
                                	double t=pop.collisionTime(qtree[current].getPopIndex(),qtree[center].getPopIndex());
                                	if( (t>=0.0) && (t<pop.getTimeStep()) ) 
					{
				        	if(ignore == qtree[current].getPopIndex()) cout << "Ignore and current match!" << endl;
                                	        cm.addPotentialCollision(qtree[current].getPopIndex(),qtree[center].getPopIndex(),t);
                                	}
				}
					
			}
			//else if it is a leaf
			else if( (current >= 0) && (center != current))
			{
				if(qtree[current].getPopIndex() !=  ignore)
				{
				//	if(ignore == qtree[current].getPopIndex()) cout << "Ignore and current match!" << endl;
                                	double t = pop.collisionTime(qtree[current].getPopIndex(),qtree[center].getPopIndex());
                                	if( (t>=0.0) && (t<pop.getTimeStep()) ) 
                                	        cm.addPotentialCollision(qtree[current].getPopIndex(),qtree[center].getPopIndex(),t);
				}
			}
		}

		//  Done
		template<class Population>
		int findNode(Population &pop, int popIndex, int current)
		{
			return popIndex;
				if(qtree[current].getPopIndex() == popIndex)
				{
					return current;
				}
				
				vector<double> p;
				p.resize(2);
	
				p[0] = pop.getx(popIndex);
					p[1] = pop.gety(popIndex);
				
				int Q = qtree[current].getQuadrant(p);
			

				return findNode(pop, popIndex, qtree[current].getChildValue(Q));
		}
		
		int getFreeCell(void)
		{
			int tmp = FreeNode;	
			FreeNode = FreeNode + 1;
			if(FreeNode >= qtree.size()) 
			{
				FreeNode = m_iRoot;
			}
			return tmp;
		};

		void getCoordFile(int buildNo, int num)
		{
	
			string name = "Coord";
			//itoa(buildNo);
			ostringstream oss;
			oss << buildNo;
			string no(oss.str());
			name = name+no;
			ofstream coordfp(name.c_str());
			if(coordfp.is_open())
			{
				coordfp << "#X\t Y" << endl;
				for(int i = 0; i < num ; i++)
					coordfp << qtree[i].getDivisionX() << "\t " << qtree[i].getDivisionY() << endl;
				coordfp.close();
			}
	
		}
		int builds;
		double radius;
		int FreeNode;
		vector<Node> qtree;
		int m_iRoot;

		// debugging variable
		int stop;
};

#endif

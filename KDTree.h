#ifndef KDTREE
#define KDTREE

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;

class Node {
	
	public:	
		Node()
		{
			m_iPopIndex = -1;
			m_iParentIndex = -1;
			m_vChildren.resize(2);
			m_vChildren[0] = -1;
			m_vChildren[1] = -1;
			m_iDivX = 0;
			m_iDivY = 0;
			m_iSplitDim = 0;
		}

		
		int getQuadrant(vector<double> &point) {
			int ret=0;
			if(m_iSplitDim == 1) // Compare X
			{
			 	if(point[0]>m_iDivX) 
					ret = 1;
				else  ret = 0; 
			}
			else if(m_iSplitDim == 0) // Compare Y
			{
				if(point[1]>m_iDivY)
					ret = 1;
				else ret = 0;
			}
			
			return ret;
		}

		int getPopIndex() const { return m_iPopIndex; }
		int getParentIndex() const { return m_iParentIndex; } 
		double getDivisionX() { return m_iDivX; }
		double getDivisionY() { return m_iDivY; }
		double getSplitDimension() { return m_iSplitDim; }

		int getChildValue(int index) { return m_vChildren[index]; }

		bool isEmpty() const 
		{ 
			if(m_iPopIndex == -1)
				return true;
			else
				return false;
		}

		template <class Population>
		void setPopIndex(Population &pop,int i) { m_iPopIndex = i; m_iDivX=pop.getx(i); m_iDivY=pop.gety(i); }
		void setParentIndex(int i) { m_iParentIndex = i; }

		void setChildValue(int index,int val) {
			m_vChildren[index]=val;
		}

		bool isLeaf()
		{
			return (m_vChildren[0] == -1) && (m_vChildren[1] == -1); 
		}

	private:

		int m_iPopIndex;
		double m_iDivX;
		double m_iDivY;
		int m_iSplitDim;
		int m_iParentIndex;
		bool m_bEmpty;
		vector<int> m_vChildren;
};
	
class KDTree {

	public:
		KDTree()
		{
			builds = 0;
			radius = 0;
			stop = 0;
		}
		

	// Build a KDTree representation of the population
	// This will call virtually all the functions in this
	// class.
	template<class Population>
	void build(Population &SamplePop)
	{
	        int num = SamplePop.getNumBodies();
		kdtree.clear();
	        kdtree.resize( (num * 2) );
		m_iRoot = 0;
		FreeNode = 0;

		if( (builds%100) == 0 )
			radius = (2 * SamplePop.getMaxParticleRadius() ) + (3 * getRMSVelocity(SamplePop) * SamplePop.getTimeStep() ) ;
		builds = builds + 1;

		vector<double> coord(3);
	        for(int i = 0; i < num; i++)
	        {
			coord[0] = SamplePop.getx(i);
			coord[1] = SamplePop.gety(i);
			coord[2] = SamplePop.getz(i);
			insert(SamplePop, coord, i);
	        }
	};

	template<class Population>
	double getRMSVelocity(Population &pop) 
	{
		double vms, temp;
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
			while( kdtree[index].getPopIndex() != -1 );

			if(index == 0) index = 1;
			parent = kdtree[index].getParentIndex();
			temp = pow((pop.getvx(parent) - pop.getvx(index)),2);
			temp += pow((pop.getvy(parent) - pop.getvy(index)),2);
			temp += pow((pop.getvz(parent) - pop.getvz(index)),2);
			vms = 0.5 * sqrt(temp);
			sum += vms;
			
		}
	
		double average = sum/n;
	
		
		return average; 
	}

	template<class TreeForcing,class Population>
	void findAllCollisions(Population &pop, TreeForcing &tf)
	{
		findInitialPotentialCollisions(pop,tf);
	}
	
		
	template<class TreeForcing,class Population>
	void findCollisions(Population &pop,TreeForcing &tf,int particle)
	{
		findPotentialCollisionsForParticle(pop,tf, particle);
	}

	template<class TreeForcing,class Population>
	void findInitialPotentialCollisions(Population &pop,TreeForcing &tf)
	{
		for(int i = 0; i<pop.getNumBodies();i++)
		{
			findPotentialCollisionsForParticle(pop,tf,i);
		}
	 }
	
	template<class TreeForcing,class Population>
	void findPotentialCollisionsForParticle(Population &pop,TreeForcing &tf,int particle)
	{ 
		// radius = 2 * particleRadius + 3*RMS*dt
		vector<int> list =  getCollisions(pop, tf, particle, radius);
		if(list.size() > 0) stop = stop + 1;
		////cout << "******************************** COLLISIONS RESULTS **********************" << endl;
		////cout << "stop= " << stop << endl;
		//////cout << " Possible collisions: " << list.size() << endl;
		////cout << "**************************************************************************" << endl;
	}
	
		

	private:

		// Simple wrapper function for an iterative or 
		// recursive insertation method.
		template<class Population>
		void insert(Population &pop, vector<double> &point, int &popIndex)
		{
		//	////////cout << "Entering special insert function" << endl;
		//	////////cout << "Qtree.size = " << kdtree.size() << endl;
			recursiveInsert(pop, popIndex, point, m_iRoot);
		};

		// Recursive Insertation of data into the KDTree
		// I might later replace this with an iterative insertation
		// method to get a slight speed increase, but we'll see.
		// For testing purposes this will do.
		template<class Population>
		bool recursiveInsert(Population &pop, int &popIndex, vector<double> &point, int parent)
		{
		//	////////cout << "root=" << m_iRoot << endl;
			int Q = kdtree[parent].getQuadrant(point);
		//	////////cout << "Quadrant is " << Q << endl;

			if(kdtree[parent].getChildValue(Q) != -1)
			{
				recursiveInsert(pop, popIndex, point, kdtree[parent].getChildValue(Q) );
			}
			else
			{
//				////////cout << "Adding element " << endl;
				int index = getFreeCell();
				if(parent != index)
				{
					kdtree[parent].setChildValue(Q,index);
					kdtree[index].setPopIndex(pop,popIndex);
					kdtree[index].setParentIndex(parent);
					return true;
				}
				// Special case for the root
				if(index == m_iRoot)
				{
			//		kdtree[parent].setChildValue(Q,index);
					kdtree[index].setPopIndex(pop,popIndex);
					kdtree[index].setParentIndex(-1);
					return true;
				}
				
				else return false;
			}
		//	////////cout << "Warning!!! Somehow the KDTree vector is full!!!" << endl;
			return false;
		};


		template<class TreeForcing, class Population>
		vector<int> getCollisions(Population &pop, TreeForcing &cm, int popIndex, double radius)
		{
			//Find the Node
			int index = findNode(pop, popIndex, m_iRoot);
			vector<int> a;
			findOverlaps(pop,cm,m_iRoot, index, radius, a);
			return a;
		}
	
		// UNDER CONSTRUCTION
		template<class TreeForcing,class Population>
		void findOverlaps(Population &pop, TreeForcing &cm, int current, int center, double radius, vector<int> & list)
		{
			if(current >= 0 && (kdtree[current].isLeaf() == false) )
			{
				double diffSquared, rSquared;
				int quad = 0;
				rSquared = radius * radius;

				diffSquared = kdtree[current].getDivisionX() - kdtree[center].getDivisionX();
				diffSquared = diffSquared * diffSquared;

				if(diffSquared <= rSquared) // 0011
					quad = quad ^ 3;
				else // 0001
					quad = quad ^ 1;
				
				diffSquared = kdtree[current].getDivisionY() - kdtree[center].getDivisionY();
				diffSquared = diffSquared * diffSquared;
				
				if(diffSquared <= rSquared) // 1100
					quad = quad ^ 12;
				else // 0100
					quad = quad ^ 4;
				
				int kid;
				if(quad == 15)
				{
					// recurse down all of children
					kid = kdtree[current].getChildValue(0);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);

					kid = kdtree[current].getChildValue(1);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);

					kid = kdtree[current].getChildValue(2);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);

					kid = kdtree[current].getChildValue(3);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);
				}

				if(quad == 7)
				{
					vector<double> a(2);
					a[0] = kdtree[center].getDivisionX();
					a[1] = kdtree[center].getDivisionY();
					int Q = kdtree[current].getQuadrant(a);
					// Recurse along XN,YN, XO
					kid = kdtree[current].getChildValue(Q);
					if(kid != -1)
						findOverlaps(pop, cm,kid, center, radius, list);

					Q = Q ^ 1;

					kid = kdtree[current].getChildValue(Q);
					if(kid != -1)
						findOverlaps(pop,cm, kid, center, radius, list);
				//	findOverlaps(kdtree[parent].getChildValue(Q), center, radius, list);
				}
				
				if(quad == 13)
				{
					vector<double> a(2);
					a[0] = kdtree[center].getDivisionX();
					a[1] = kdtree[center].getDivisionY();
					int Q = kdtree[current].getQuadrant(a);
					//Recurse along XN, YN, YO
					kid = kdtree[current].getChildValue(Q);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);
					Q = Q ^ 2;
					kid = kdtree[current].getChildValue(Q);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);
				}

				if(quad == 5)
				{
					vector<double> a(2);
					a[0] = kdtree[center].getDivisionX();
					a[1] = kdtree[center].getDivisionY();
					int Q = kdtree[current].getQuadrant(a);
					// Recurse down XN, YN
					kid = kdtree[current].getChildValue(Q);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);
				}

                                double t=pop.collisionTime(kdtree[current].getPopIndex(),kdtree[center].getPopIndex());
                                if(t>=0.0 && t<pop.getTimeStep()) {
                                        cm.addPotentialCollision(kdtree[current].getPopIndex(),kdtree[center].getPopIndex(),t);
                                }
					
			}
			//else if it is a leaf
			else if( (current >= 0) && (center != current))
			{
                                double t=pop.collisionTime(kdtree[current].getPopIndex(),kdtree[center].getPopIndex());
                                if(t>=0.0 && t<pop.getTimeStep()) {
                                        cm.addPotentialCollision(kdtree[current].getPopIndex(),kdtree[center].getPopIndex(),t);
                                }
			}
		}

		bool withinRadius(int c1, int c2, double radius)
		{
			double radiusSquared = radius * radius;
			float p1[2];
			float p2[2];
			// THIS STUFF USES DIVISIONS!
			p1[0] = kdtree[c1].getDivisionX();
			p1[1] = kdtree[c1].getDivisionY();
			p2[0] = kdtree[c2].getDivisionX();
			p2[1] = kdtree[c2].getDivisionY();
			float distSquared = pow( (p1[0] - p2[0]), 2) + pow( (p1[1] - p2[1]), 2);
			distSquared = distSquared * distSquared;
 
			if(radiusSquared >= distSquared)
				return true;
			else
				return false;
		}

		//  Done
		template<class Population>
		int findNode(Population &pop, int popIndex, int current)
		{
			return popIndex;
			////cout << "************** IN FINDNODE ********************" << endl;
			////cout << "current = " << current << endl;
			////cout << "current's index= " << kdtree[current].getPopIndex() << endl;
			////cout << "popIndex= " << popIndex << endl; }
		//	if(current != -1)
		//	{
				if(kdtree[current].getPopIndex() == popIndex)
				{
		//			////////cout << "******************* RETURNING current: " << current << "****************************" << endl;
					return current;
				}
				
				vector<double> p;
				p.resize(2);
	
				p[0] = pop.getx(popIndex);
					p[1] = pop.gety(popIndex);
				
				////cout << "(x,y)= (" << p[0] << "," << p[1] << ")" << endl;
				int Q = kdtree[current].getQuadrant(p);
			
				////cout << "Q= " << Q << endl;
				////cout << "***************** RECURSING DOWN FINDNODE **********\n" << endl;

				return findNode(pop, popIndex, kdtree[current].getChildValue(Q));
		//	}
		//	else return -1;
		}
		
		int getFreeCell(void)
		{
			int tmp = FreeNode;	
			FreeNode = FreeNode + 1;
			if(FreeNode >= kdtree.size()) 
			{
				FreeNode = m_iRoot;
				//////cout << "Wrapping around" << endl;
			}
			//for(int i = tmp + 1; i < kdtree.size(); i++)
			//	if(kdtree[i].isEmpty() )
			//	{
			//		FreeNode = i;
			//		break;
			//	}
			return tmp;
		};

		int builds;
		double radius;
		int FreeNode;
		vector<Node> kdtree;
		int m_iRoot;

		// debugging variable
		int stop;
};

#endif

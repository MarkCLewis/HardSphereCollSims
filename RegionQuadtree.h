#ifndef REGIONQUADTREE
#define REGIONQUADTREE

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
			////////cout << "Creating a node" << endl;
			m_iPopIndex = -1;
			m_iParentIndex = -1;

			m_vChildren.resize(4);

			m_vChildren[0] = -1;
			m_vChildren[1] = -1;
			m_vChildren[2] = -1;
			m_vChildren[3] = -1;

			divX = 0;
			////////cout << divX << endl;	
			divY = 0;
		}

		
		int getQuadrant(vector<double> &point) {
			int ret=0;
		//	cout << "point[0] = " << point[0] << endl;
		///	cout << "point[1] = " << point[1] << endl;
		//	cout << "*******************" << endl;
		//	cout << "pop = " << m_iPopIndex << endl;
		//	cout << "ret = " << ret << endl;
		//	cout << "divX = " << divX << endl;
			//cout << "ret = " << ret << endl;
			if(point[0]>divX) { ret|=1; }
			////////cout << "point[1] = " << point[1] << endl;
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
	
class RegionQuadTree {

	public:
		RegionQuadTree()
		{
			builds = 0;
			radius = 0;
			stop = 0;
		}
		

	// Build a RegionQuadTree representation of the population
	// This will call virtually all the functions in this
	// class.
	template<class Population>
	void build(Population &SamplePop)
	{
		////////cout << "Building tree" << endl;
//		////////cout << "root = " << m_iRoot << endl;
//		////////cout << "qtree.size = " << qtree.size() << endl;
	        int num = SamplePop.getNumBodies();
		qtree.clear();
	        qtree.resize( (num * 2) );
		m_iRoot = 0;
		FreeNode = 0;
		// radius = 2 * particleRadius + 3*RMS*dt
		//////cout << "Builds= " << builds << endl;
		if( (builds%100) == 0 )
			radius = (2 * SamplePop.getMaxParticleRadius() ) + (3 * getRMSVelocity(SamplePop) * SamplePop.getTimeStep() ) ;
		builds = builds + 1;
		////////cout << "Built an array of " << (num * 2) << " elements"<< endl;
		////////cout << "There are " << num << " particles. " << endl;
		//////////cout << "root's divX = " << qtree[m_iRoot].getDivisionX() << endl;
		vector<double> coord(3);
	        for(int i = 0; i < num; i++)
	        {
			coord[0] = SamplePop.getx(i);
			coord[1] = SamplePop.gety(i);
			coord[2] = SamplePop.getz(i);
			//////cout <<"********************************" << endl;
	//		////////cout << "Adding a coordinate into tree " << endl;
			//////cout << "(" << coord[0] << "," << coord[1] << ")" << endl;
//			////////cout << "i= " << i << endl;
			insert(SamplePop, coord, i);
	        }
	};

	// radius = 2 * particleRadius + 3*RMS*dt
	// Done Coding.. need testing
	template<class Population>
	double getRMSVelocity(Population &pop) 
	{
		//////cout << "In RMS function" << endl;
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
					//////cout << "index= " << index << endl;
				}
			while( qtree[index].getPopIndex() != -1 );

			if(index == 0) index = 1;
			parent = qtree[index].getParentIndex();
			temp = pow((pop.getvx(parent) - pop.getvx(index)),2);
			temp += pow((pop.getvy(parent) - pop.getvy(index)),2);
			temp += pow((pop.getvz(parent) - pop.getvz(index)),2);
			vms = 0.5 * sqrt(temp);
			sum += vms;
			//////cout << "In loop, i = " << i << endl;
			
		}
	
		double average = sum/n;
	
		
		//////cout << "************************************" << endl;
		//////cout << "Average Velocity: " << average  << endl;
		//////cout << "************************************" << endl;
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
//		////////cout << "\n\n********ENTERING FINDINITIALPOTENTIALCOLLISIONS********" << endl;
		// Fix this so that this function does not go ahead and scan nodes which do not match
		// to particles in the simulation
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
		//	////////cout << "Qtree.size = " << qtree.size() << endl;
			recursiveInsert(pop, popIndex, point, m_iRoot);
		};

		// Recursive Insertation of data into the RegionQuadTree
		// I might later replace this with an iterative insertation
		// method to get a slight speed increase, but we'll see.
		// For testing purposes this will do.
		template<class Population>
		bool recursiveInsert(Population &pop, int &popIndex, vector<double> &point, int parent)
		{
		//	////////cout << "root=" << m_iRoot << endl;
			int Q = qtree[parent].getQuadrant(point);
		//	////////cout << "Quadrant is " << Q << endl;

			if(qtree[parent].getChildValue(Q) != -1)
			{
				recursiveInsert(pop, popIndex, point, qtree[parent].getChildValue(Q) );
			}
			else
			{
//				////////cout << "Adding element " << endl;
				int index = getFreeCell();
				if(parent != index)
				{
		//		////////cout << "index = " << index << endl;
		//		////////cout << "parent = " << parent << endl;
//					if(parent == 0)
////					{
						//////cout << " | Q = " << Q ;
						//////cout << " | index = " << index ; 
						//////cout << " | popIndex = " << popIndex;
						//////cout << " | parent = " << parent << endl;
//					}
					qtree[parent].setChildValue(Q,index);
					qtree[index].setPopIndex(pop,popIndex);
					qtree[index].setParentIndex(parent);
		//			////////cout << "Done adding element" << endl;
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
		//	////////cout << "Warning!!! Somehow the RegionQuadTree vector is full!!!" << endl;
			return false;
		};


		template<class TreeForcing, class Population>
		vector<int> getCollisions(Population &pop, TreeForcing &cm, int popIndex, double radius)
		{
			//Find the Node
//			////////cout << " IN GETCOLLISIONS" << endl;
//			////////cout << "popIndex= " << popIndex << endl; 
			int index = findNode(pop, popIndex, m_iRoot);
			vector<int> a;
	//		////////cout << "parent = " << m_iRoot << endl;
	//		////////cout << "center = " << index << endl;
	//		////////cout << "radius = " << radius << endl;
	//		////////cout << "root's kids = " << endl;
	//		////////cout <<"	* " << qtree[m_iRoot].getChildValue(0) << endl;
	//		////////cout <<"	* " << qtree[m_iRoot].getChildValue(1) << endl;
	//		////////cout <<"	* " << qtree[m_iRoot].getChildValue(2) << endl;
	//		////////cout <<"	* " << qtree[m_iRoot].getChildValue(3) << endl;
			findOverlaps(pop,cm,m_iRoot, index, radius, a);
			////////cout << "a.size= " << a.size() << endl;
			////////cout << "continue?" << endl; cin >> stop;
			return a;
		}
	
		// UNDER CONSTRUCTION
		template<class TreeForcing,class Population>
		void findOverlaps(Population &pop, TreeForcing &cm, int current, int center, double radius, vector<int> & list)
		{
			//cout << "************ENTERING FINDOVERLAPS FUNCTION*************" << endl;
			//cout << "current =" << current << endl;
			//cout << "center =" << center << endl;
			//cout << "is current a leaf? " << qtree[current].isLeaf() << endl;
			if(current >= 0 && (qtree[current].isLeaf() == false) )
			{
				double diffSquared, rSquared;
				int quad = 0;
				rSquared = radius * radius;

				diffSquared = qtree[current].getDivisionX() - qtree[center].getDivisionX();
				diffSquared = diffSquared * diffSquared;

				//////cout << "diffSquaredX = "  << diffSquared << endl;
				//////cout << "rSquaredX = "  << rSquared << endl;
				if(diffSquared <= rSquared) // 0011
					quad = quad ^ 3;
				else // 0001
					quad = quad ^ 1;
				
				diffSquared = qtree[current].getDivisionY() - qtree[center].getDivisionY();
				diffSquared = diffSquared * diffSquared;
				//////cout << "diffSquaredY = "  << diffSquared << endl;
				//////cout << "rSquaredY = "  << rSquared << endl;
				
				if(diffSquared <= rSquared) // 1100
					quad = quad ^ 12;
				else // 0100
					quad = quad ^ 4;
				
				//cout << "quad = " << quad << endl;
				int kid;
				if(quad == 15)
				{
					// recurse down all of children
					kid = qtree[current].getChildValue(0);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);

					kid = qtree[current].getChildValue(1);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);

					kid = qtree[current].getChildValue(2);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);

					kid = qtree[current].getChildValue(3);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);
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
						findOverlaps(pop, cm,kid, center, radius, list);

					Q = Q ^ 1;

					kid = qtree[current].getChildValue(Q);
					if(kid != -1)
						findOverlaps(pop,cm, kid, center, radius, list);
				//	findOverlaps(qtree[parent].getChildValue(Q), center, radius, list);
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
						findOverlaps(pop,cm,kid, center, radius, list);
					Q = Q ^ 2;
					kid = qtree[current].getChildValue(Q);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);
				}

				if(quad == 5)
				{
					vector<double> a(2);
					a[0] = qtree[center].getDivisionX();
					a[1] = qtree[center].getDivisionY();
					int Q = qtree[current].getQuadrant(a);
					// Recurse down XN, YN
					kid = qtree[current].getChildValue(Q);
					if(kid != -1)
						findOverlaps(pop,cm,kid, center, radius, list);
				}

                                double t=pop.collisionTime(qtree[current].getPopIndex(),qtree[center].getPopIndex());
				////cout << "t = " << t << endl;
				////cout << "pop.getTimeStep = " << pop.getTimeStep() << endl;
                                if(t>=0.0 && t<pop.getTimeStep()) {
                                        cm.addPotentialCollision(qtree[current].getPopIndex(),qtree[center].getPopIndex(),t);
                                }
					
	//			if(withinRadius(center,parent,radius))
	//				if(center != parent)
	//				list.push_back(parent);
			}
			//else if it is a leaf
			else if( (current >= 0) && (center != current))
			{
				//////cout << "********************************************" << endl;
				//////cout << "CHECKING IF " << parent << " IS NEAR " << center << endl;
				//////cout << "Radius: " << radius << endl;
				
//				if(withinRadius(center,parent,radius))
//					list.push_back(parent);
                                double t=pop.collisionTime(qtree[current].getPopIndex(),qtree[center].getPopIndex());
                                if(t>=0.0 && t<pop.getTimeStep()) {
                                        cm.addPotentialCollision(qtree[current].getPopIndex(),qtree[center].getPopIndex(),t);
				////cout << "t = " << t << endl;
				////cout << "pop.getTimeStep = " << pop.getTimeStep() << endl;
                                }
			}
		}

		bool withinRadius(int c1, int c2, double radius)
		{
			double radiusSquared = radius * radius;
			float p1[2];
			float p2[2];
			// THIS STUFF USES DIVISIONS!
			p1[0] = qtree[c1].getDivisionX();
			p1[1] = qtree[c1].getDivisionY();
			p2[0] = qtree[c2].getDivisionX();
			p2[1] = qtree[c2].getDivisionY();
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
			////cout << "current's index= " << qtree[current].getPopIndex() << endl;
			////cout << "popIndex= " << popIndex << endl; }
		//	if(current != -1)
		//	{
				if(qtree[current].getPopIndex() == popIndex)
				{
		//			////////cout << "******************* RETURNING current: " << current << "****************************" << endl;
					return current;
				}
				
				vector<double> p;
				p.resize(2);
	
				p[0] = pop.getx(popIndex);
					p[1] = pop.gety(popIndex);
				
				////cout << "(x,y)= (" << p[0] << "," << p[1] << ")" << endl;
				int Q = qtree[current].getQuadrant(p);
			
				////cout << "Q= " << Q << endl;
				////cout << "***************** RECURSING DOWN FINDNODE **********\n" << endl;

				return findNode(pop, popIndex, qtree[current].getChildValue(Q));
		//	}
		//	else return -1;
		}
		
		int getFreeCell(void)
		{
			int tmp = FreeNode;	
			FreeNode = FreeNode + 1;
			if(FreeNode >= qtree.size()) 
			{
				FreeNode = m_iRoot;
				//////cout << "Wrapping around" << endl;
			}
			//for(int i = tmp + 1; i < qtree.size(); i++)
			//	if(qtree[i].isEmpty() )
			//	{
			//		FreeNode = i;
			//		break;
			//	}
			return tmp;
		};

		int builds;
		double radius;
		int FreeNode;
		vector<Node> qtree;
		int m_iRoot;

		// debugging variable
		int stop;
};

#endif

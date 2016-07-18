#ifndef elements_H
#define elements_H
#include <vector>
#include "parameter.h"
#include "frandom.h"
using namespace std;
// to describe the Edge for routing

#ifdef _edge_with_bandwidth
struct Edge 
{
  // int Source;           // First vertex in edge is implicit
    int Dest;             // Second vertex in edge
	float Cost;
	float Bandwidth;
	int Addr;
	Edge(const Edge & Rhs){
		Dest=Rhs.Dest;Addr=Rhs.Addr;Cost=Rhs.Cost;Bandwidth=Rhs.Bandwidth;
	}
	Edge(){}  
	Edge(const int &dest,const float& cost,const float &bandw,const int &addr):
		Dest(dest),Cost(cost),Bandwidth(bandw),Addr(addr){}  
	void Assign( const int &dest, const float& cost,const float &bandw,const int &addr){
		Dest=dest;Cost=cost;Bandwidth=bandw;Addr=addr;
	}
};
#else
struct Edge 
{
  // int Source;           // First vertex in edge is implicit
    int Dest;             // Second vertex in edge
	float Cost;
	int Addr;
	Edge(const Edge & Rhs){	
		Dest=Rhs.Dest;Addr=Rhs.Addr;Cost=Rhs.Cost;
	}

	Edge(){}  
	Edge(const int &dest,const float& cost,const int &addr):
		Dest(dest),Cost(cost),Addr(addr){}  
	void Assign( const int &dest, const float& cost,const int &addr){
		Dest=dest;Cost=cost;Addr=addr;
	}
};


#endif
struct Node
{
	//int Name;  // Real name
	bool Type;	// 0: server, 1:switch
	int Slot;	// number of VM slots

	//routing variables:
	int Degree;	
	float Dist;        // Cost (after running algorithm)
    int Prev;             // Previous vertex on shortest link
    bool Visited;          // Extra variable for use in algorithm

	Edge* Adj;      // Adjacent vertices
	Node(){}
	Node(const int &node_name);
};
/* tree node */
struct tNode
{
	//int Name;  // Real name
	bool Type;	// 0: server, 1:switch
	int Slot;	// number of VM slots
	int Degree;	
	int Parent;
	Edge* Adj;      // Adjacent vertices
	tNode(){}
	tNode(const int &node_name);
};
// a pair of source and destination
struct Pair
{
	int s;
	int d;
	Pair(int s1,int d1):s(s1),d(d1){}
	Pair(){}
};
// a comparable object used in heap
class Comparable 
{
public:	
    int index;        // W
    float value;   // D(W)
	Comparable(){index=0;value=0;}
    Comparable(const int &idx, const float& va) : index( idx ), value( va ) { }  
	bool operator<( const Comparable & Rhs ) const
        { return value< Rhs.value; }
	const Comparable& operator=( const Comparable & Rhs )
	{
		index=Rhs.index;
		value=Rhs.value;
		return *this;
	}
};

// to define the Virutal data center
struct Cluster
{	
	int N;
	vector<float>B;
	float Arrivaltime;
	float Holdtime;
	Cluster(){}
	Cluster(const Cluster& Rhs);//deep copy
	void random();
	void random(float load);
	void random(int nVM,float minBw,float maxBw);
	void random(int nVM1,int nVM2,float minBw,float maxBw);
};
// to store the embedding solution
class Solution 
{
public:
	float Departtime; 
	int Slot[nServer];
	float Bandwidth[Ne];      
	Solution();
	void Clear();
	Solution(const Solution &Rhs);  

	bool operator<( const Solution & Rhs ) const
        { return Departtime< Rhs.Departtime; }
	
	const Solution& operator=( const Solution & Rhs );
};

struct Performance
{
	float sucess_rate;
	float max_utilization;
	float RC;
	Performance(){}
	void clear();
};
struct allocateRange
{
	//int ARs[maxN][maxSlot][maxSlot][nServer];
	//int ARx[maxN][maxN][maxN][Nv-nServer];
	bool ***p;
	//int **split;
	int n;
	int a;
	allocateRange(){}
	void resize(int nVM,int nSlot);
	bool& operator () (int x, int y, int i);
	void release();
	~allocateRange();
};

#endif
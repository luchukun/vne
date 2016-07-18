#ifndef wavepath_H
#define wavepath_H
#include <limits>      // Used for infinity (for integers)
using namespace std;
const float Inf=numeric_limits<float>::infinity();
const int Ne=21;//Num of edges
const int Nv=14;// Num of nodes
const int Nw=32;//Num of wavelength per link
const int Nd=5;//max degree
// to maintain the wavelengths in the links
class Link 
{
public:	
	float HoldTime[Nw];	
	int Node1;
	int Node2;
	Link(const Link &x);
	Link();
//	~Link();
};

class Clink
{	public:
	bool Avai[Nw];
	Clink(){}
	Clink(const Clink &xin);
	void operator =(Link &xin);
	void operator =(bool c);
	void operator =(const Clink &xin);
	Clink operator &(const Clink &xin);
	int MinWave();
	bool IsFull();
};
// to describe the Clink//////////////
struct Edge 
{
  // int Source;                      // First vertex in edge is implicit
    int Dest;             // Second vertex in edge
	int Addr;
	float Cost;
	Edge(const Edge & E);
	Edge(){}  
	void Init( const int &d,const int &a,float c);
};

/*
class WCR
{
public:
	int InWave;
	bool OutWave[Nw];
	bool copy;
	WCR(){}
	Init(const int & i,const bool & o,const bool & c);
};
WCR::Init(const int & i,const bool & o,const bool & c)
{
	InWave=i,OutWave=o;copy=c;
}
*/
class Node
{
public:
	int Name;        // Real name
	bool MC;
	//bool Avai;	
	int Degree;	
	float Dist;        // Cost (after running algorithm)
    int Prev;             // Previous vertex on shortest link
    bool Mark;          // Extra variable for use in algorithm
	//bool OnTree;
	//bool IsBranch;
	bool OnTree;
	bool OnWave;	
	int nBranch;
	int Wavelength;
	Edge Adj[Nd];      // Adjacent vertices
	int WaveRange[Nw];
	Node(){}
	Node(const int &node_name);
	Clink WaveList;
	void Init(const int &node_name);

//	~Node();
};

// for routing assignment
struct Path 
{
    int Vdest;        // W
    float Distance;   // D(W)
    Path( int D = 0, float C = 0 ) : Vdest( D ), Distance( C ) { }  
	int operator<( const Path & Rhs ) const
        { return Distance < Rhs.Distance; }
};

class MidNode
{
public:
	int inode;
	int iwave;
	MidNode(int n, int w){inode=n; iwave=w;}
	//MidNode(const MidNode &x){inode=x.inode; iwave=x.iwave;}
};
class request
{	public:
	bool Type;//unicast or multicast
	int Snode;
	bool Dnodes[Nv];
	bool Member[Nv];
	int MemNumber;
	int Mem_accept;
	float HoldTime;
	int ConvTimes;
	request();
};
enum sim_type{fWC_down2up,fWC_up2down,MWC_down2up,sWC_up2down,sWC_up2down_v1,MWC_sparse};

#endif
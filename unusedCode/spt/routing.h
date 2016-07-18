#ifndef routing_H
#define routing_H
#include <iostream>      // for I/O
#include <strstream>
#include <fstream>       // for file I/O
#include <string>        // for strings
#include <Queue>
#include <list>
#include <vector>
//#include <Hash>
#include "Heap.h"
#include "wavepath.h"
using namespace std;

///////////////////////////////
class Graph
{
public:
	Link AllLink[Ne];
	request Req; 
	Node Table[Nv];//store the index of its Adj
	Graph();
	void AddEdge(const int & Source,const int & Dest, const float &Cost,const int &N_link);
	int Gen_gragh();
	void ClearTable(int type);
	void InitTable( );
	int Dijkstra();
	int Member_only( );//muticast routing algorithm
	float  WAA_down2up_MWCn();
	float  WAA_down2up_Sparse_MWCn();
	float  WAA_down2up();
	void oneWAA(int Snode);
	void oneWAA_v1(int Snode);
	void anyWAA(int Snode);
	float  WAA_up2down(sim_type WAA_flag);
	bool AssignWave(int &Vm,int &Vs,queue<MidNode> &RevNode);
	void PrintPathRec(int DestNode) const;
	void PrintPath() const;
	void LinkUpdate(float & arrive_time);
	float ProcessRequest(int SimTime,float  load,float  mem_ratio,sim_type WAA_flag);	
	int ProcessRequest(int  SimTime,float  load,float mul_ratio,float &unicast_bp,float &multicast_bp);
	int ProcessRequest();
	int UnicastRequest(const int SimTime,float load,float &unicast_bp);
	int MulticastRequest(const int SimTime,float load,float &multicast_bp);
	void WaveBlock(int &Vm,queue<MidNode> &RevNode);
	void LinkRandom(const float & link_ratio);
	void LinkClear();
	void MCRandom(const float & MC_ratio);
	void MCMax(int *NodePriority);
	void fMC();
	//~Graph();
};


#endif
#ifndef trees_H
#define trees_H
#include <iostream>      // for I/O
#include <strstream>
#include <fstream>       // for file I/O
#include <string>        // for strings
#include <queue>
#include <list>
#include <vector>
#include <stack>
#include <algorithm>
#include <functional>      // For greater<int>( )
#include "binaryheaps.h"
#include "elements.h"
#include "parameter.h"
using namespace std;
///////////////////////////////


class Tree
{
private:
	int root;
	tNode Table[Nv];//store the index of its Adj
public:
	
	float resBandwidth[Ne];
	float Bandwidth[Ne];

	//int ARs[maxN][maxSlot][maxSlot][nServer];
	//int ARx[maxN][maxN][maxN][Nv-nServer];
	bool cut[Ne][nServer];

	Tree();
	~Tree();
	void drawTree();
	void AddEdge(const int & s, const int & d, const float &cost,const float &bandw,const int&addr);
	void ClearNetwork();
	void RandomNetwork(float p,float p_minResBw,float p_maxResBw);

	void reserveResource(const Solution &map);
	void releaseResource(const Solution &map);
	float ProcessRequest(char algorithm,int numOfreq,float load,
						float &max_utilization,float &success_rate,float &bandwidth_cost,float &RC);
	float SingleRequest(char algorithm,int numOfreq,
		float prob,float p_minResBw,float p_maxResBw,int  numOfVM,float minBw,float maxBw,
		float &max_utilization,float &success_rate,float &bandwidth_cost);

	int findBottleneck(float* hostBw,const Solution &map);
	int findBottleneck(float *hostBw,float *TS,float *TD);
	void updateBottlenecks(float *hostBw,const Solution&map,float *TS,float *TD);

	bool CongestDetect(int s,Cluster& req,float* hostBw,float sumB,Solution& map);
	bool QuickFail(const Cluster& req,float& sumB,float* res_port_B);
	// the embedding algorithm
	bool Pertubation(Cluster &req,Solution &map);
	bool FirstFit(Cluster &req,Solution &map);
	bool BackTracking(Cluster &req,Solution &map);
	bool recursivePlacement(Cluster &req, Solution &map,vector<int>& assignment,float sumB,float* res_port_B,
						float* hostBw,int &numVMembedded);
};

#endif
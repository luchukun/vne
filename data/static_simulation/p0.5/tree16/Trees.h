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
#include "matrix.h"
using namespace std;
///////////////////////////////

class Tree
{
private:
	int root;
	tNode Table[Nv];//store the index of its Adj
	//int ARs[maxN][maxSlot][maxSlot][nServer];
	//int ARx[maxN][maxN][maxN][Nv-nServer];
	allocateRange AR[Nv];
	Matrix<int> split;
public:
	
	float resBandwidth[Ne];
	float Bandwidth[Ne];	
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
	bool NextFit(Cluster &req,Solution &map);
	bool BestFit(Cluster &req,Solution &map);
	bool BackTracking(Cluster &req,Solution &map);
	bool recursivePlacement(Cluster &req, Solution &map,vector<int>& assignment,float sumB,float* res_port_B,
						float* hostBw,int &numVMembedded);
	int AddBasic(int x,int y, int v,int v_exclude);
	void calculateAR(Cluster &req,float& sumB, int x, int y, int v);
	bool subFind(Cluster &req,int x, int y, int v,int v_exclude,Matrix<bool> &Q,Solution &map);
	bool subAssign(Cluster &req,int x, int y, int v,Matrix<bool> &Q,Solution &map);
	bool GreedyARAllocation(Cluster &req,Solution &map);
	int getSplit(Cluster &req,float& sumB, int x, int y);
};

#endif
#ifndef graphs_H
#define graphs_H
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
#include <ctime>
#include "binaryheaps.h"
#include "elements.h"
#include "parameter.h"

using namespace std;
//#define _server2server
#define _edge_rank
//#define _SQRT_LB
#define _ECMP_nonzero
//#define _cut_check
//#define _BottleckXHostBw
#define _newBottleckCriteria
///////////////////////////////
class Graph
{
private:
	Node Table[Nv];//store the index of its Adj
	float resBandwidth[Ne];
	float Bandwidth[Ne];
	float f[nServer][nServer][Ne];//f[s][d][i][j]表示f-sd在e(i,j)上的分配比例
	bool b[nServer][nServer][Ne];
	Pair pass [Ne][nPair];//pass[i][j][p]表示p(s,d)经过e(i,j)
	int n_pair[Ne];//number of pairs that passes e(i,j);
	int opposite[Ne];
	int costMatrix[nServer][nServer];
	vector<vector<int>> serverCluster;//物理机分组 
	
public:
	
	int Kmax;
	T_Routing routingOption;
	unsigned short Kpaths[nServer][nServer][Kspt][Nv+1];//Kpaths[][][][0]存储路径长度
	unsigned short n_path[nServer][nServer];
	Graph();
	~Graph();
	void printCostMatrix();
	void printserverCluster();
	void drawVL2();
	void drawVL2(float oversub_ToR,float oversub_AGG);
	void drawFattree();
	void drawBcube();
	void drawL2VL2();
	void AddEdge(const int & s, const int & d, const float &cost,const float &bandw,const int&addr);
	int drawGragh(string &file_name);
	void ClearTable();
	void ClearNetwork();
	void RandomNetwork(float p,float p_minResBw,float p_maxResBw);
	void PrintPathRec(int DestNode) const;
	void PrintPath(int s,int d) const;
	void reserveResource(const Solution &map);
	void releaseResource(const Solution &map);

	float ProcessRequest(char algorithm,bool enLProuting,int numOfreq,float load,
						float &max_utilization,float &success_rate,float &bandwidth_cost,float &RC);
	float SingleRequest(char algorithm,bool enLProuting,int numOfreq,
			float p,float p_minResBw,float p_maxResBw,int numberOfGroup,int min_numberOfVM, int max_numberOfVM,float minBw,float maxBw,
			float& max_utilization,float& success_rate,float& bandwidth_cost);
	float SingleRequest(char algorithm,bool enLProuting,int numOfreq,
		float prob,float p_minResBw,float p_maxResBw,int  min_numOfVM,int  max_numOfVM,float minBw,float maxBw,
		float &max_utilization,float &success_rate,float &bandwidth_cost);
	bool oversubscribedPertubationVmplacement(OversubscriptionCluster &req,Solution& map,int maxLoop);
	float HomoRequest(char algorithm,bool enLProuting,int numOfreq,
		float prob,float p_minResBw,float p_maxResBw,int  numOfVM,float minBw,float maxBw,
		float &max_utilization,float &success_rate,float &bandwidth_cost);

	int Dijkstra(int s);
	float Dijkstra(int s, int d);
	bool calCostMatrix(int numberofCluster,int numberofVm);
	void YenKSP(int s,int d);	
	float LoadBalance(int s,int d);
	float ECMP(int s,int d);
	float MaxFlowRouting(int s,int d);
	int Sort_fsd(int e);//return the number of pairs that passes e[i][j]
	float maxTraffic(int e,float* hostBw);
	int findBottleneck(int e,float* hostBw);
	int findBottleneck(float *hostBw,const Solution&map);
	int findBottleneck(float *hostBw,float *TS,float *TD);
	void updateBottlenecks(float *hostBw,const Solution&map,float *TS,float *TD);
	bool CongestDetect(int x,const OversubscriptionCluster& req,float* hostBw,float (*sum_capacity)[nServer],Solution& map);
	bool CongestDetect(int x,const Cluster& req,float* hostBw,float (*sum_capacity)[nServer],Solution& map);
	float CalcMaxLinkUt(const Cluster& req,float* hostBw,Solution&map);//dual LP
	float assignBandwidth(Cluster& req,float* hostBw,Solution& map);//|E|*min_cost_flow LP
	float LPmaxTraffic(int e,float* hostBw);
	float LPmaxTrafficUnderValidtraffic(OversubscriptionCluster& req,int e,float * hostBw,vector<int>& assignment);
	float LPRouting(const Cluster& req,float* hostBw,Solution&map);//LP and calc map.bandwidth
	float OptimalRouting(const Cluster& req,float* hostBw,Solution&map);//LP
	bool QuickFail(Cluster& req,Solution&map,float(*sum_capacity)[nServer],float& sumB,float* res_port_B);
#ifdef _cut_check
	bool CutCheck(const Cluster& req,float& sumB,float* hostBw);
#endif
	bool oversubscribedCongestDetect(int x,OversubscriptionCluster& req,float* hostBw,float (*sum_capacity)[nServer],Solution& map,vector<int>& assignment);
	bool oversubscribedVmpalcement(OversubscriptionCluster& req,Solution& map,int maxLoop,vector<int>& assignment);
	bool GroupAllocate(bool enLProuting,vector<int>& servercluster,OversubscriptionCluster& req,Solution& map,vector<int>& groupassignment,float (*sum_capacity)[nServer]);
	// the embedding algorithm
	bool oversubscribedQuickFail(OversubscriptionCluster& req,Solution&map,float(*sum_capacity)[nServer],float& sumB,float* res_port_B);
	bool Pertubation(Cluster &req,bool enLProuting,Solution &map);
	bool oversubscribedFirstFit(OversubscriptionCluster &req, Solution &map);
	bool oversusbcribedVmpalcement(OversubscriptionCluster& req,Solution&map,int maxLoop,vector<int>& assignment);
	bool PertubationVmplacement(OversubscriptionCluster& req,Solution&,int maxpertubation,vector<int>& assignment);
	void findserver(int& thelink,int& serverfrom,int&serverto,Solution &map,float * hostBw);
	bool randomDrop(Cluster &req,bool enLProuting,Solution &map);
	bool LocalSearch(Cluster& req,bool enLProuting,Solution& map);
	bool FirstFit(Cluster &req,bool enLProuting,Solution &map);
	bool NextFit(Cluster &req,bool enLProuting,Solution &map);
	bool BestFit(Cluster &req,bool enLProuting,Solution &map);

	bool recursivePlacement(Cluster &req,bool enLProuting,int max_backtrack,int &n_backtrack,Solution &map,vector<int>& assignment,
		float sumB,float* res_port_B,float (*sum_capacity)[nServer],float* hostBw,int &numVMembedded);
	bool BackTracking(Cluster &req,bool enLProuting,Solution &map);
	bool recursivePlacement(Cluster &req,bool enLProuting,Solution &map,vector<int>& assignment,
		float sumB,float* res_port_B,float (*sum_capacity)[nServer],float* hostBw,int &numVMembedded);
	float VC_ACE(int N,float B,int star, float cost_factor,Solution &map);
	bool HVC_ACE(Cluster &req,Solution &map);

	//limited number of backtrack/pertubation version
	bool Pertubation(Cluster &req,bool enLProuting,int max_pertubation,Solution &map);
	bool BackTracking(Cluster &req,bool enLProuting,int max_backtrack,Solution &map);
	//test running time of RA
	float SingleRequest(char algorithm,bool enLProuting,int numOfreq,
		float prob,float p_minResBw,float p_maxResBw,int  numOfVM,float minBw,float maxBw,
		float &max_utilization,float &success_rate,float &bandwidth_cost,float &t_RA);
	bool Pertubation(Cluster& req,bool enLProuting,Solution& map,float& t_RA,int& n_RA);

#ifdef _enumeration
	bool ExhaustiveSearch(Cluster &req,bool enLProuting, Solution &map);
	bool Enumeration(Cluster &req, Solution &map,vector<int>& assignment,float sumB,float* res_port_B,
							   float* hostBw,int &numVMembedded);
#endif	
};


#endif
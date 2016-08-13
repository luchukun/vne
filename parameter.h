#ifndef parameter_H
#define parameter_H
#include <limits>      // Used for infinity (for integers)
using namespace std;
const float Infinity=numeric_limits<float>::infinity();
const float min_error=1e-5;
enum T_Routing{_KshortestLB,_ECMP};

#define _VL2

#ifdef _Bcube
const int H=4;
const int maxDegree=H;
const int n_server_port=2;
const int nServer=H*H;
const int nToR=H;
const int nAggregate=H;
const int nServerInRack=H;
const int Nv=nServer+nToR+nAggregate;
const int nSwitch=Nv-nServer;
const int Ne=nServer*4;
const float gbpsCommodity=1e3;

// DCN parameters
const int Kspt=4;//K shortest paths for K shortestLB routing 
const int Kwidest=4;
const int maxSlot=4;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=2;
const int maxN=10;
const float minB=100;  
const float maxB=700;
#endif
#ifdef _VL2
//const int DA=8;
//const int DI=4;
const int DA=4;
const int DI=4;
const int nToR=DA*DI/4;
const int nAggregate=DI;
const int nIntermediate=DA/2;
const int nServerInRack=4;
//const int nServerInRack=2;
const int nServer=nServerInRack*nToR;
const int Nv=nServer+nToR+nAggregate+nIntermediate;
const int nSwitch=Nv-nServer;
const int Ne=DA*DI+DA*DI+nServer*2;
const int maxDegree=4;
const int n_server_port=1;
const float oversubstription=1;
//const float oversub_ToR=1;
//const float oversub_AGG=0.5;

const float gbpsCommodity=1e3;
const float gbpsServer2ToR=gbpsCommodity;//1Gbps
//const float gbpsToR2Agg=gbpsServer2ToR*nServerInRack*oversub_ToR/2;//10Gbps
//const float gbpsAgg2Int=gbpsToR2Agg*(DA/2)*oversub_AGG/(DA/2);//10Gbps
const float gbpsToR2Agg=gbpsServer2ToR*nServerInRack*oversubstription/2;//10Gbps
const float gbpsAgg2Int=gbpsToR2Agg*oversubstription;//10Gbps
// DCN parameters
const int Kspt=8;//K shortest paths
const int Kwidest=1;
const int maxSlot=2;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=2;
const int maxN=10;
const float minB=100;  
const float maxB=700;
#endif

#ifdef _FatTree
const int H=4;
const int maxDegree=H;
const int n_server_port=1;
const int nServer=H*(H/2)*(H/2);
const int nToR=H*(H/2);
const int nAggregate=H*(H/2);
const int nServerInRack=H/2;
const int nCore=(H/2)*(H/2);
const int Nv=nServer+nToR+nAggregate+nCore;
const int nSwitch=Nv-nServer;
const int Ne=H*H*H*3/2;//H*H*H/2+(H/2)*(H/2)*H*2+(H/2)*(H/2)*H*2;
const float gbpsCommodity=1e3;
const float oversubstription=1;
// DCN parameters
const int Kspt=4;//K shortest paths
const int Kwidest=3;
const int maxSlot=4;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=2;
const int maxN=10;
const float minB=100;  
const float maxB=700;

#endif
#ifdef _Tree
const int H=2;
const int maxDegree=H+1;
const int n_server_port=1;
const int nToR=(maxDegree-1)*(maxDegree-1);
const int nAggregate=maxDegree-1;
const int nServerInRack=4;
const int nServer=(maxDegree-1)*(maxDegree-1)*nServerInRack;
const int nSwitch=1+(maxDegree-1)+(maxDegree-1)*(maxDegree-1);
const int Nv=nServer+nSwitch;
const int Ne=Nv-1;
const float gbpsCommodity=1e3;
const int nLayer=3;
const float oversubstription=0.5;//;
// DCN parameters
const int Kspt=4;//K shortest paths
const int Kwidest=2;
const int maxSlot=4;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=4;
const int maxN=10;
const float minB=100;  
const float maxB=400;
#endif

#ifdef _L2VL2
const int nToR=3;
const int nEdge=2;
const int nServerInRack=2;
const int nServer=nServerInRack*nToR;
const int Nv=nServer+nToR+nEdge;
const int nSwitch=Nv-nServer;
const int Ne=2*nToR*nEdge+nServer*2;
const int nEdgePort=nToR;
const int maxDegree=max(nServerInRack+nEdge,nEdgePort);
const int n_server_port=1;
const float gbpsCommodity=1e3;
const float gbpsServer2ToR=gbpsCommodity;//1Gbps
const float gbpsToR2Edge=gbpsCommodity;//10Gbps

// DCN parameters
const int Kspt=2;//K shortest paths
const int Kwidest=2;
const int maxSlot=2;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=4;
const int maxN=8;
const float minB=100;  
const float maxB=400;
#endif
const T_Routing default_routing=_KshortestLB;
const float randBw=0.5;
const int Nvdc=20; // number of VDCs allocated to the substrate DCN
const int nPair=nServer*(nServer-1);
const float kv=0.04/60; //$/min
const float kb=0.00016/(60*1e3);// $/MB/min
const bool en_limited_pertubation=false;
const bool en_limited_backtrack=false;
const int _max_pertubation=1;
const int _max_backtrack=1;
const bool en_RA_runningtime=1;
#endif
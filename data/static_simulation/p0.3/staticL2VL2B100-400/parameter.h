#ifndef parameter_H
#define parameter_H
#include <limits>      // Used for infinity (for integers)
using namespace std;
const float Infinity=numeric_limits<float>::infinity();
#define _L2VL2
#ifdef _VL2
const int DA=4;
const int DI=4;
const int nToR=DA*DI/4;
const int nAggregate=DI;
const int nIntermediate=DA/2;
const int nServerInRack=4;
const int nServer=nServerInRack*nToR;
const int Nv=nServer+nToR+nAggregate+nIntermediate;
const int nSwitch=Nv-nServer;
const int Ne=DA*DI+DA*DI+nServer*2;
const int maxDegree=4;
const int n_server_port=1;
const float gbpsServer2ToR=1e4;//1Gbps
const float gbpsToR2Agg=1e4;//10Gbps
const float gbpsAgg2Int=1e4;//10Gbps
// DCN parameters
const int K=4;//K shortest paths
const int Kmax=2;
const int maxSlot=2;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=4;
const int maxN=10;
const float minB=100;  
const float maxB=800;
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

// DCN parameters
const int K=4;//K shortest paths
const int Kmax=2;
const int maxSlot=2;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=4;
const int maxN=10;
const float minB=100;  
const float maxB=800;

#endif
#ifdef _Tree
const int H=2;
const int maxDegree=H+1;
const int n_server_port=1;
const int nToR=(maxDegree-1)*(maxDegree-1);
const int nAggregate=maxDegree-1;
const int nServerInRack=5;
const int nServer=(maxDegree-1)*(maxDegree-1)*nServerInRack;
const int nSwitch=1+(maxDegree-1)+(maxDegree-1)*(maxDegree-1);
const int Nv=nServer+nSwitch;
const int Ne=Nv-1;
const float gbpsCommodity=1e3;
const int nLayer=3;
const float oversubstription=0.5;
// DCN parameters
const int K=4;//K shortest paths
const int Kmax=2;
const int maxSlot=5;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=4;
const int maxN=10;
const float minB=100;  
const float maxB=400;
#endif

#ifdef _L2VL2
const int nToR=4;
const int nEdge=2;
const int nServerInRack=4;
const int nServer=nServerInRack*nToR;
const int Nv=nServer+nToR+nEdge;
const int nSwitch=Nv-nServer;
const int Ne=2*nToR*nEdge+nServer*2;
const int nEdgePort=4;
const int maxDegree=max(nServerInRack+2,nEdgePort);
const int n_server_port=1;
const float gbpsServer2ToR=1e3;//1Gbps
const float gbpsToR2Edge=1e3;//10Gbps

// DCN parameters
const int K=2;//K shortest paths
const int Kmax=2;
const int maxSlot=2;
// VDC parameters
const float muArrivaltime=2;//minute
const float muHoldtime=60; //minute
const int minN=4;
const int maxN=10;
const float minB=100;  
const float maxB=400;
#endif

const int Nvdc=20; // number of VDCs allocated to the substrate DCN
const int nPair=nServer*(nServer-1);
#endif
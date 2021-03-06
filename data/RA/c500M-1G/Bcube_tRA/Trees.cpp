#include <ilcplex/ilocplex.h>
#include<ctime>
#include "frandom.h"
#include "elements.h"
#include "Trees.h"
ILOSTLBEGIN
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVarCubic;
#define _server2server
#define _edge_rank
using namespace std;
Tree::Tree()
{
	for( int i = 0; i < Nv; i++ )
	{	
		Table[i].Type=1;
		Table[i].Slot=-1;
		Table[i].Parent=-1;
		Table[i].Degree=0;
	}

}
Tree::~Tree()
{
	for( int i = 0; i < Nv; i++ )
		delete []Table[i].Adj;
}

void Tree::drawTree()
{	
#ifdef _Tree
	float cost=1,line_rate=0;
	//allocate *Adj according to node degree
	// servers
	for (int i=0;i<nServer;i++){
		Table[i].Type=0;  Table[i].Slot=maxSlot;
		Table[i].Adj=NULL;//new Edge[n_server_port];
	}	
	// switches
	for (int i=0;i<nSwitch;i++){
		int switch_index=i+nServer;
		Table[switch_index].Type=1;	Table[switch_index].Slot=0;
		if (i<nToR)
			Table[switch_index].Adj=new Edge[nServerInRack];
		else
			Table[switch_index].Adj=new Edge[H];
	}	
	root=Nv-1;

	//layer 1
	for (int i=0;i<nToR;i++)
	{	int ToR_index=nServer+i;
		for (int j=0;j<nServerInRack;j++)
		{	int Server_index=nServerInRack*i+j;
			AddEdge(ToR_index,Server_index,cost,gbpsCommodity,Server_index);// Edge->Server
			Table[Server_index].Parent=ToR_index; // parent
		}
	}

	//layer 2
	line_rate=gbpsCommodity;//nServerInRack*gbpsCommodity*oversubstription;
	for (int i=0;i<nAggregate;i++)
	{	
		int Agg_index=nServer+nToR+i;
		for (int j=0;j<H;j++)
		{	
			int ToR_index=nServer+i*H+j;
			AddEdge(Agg_index,ToR_index,cost,line_rate,ToR_index);// Aggragate->Edge
			Table[ToR_index].Parent=Agg_index;	 // parent
			
		}
	}
	//layer 3
	//line_rate=line_rate*H*oversubstription;
	line_rate=gbpsCommodity;
	for (int i=0;i<nAggregate;i++)
	{   
		int Agg_index=nServer+nToR+i;
		AddEdge(root,Agg_index,cost,line_rate,Agg_index); // root->Aggragate
		Table[Agg_index].Parent=root; // parent

	}
	Table[root].Parent=-1;

	// cal the cut set
	int child=-1;
	for(int e=0;e<Ne;e++)
		for (int j=0;j<nServer;j++)
			cut[e][j]=0;
	//the servers  is in the subtree of its uplink
	for(int e=0;e<nServer;e++)
		cut[e][e]=1;
	//the server in the subtree of an edge switch
	for(int e=nServer;e<root;e++)// e: a switch and its uplink
	{	
		for (int j=0;j<nServer;j++)	for (int p=0;p<Table[e].Degree;p++){
			child=(Table[e].Adj)[p].Addr;
			if (child!=e)
				cut[e][j]=cut[e][j]||cut[child][j];
		}
	}
	

#endif
}

// Add the edge ( Source, Dest, cost, Bandwidth ) to the Tree
void Tree::AddEdge( const int & s, const int & d, const float &cost,const float &bandw,const int&addr)
{	
	int i;
	//	e(s,d)
	i=Table[s].Degree;	
	(Table[s].Adj)[i].Assign(d,cost,addr);
	Bandwidth[addr]=bandw;	resBandwidth[addr]=bandw;
	Table[s].Degree++;	
	//	e(d,s)
	/*
	i=Table[d].Degree;
	(Table[d].Adj)[i].Assign(s,cost,addr+Ne/2);
	Bandwidth[addr+Ne/2]=bandw;	resBandwidth[addr+Ne/2]=bandw;
	Table[d].Degree++;
	*/
	//opposite[addr]=addr+Ne/2;	opposite[addr+Ne/2]=addr;
}




void Tree::ClearNetwork()
{
	//renew server space
	for (int i=0;i<nServer;i++){
		Table[i].Slot=maxSlot;
	}
	// renew link capacity
	for (int e=0;e<Ne;e++){	
		resBandwidth[e]=Bandwidth[e];
	}		
}
void Tree::RandomNetwork(float p,float p_minResBw,float p_maxResBw)
{
	//renew server space
	for (int i=0;i<nServer;i++){
		if (rand_b01(p))
			Table[i].Slot=unif_int(0,maxSlot);
		else
			Table[i].Slot=maxSlot;
	}
	// renew link capacity
	for (int e=0;e<Ne;e++)
	{	
		if (rand_b01(p))
			resBandwidth[e]=unif_int(Bandwidth[e]*p_minResBw,Bandwidth[e]*p_maxResBw)/100*100;
		else
			resBandwidth[e]=Bandwidth[e];

	}
}
void Tree::reserveResource(const Solution&map)
{	
	// bandwidth reservation
    for (int e=0;e<Ne;e++)	
		resBandwidth[e]-=map.Bandwidth[e];
	// VM slot reservation
	for (int i=0;i<nServer;i++)
		Table[i].Slot-=map.Slot[i];
}

void Tree::releaseResource(const Solution&map)
{
	// bandwidth release
    for (int e=0;e<Ne;e++)	
		resBandwidth[e]+=map.Bandwidth[e];

	// VM slot release
	for (int i=0;i<nServer;i++)
		Table[i].Slot+=map.Slot[i];
}

float Tree::ProcessRequest(char algorithm,int numOfreq,float load,
						float &max_utilization,float &success_rate,float &bandwidth_cost,float &RC)
{
	ClearNetwork();
	init_genrand(0);
	//RandomNetwork(0.5);
	//init_genrand( (unsigned long)time( NULL ));				
	float t=0;

	Cluster req;	// embedding input
	Solution map;	//embedding results

	//performance variables:	
	int n_success=0,n_req=0;
	float revenue=0;
	success_rate=0,max_utilization=0,bandwidth_cost=0,RC=0;
	BinaryHeap<Solution> depart_list(Nvdc);
	while(n_req<numOfreq){
		req.random(load);
		n_req++;
		t=t+req.Arrivaltime;
		// process depart event before next arrival
		while (!depart_list.isEmpty())
		{
			map=depart_list.findMin();
			if (map.Departtime>t)
				break;
			releaseResource(map);
			depart_list.deleteMin();
		}
		bool accept;
		switch(algorithm)
		{
		case 'P':
			accept=Pertubation(req,map);
			break;
		case 'F':
			accept=FirstFit(req,map);
			break;
		case 'B':
			accept=BackTracking(req,map);
			break;
		default:
			break;
		}
		if (accept)
		{   
			float virtual_bandwidth=0,bw_used=0,max_ut=0;
			reserveResource(map);
			map.Departtime=t+req.Holdtime;
			depart_list.insert(map);
			//performance evaluation
			n_success++;
			
			for (int j=0;j<req.N;j++)
				virtual_bandwidth+=req.B[j];
			revenue+=(kv*req.N+kb*virtual_bandwidth)*req.Holdtime;
			for (int j=0;j<Ne;j++)
				bw_used+=map.Bandwidth[j];
			
			for (int e=0;e<Ne;e++)	
			{
				max_ut=max(max_ut,resBandwidth[e]/Bandwidth[e]);
				//max_ut+=(1-resBandwidth[e]/(Table[v].Adj)[p].Bandwidth);
			}
			bandwidth_cost+=bw_used;
			max_utilization+=max_ut;
		}	
	}
	success_rate=(float)n_success/(float)n_req;

	if(n_success>0){
		max_utilization=max_utilization/(float)n_success;
		bandwidth_cost=bandwidth_cost/(float)n_success;
		RC=revenue/t;//only revenue is considered
	}
    return success_rate;
}

float Tree::SingleRequest(char algorithm,int numOfreq,
		float prob,float p_minResBw,float p_maxResBw,int numOfVM,float minBw,float maxBw,
		float &max_utilization,float &success_rate,float &bandwidth_cost)
{	
	ClearNetwork();
	init_genrand(0);//init_genrand( (unsigned long)time( NULL ));	
	
	// memory allocation for ARs ARx
	if(algorithm=='A'){	
		for(int v=0;v<Nv;v++)
			if(Table[v].Type==0)
				AR[v].resize(numOfVM,maxSlot);//a=-1 it's a switch
			else
				AR[v].resize(numOfVM,-1);
		split.resize(numOfVM,numOfVM);
	}
	Cluster req;	//embedding input
	Solution map;	//embedding results
	//performance variables:	
	int n_success=0,n_req=0;
	success_rate=0,max_utilization=0,bandwidth_cost=0;
	//BinaryHeap<Solution> depart_list(Nvdc);
	while(n_req<numOfreq){
		RandomNetwork(prob,p_minResBw,p_maxResBw);
		req.random(numOfVM,minBw,maxBw);
		n_req++;
		bool accept;
		switch(algorithm)
		{
		case 'P':
			accept=Pertubation(req,map);
			break;
		case 'A':
			accept=GreedyARAllocation(req,map);
			break;
		case 'B':
			accept=BackTracking(req,map);
			break;	
		case 'F':
			accept=FirstFit(req,map);
			break;
		case 'G':
			accept=BestFit(req,map);
			break;
		case 'N':
			accept=NextFit(req,map);
			break;
		default:
			break;
		}
		if (accept)
		{   
			//performance evaluation
			n_success++;
			float max_ut=0,bw_used=0;
			for (int e=0;e<Ne;e++)	
			{
				if(resBandwidth[e]>0)
				max_ut=max(max_ut,map.Bandwidth[e]/resBandwidth[e]);
			}
			max_utilization+=max_ut;
			
			for (int e=0;e<Ne;e++){
				bw_used+=map.Bandwidth[e];
				
			}
			bandwidth_cost+=bw_used;
		}	
	}
	success_rate=(float)n_success/(float)n_req;
	if(n_success>0){
		max_utilization=max_utilization/(float)n_success;
		bandwidth_cost=bandwidth_cost/(float)n_success;
	}
	if(algorithm=='A'){	
		for(int v=0;v<nServer;v++)	
			AR[v].release();//a=-1 it's a switch
		for(int v=nServer;v<Nv;v++)	
			AR[v].release();//a=-1 it's a switch
	}
    return success_rate;
}


bool Tree::QuickFail(const Cluster& req,float& sumB,float* res_port_B)
{
	bool resource_satisfied=true;
	float res_sumB=0;
	int res_slot=0;
	for (int i=0;i<req.N;i++)// total suscribed bandwidth
		sumB+=req.B[i];
	for(int i=0;i<nServer;i++){
		res_port_B[i]+=resBandwidth[i];
		if (Table[i].Slot>0&&res_port_B[i]>0){
			res_slot+=Table[i].Slot;
			res_sumB+=res_port_B[i];
		}
	}
	if (res_slot<req.N||res_sumB<sumB)
		resource_satisfied=false;
	return resource_satisfied;
}
int Tree::findBottleneck(float *hostBw,const Solution& map)
{	
	int the_server=-1;
	float sum_hostBw=0;
	for (int i=0;i<nServer;i++)
		sum_hostBw+=hostBw[i];
	if (sum_hostBw==0)
		return the_server;
	bool selected[Ne]={0};
	for(int n_selected=0;n_selected<Ne;n_selected++)
	{
		float mlu=0;	int the_link=0;
		for (int e=0;e<Ne;e++){
			if (mlu<map.Bandwidth[e]&&selected[e]==false)
			{
				mlu=map.Bandwidth[e];
				the_link=e; //find the most congested link
			}
		}
		selected[the_link]=true;
		
		float max_traffic=0;
		for (int j=0;j<nServer;j++)
			if (cut[the_link][j]&&hostBw[j]> max_traffic){
				the_server=j;
				max_traffic=hostBw[j];
			}	
		
	}
	return the_server;
}

void Tree::updateBottlenecks(float *hostBw,const Solution&map,float *TS,float *TD)
{
	float mlu=0;	int the_link=0;
	for (int e=0;e<Ne;e++){
		if (mlu<map.Bandwidth[e]/resBandwidth[e])
		{
			mlu=map.Bandwidth[e]/resBandwidth[e];
			the_link=e; //find the most congested link
		}
	}

	int nserver_in_subtree=0,nserver_out_subtree=0;
	for (int j=0;j<nServer;j++){
		if (cut[the_link][j]&&hostBw[j]> 0){
			nserver_in_subtree++;
		}
		if (cut[the_link][j]==0&&hostBw[j]> 0){
			nserver_out_subtree++;
		}
	}
	for (int j=0;j<nServer;j++){
		if (cut[the_link][j]&&hostBw[j]> 0){
			TS[j]+=nserver_out_subtree;
		}
		if (cut[the_link][j]==0&&hostBw[j]> 0){
			TD[j]+=nserver_in_subtree;
		}
	}
}

int Tree::findBottleneck(float *hostBw,float *TS,float *TD)
{
	int the_server=-1;
	Pair bottleneck(-1,-1);
	float max_TS=-1,max_TD=-1;
	for(int x=0;x<nServer;x++)
	{	if (max_TS<TS[x]&&hostBw[x]>0)
		{	max_TS=TS[x];
			bottleneck.s=x;
		}
		if (max_TD<TD[x]&&hostBw[x]>0)
		{	max_TD=TD[x];
			bottleneck.d=x;
		}
	}

	if (rand_b01(0.5)){
		the_server=bottleneck.s;
		TS[the_server]=0;
	}
	else{
		the_server=bottleneck.d;
		TD[the_server]=0;
	}	
/*
	float max_TS=-1;
	for(int x=0;x<nServer;x++)
	{	if (max_TS<TS[x]+TD[x]&&hostBw[x]>0)
		{	max_TS=TS[x]+TD[x];
			the_server=x;
			
		}
		
	}
	if(the_server>0&&the_server<nServer){
		TS[the_server]=0;TD[the_server]=0;
	}*/
	return the_server;

}

bool Tree::CongestDetect(int s,Cluster& req,float* hostBw,float sumB,Solution& map)
{		
	//Solution oldmap(map);
	bool congest=false;
	
	//whether e(v,w) is congested?
	/*int parent=s,e; // the up_link edge of x
	while(1)
	{	e=parent;
		float B=0; // 
		for (int j=0;j<nServer;j++)
			if (cut[e][j]) B+=hostBw[j];
		map.Bandwidth[e]=min(B,sumB-B);

		if (map.Bandwidth[e]>resBandwidth[e]){
			congest=true;
			//map=oldmap;
			break;
		}
		parent=Table[parent].Parent;	
		if (parent==root) 
			break;
	}	//end for: congestion test for concurrent traffic
	*/

	int _rank[Ne];	float _load[Ne];bool caculated[Ne]={0};
	for( int e=0;e<Ne;e++)
	{	_rank[e]=e;				
		_load[e]=resBandwidth[e]-map.Bandwidth[e];				
	}
	for(int a=0;a<Ne-1;a++)
		for(int b=a;b<Ne;b++)
		{					
			if(_load[a] >_load[b]) //f-sd-e(i,j)
			{	//交换a,b位置
				swap(_load[a],_load[b]);
				swap(_rank[a],_rank[b]);						
			}
		}
	
	for (int j=0;j<Ne;j++)
	{	
		int e=_rank[j];
		if (caculated[e])
			continue;
		float B1=0,B2=0; // 
		for (int j=0;j<nServer;j++)
			if (cut[e][j]) 
				B1+=hostBw[j];
			else
				B2+=hostBw[j];
		
		map.Bandwidth[e]=min(B1,B2);
		caculated[e]=true; 
		if (map.Bandwidth[e]>resBandwidth[e])
		{	//congestion detected!
			congest=true;
			break;
		}		
		
	}	
	return congest;
}
// pertubation only when a VM can not placed
bool Tree::Pertubation(Cluster& req,Solution& map)
{	
	sort(req.B.begin(),req.B.end(),greater<float>());	// sort B1,...BN in descending order		
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	vector<vector<bool>> tabu(req.N,vector<bool>(nServer,false));
	//intiate the map variables:
	map.Clear();
	
	//quick fail 
	float res_port_B[nServer]={0},sumB=0;
	if (QuickFail(req,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!
	int numVMembedded=0,maxLoop=req.N*nServer;
	float TS[nServer]={0},TD[nServer]={0};
	for(int t=0;t<maxLoop+1;t++)
	{	
		// find VM(Bi) with max B, and try to place to server x
		int i,x;	//i: visit the VMs, x: visit the servers
		for (i=0;i<req.N;i++)
			if (assignment[i]==-1)
				break;
		//to place VM(Bi)into the servers, {0,1,...tabu[i]} are not considered		
		int n_usefull_server=0;
		for (x=0;x<nServer;x++)
		{	
			if (tabu[i][x])// this server x is in the tabu list of VM(Bi)
				continue;
			if (Table[x].Slot<=map.Slot[x])// there is no space in server m
				continue;	//try the next server for VM(Bi)
			if(res_port_B[x]<min(hostBw[x]+req.B[i],sumB-hostBw[x]-req.B[i]))
				continue;	 //less ingress bandwidth!!  		
			n_usefull_server++;
			map.Slot[x]+=1;hostBw[x]+=req.B[i];assignment[i]=x;numVMembedded++;
			if (CongestDetect(x,req,hostBw,sumB,map)){
				map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;	
				updateBottlenecks(hostBw,map,TS,TD);
			}
			else	
				break;
			//else, try server x+1
		}	//end for server x

		if (numVMembedded==req.N)	
			return true;
		if (x>=nServer)
		{
			if(n_usefull_server==0||numVMembedded==0)	 // running out of VM slots!
				return false;
			// fail due to conestion,try pertubation
			
			// fail due to conestion,try pertubation						
			int the_server=-1,the_vm;

			the_server=findBottleneck(hostBw,TS,TD);
			for(the_vm=req.N-1;the_vm>=0;the_vm--)
				if (assignment[the_vm]==the_server)
					break;	
			if (the_server<0||the_server>=nServer||the_vm<0||the_vm>req.N) 
				return false; //fail due to bisection traffic, rather than concurrent congestion.
			//n_bottle[the_server]=0;
			//unload a VM from the bottleneck server
			map.Slot[the_server]-=1;	hostBw[the_server]-=req.B[the_vm];
			assignment[the_vm]=-1;		numVMembedded--;
			tabu[the_vm][the_server]=1;


		}
		//else, Let's place the next VM(Bi+1)
	}
	return false;
}


// greedy pack a server with any fitable VMs
bool Tree::FirstFit(Cluster& req,Solution& map)
{	
	sort(req.B.begin(),req.B.end(),greater<float>());	// sort B1,...BN in descending order	
	
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in

	//intiate the map variables:
	map.Clear();
	
	//quick fail 
	float res_port_B[nServer]={0},sumB=0;
	if (QuickFail(req,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!

	//i: visit the VMs, j: visit the servers
	int numVMembedded=0,maxLoop=req.N*nServer;;
	for (int i=0;i<req.N;i++)
	{	int x;
		for (x=0;x<nServer;x++)
		{	// packing VMs( i=1:N )into server x
			
			if (Table[x].Slot<=map.Slot[x])// there is no space in server x
				continue;	//open the next server 
			if(res_port_B[x]<min(hostBw[x]+req.B[i],sumB-hostBw[x]-req.B[i]))
				continue;	 //less ingress bandwidth!!  			
			map.Slot[x]+=1;hostBw[x]+=req.B[i];assignment[i]=x;numVMembedded++;
			if (CongestDetect(x,req,hostBw,sumB,map))
			{	map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;	}
			else
			{
				if (numVMembedded==req.N)	
				{
					return true;
				}
				else
					break;
			}
			//else continue;	//to try the next server
					
		}	//end for server x
		if (x>=nServer)
			return false;
	}	// end for VM(Bi)

	
}


bool Tree::NextFit(Cluster& req,Solution& map)
{	
	sort(req.B.begin(),req.B.end(),greater<float>());	// sort B1,...BN in descending order		
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	//intiate the map variables:
	map.Clear();
	
	//quick fail 
	float res_port_B[nServer]={0},sumB=0;
	if (QuickFail(req,sumB,res_port_B)==false)
		return false;

	//i: visit the VMs, j: visit the servers
	int numVMembedded=0,lastHost=0;
	for (int i=0;i<req.N;i++)
	{	int x=lastHost,j=0;
		for (;j<nServer;j++)
		{	// packing VMs( i=1:N )into server x
			x=(x+j)%nServer;
			if (Table[x].Slot<=map.Slot[x])// there is no space in server x
				continue;	//open the next server 
			if(res_port_B[x]<min(hostBw[x]+req.B[i],sumB-hostBw[x]-req.B[i]))
				continue;	 //less ingress bandwidth!!  			
			map.Slot[x]+=1;hostBw[x]+=req.B[i];assignment[i]=x;numVMembedded++;
			if (CongestDetect(x,req,hostBw,sumB,map))
			{	//undo the placement and try the next server
				map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;
			}
			else{
				if (numVMembedded==req.N)	
				{	// already found the host for the last VM
					return true;
				}
				else // already found the host for the current VM
				{	lastHost=x;
					break;
				}
			}			
		}	//end for server x
		if (j>=nServer) // fail to find a suitable host
			return false;
	}	// end for VM(Bi)

	
}


bool Tree::BestFit(Cluster& req,Solution& map)
{	
	sort(req.B.begin(),req.B.end(),greater<float>());	// sort B1,...BN in descending order	
	float sum_capacity[nServer][nServer];
	
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	//intiate the map variables:
	map.Clear();
	
	//quick fail 
	float res_port_B[nServer]={0},sumB=0;
	if (QuickFail(req,sumB,res_port_B)==false)
		return false;
	//i: visit the VMs, j: visit the servers
	int numVMembedded=0;
	for (int i=0;i<req.N;i++)
	{	int x,bestHost=-1;
		float min_mlu=1;
		for (x=0;x<nServer;x++)
		{	// packing VMs( i=1:N )into server x
			
			if (Table[x].Slot<=map.Slot[x])// there is no space in server x
				continue;	//open the next server 
			if(res_port_B[x]<min(hostBw[x]+req.B[i],sumB-hostBw[x]-req.B[i]))
				continue;	 //less ingress bandwidth!!  			
			hostBw[x]+=req.B[i];
			CongestDetect(x,req,hostBw,sumB,map);
			float mlu=0;	
			for (int e=0;e<Ne;e++){
				if (mlu<map.Bandwidth[e]/resBandwidth[e]){
					mlu=map.Bandwidth[e]/resBandwidth[e];
				}
			}
			if (mlu<min_mlu) //better host is found!
			{
				min_mlu=mlu;
				bestHost=x;
			}
			//undo the placement and try the next server
			hostBw[x]-=req.B[i];

		}	//end for server x

		if (min_mlu>=1) // fail to find a suitable host
			return false;
		map.Slot[bestHost]+=1;hostBw[bestHost]+=req.B[i];assignment[i]=bestHost;numVMembedded++;

		if (numVMembedded==req.N)	
		{	// already found the host for the last VM
			return true;
		}

	}	// end for VM(Bi)	
}


// exhaustive search with load-balanced routing
bool Tree::BackTracking(Cluster &req, Solution &map)
{
	sort(req.B.begin(),req.B.end(),greater<float>());	// sort B1,...BN in descending order	
	float sum_capacity[nServer][nServer];
	
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in

	//vector<vector<bool>> tabu(req.N,vector<bool>(nServer,false));
	//intiate the map variables:
	map.Clear();
	
		//quick fail 
	float res_port_B[nServer]={0},sumB=0;
	if (QuickFail(req,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!

	//i: visit the VMs, j: visit the servers
	int numVMembedded=0;
	if (recursivePlacement(req, map,assignment,sumB,res_port_B,hostBw,numVMembedded))
	{
		return true;
	}
	else
		return false;

}
// used in backtracking algorihtm
bool Tree::recursivePlacement(Cluster &req, Solution &map,vector<int>& assignment,float sumB,float* res_port_B,
							   float* hostBw,int &numVMembedded)
{
// find VM(Bi) with max B, and try to place to server x
	int i,x;
	for (i=0;i<req.N;i++)
		if (assignment[i]==-1)
			break;
	//to place VM(Bi)into the servers, {0,1,...tabu[i]} are not considered
	for (x=0;x<nServer;x++)
	{	
		if (Table[x].Slot<=map.Slot[x])// there is no space in server m
			continue;	//try the next server for VM(Bi)
		if(res_port_B[x]<min(hostBw[x]+req.B[i],sumB-hostBw[x]-req.B[i]))
			continue;	 //less ingress bandwidth!! 
		map.Slot[x]+=1;hostBw[x]+=req.B[i];assignment[i]=x;numVMembedded++;
		if (CongestDetect(x,req,hostBw,sumB,map))
		{
			//routing failure,try the next server x+1 for VM(Bi)!!
			map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;	
		}
		else //not congestion
			if (numVMembedded==req.N)
			{
				return true;
			}
			else  //place the next VMs recursively
				if (recursivePlacement(req, map,assignment,sumB,res_port_B,
					hostBw,numVMembedded))
					return true;
				else // back tracking!!
				{	map.Slot[x]-=1;	hostBw[x]-=req.B[i];	assignment[i]=-1; numVMembedded--;
					continue;//try the next server x+1 for VM(Bi)!!
				}

	}//end for server x
	if (x>=nServer)	
	return false;	

}






// alogrithm based on allocation range
void Tree::calculateAR(Cluster &req, float& sumB, int x, int y, int v)
{
	//x=0,1,2,...N-1
	if(x>req.N-1||x<0||y>req.N-1||y<x||v>Nv-1||v<0)
		return;
	// suppose B is ordered ascendingly
	//reset ARs and ARx
	if(Table[v].Type==0) //server
	{	
		for(int n=0;n<=min(AR[v].a,y-x+1);n++)
			AR[v](x,y,n)=0;		
	}
	else //switch
	{
		for(int n=0;n<=y-x+1;n++)
			AR[v](x,y,n)=0;	
	}

	//layer 0: servers
	
	vector<bool> AR1(y-x+2,0);
	vector<bool> AR2(y-x+2,0);

	if(Table[v].Type==0) {//server
		for(int n=0;n<=min(Table[v].Slot,y-x+1);n++)
			AR1[n]=1;		
	}
	else{
		int es_low=-1,es_up=-1;
		for(int j=0;j<Nv;j++){
			if (Table[j].Parent==v)//it's a child
			{	//have not extended section
				if(split(x,y)==-1)
					continue;
				es_low=split(x,y);	
				es_up=-1;				
				for(int m1=es_low;m1<=y-x+1;m1++)
					if(AR[j](x,y,m1)==1){
						es_low=m1;//start of an ES
						for(int m2=m1;m2<=y-x+1;m2++){						
							if(AR[j](x,y,m2)==0){												
								m1=es_up;								
								break;//end of an ES!
							}
							else
								es_up=m2;
						}
						// for the found es
						if(es_up<y-x+1)
							es_up+=AddBasic(x+es_low,y,v,j);
						es_up=min(es_up,y-x+1);
						for(int k=es_low;k<=es_up;k++)AR1[k]=true;
					}
			}
		}
		es_low=0;
		es_up=AddBasic(x,y,v,-1);
		for(int k=es_low;k<=es_up;k++)AR1[k]=true;

	}
	// AR2 constrained by uplink
	int up=-1,low=-1;
	float c=0;
	if(v!=root){
		c=resBandwidth[v];		
	}
	else
		c=Infinity;

	float b=0;
	for(int j=0;j<=y-x;j++){
		b+=req.B[y-j];
		if(b>c)
			break;
		else
			up=j+1;
	}
	b=0;
	for(int j=0;j<=y-x;j++){
		b+=req.B[x+j];
		if(b>=sumB-c){			
			low=j+1;
			break;
		}
	}

	for(int n=0;n<=y-x+1;n++){
		if((up>=1&&n<=up)||(low>=1&&n>=low))
			AR2[n]=true;
	}

	
	// intersections
	AR[v](x,y,0)=1;
	if(Table[v].Type==0) {//server
		for(int n=0;n<=min(AR[v].a,y-x+1);n++)
			AR[v](x,y,n)=AR1[n]&AR2[n];		
	}
	else{
		for(int n=0;n<=y-x+1;n++)
			AR[v](x,y,n)=AR1[n]&AR2[n];
	}

}
int Tree::getSplit(Cluster &req,float& sumB, int x, int y)
{
	// split point
	float sumx2y=0;
	int split_point=0;
	for(int j=x;j<=y;j++)
		sumx2y+=req.B[j];
	if(sumx2y>sumB/2)
	{
		float b=0;
		for(int j=0;j<=y-x;j++){
			b+=req.B[y-j];
			if(b>sumB/2){								
				split_point=j+1;break;
			}
				
		}
		return split_point;
	}
	else
		return-1;
}
int Tree::AddBasic(int x, int y, int v,int v_exclude)
{
	int n_children=Table[v].Degree;
	vector<int>children(n_children,0);
	bool exclude[Nv]={false};
	if (v_exclude>=0)
		exclude[v_exclude]=true;

	// list of children j:index of children
	for(int i=0,j=0;i<Nv;i++){
		if (Table[i].Parent==v)//it's a child
		{	//have not extended section
			children[j]=i;
			j++;
		}
	}
	int x1=x,y1=y,up_added=0;
	for (int i=n_children;i>=1&&x1<=y1;i--) 
	{
		int max_up=0,max_child=-1,the_up=0;
		for(int j=0;j<n_children;j++){
			int child=children[j];
			if(exclude[child]==true) continue;
			int len=(Table[child].Type==0)?min(AR[child].a+1,y1-x1+2):y1-x1+2;
			for(int n=0;n<len;n++)
				if(AR[child](x1,y1,n)==0)
					break;
				else
					the_up=n;
			if(the_up<=0){
				exclude[child]=true;// basic section is null
				continue;
			}
			if(the_up>max_up){
				max_up=the_up;
				max_child=child;
			}
		}
		if(max_up>0&&max_child>=0&&max_child<Nv){
			up_added+=max_up;
			y1=y1-max_up;
			exclude[max_child]=true;
		}
		else
			break;
		
	}
	
	return up_added;


}

bool Tree::subAssign(Cluster &req,int x, int y, int v,Matrix<bool> &Q,Solution &map)
{
	int n=y-x+1,es_low=-1,es_up=-1;
	for(int j=0;j<Nv;j++){
		if (Table[j].Parent==v)//it's a child
		{	
			if(split(x,y)<=0)//have not extended section
				continue;
			es_low=split(x,y);	
			es_up=-1;				
			for(int m1=es_low;m1<=y-x+1;m1++){
				if(AR[j](x,y,m1)==1){
					es_low=m1;//start of an ES
					for(int m2=m1;m2<=y-x+1;m2++){							
						if(AR[j](x,y,m2)==0){												
							m1=es_up;break;//end of an ES!
						}
						else	
							es_up=m2;
					}
					// for the found es
					float b=0;
					if(n>es_up){
						if(subFind(req,x,y-es_up,v,j,Q,map)==true){							
							for(int m=y-es_up+1;m<=y;m++){
								Q(j,m)=1;	
								b+=req.B[m];
							}
							if(Table[j].Type==1)//it's a switch
								subAssign(req,y-es_up+1,y,j,Q,map);
							return true;
						}						
					}
					else
					{
						for(int m=x;m<=y;m++){
							Q(j,m)=1;	
							b+=req.B[m];
						}
						if(Table[j].Type==1)//it's a switch
							subAssign(req,x,y,j,Q,map);
						return true;
					}
					int uplink=j;
					if(uplink>=0&&uplink<Ne&&b>0){
						map.Bandwidth[uplink]=b;
					}
				}//end if start ES
			}//end for m1	
		}//end if j=it's a child
	}//end for j
	return subFind(req,x,y,v,-1,Q,map);
}

bool Tree::subFind(Cluster &req, int x, int y, int v, int v_exclude, Matrix<bool> &Q, Solution &map)
{
	int n_children=Table[v].Degree;
	vector<int>children(n_children,0);
	bool exclude[Nv]={false};
	if (v_exclude>=0)
		exclude[v_exclude]=true;

	// list of children j:index of children
	for(int i=0,j=0;i<Nv;i++){
		if (Table[i].Parent==v)//it's a child
		{	//have not extended section
			children[j]=i;
			j++;
		}
	}
	int x1=x,y1=y,vm_added=0;
	for (int i=n_children;i>=1&&x1<=y1;i--) 
	{
		int max_up=0,max_child=-1,the_up=0;
		for(int j=0;j<n_children;j++){
			int child=children[j];
			if(exclude[child]==true) continue;		
			int len=(Table[child].Type==0)?min(AR[child].a+1,y1-x1+2):y1-x1+2;
			for(int n=0;n<len;n++)
				if(AR[child](x1,y1,n)==0)
					break;
				else
					the_up=n;
			if(the_up<=0){
				exclude[child]=true;// basic section is null
				continue;
			}
			if(the_up>max_up){
				max_up=the_up;
				max_child=child;
			}
		}
		float b=0;
		if(max_up>0&&max_child>=0&&max_child<Nv){
			int vm_left=y1-x+1;
			if(max_up<=vm_left){
				vm_added+=max_up;
				for(int m=y1-max_up+1;m<=y1;m++){
					Q(max_child,m)=1;	
					b+=req.B[m];
				}
				if(Table[max_child].Type==1)//it's a switch
					subAssign(req,y1-max_up+1,y1,max_child,Q,map);

			}
			else
			{
				vm_added+=vm_left;
				for(int m=x1;m<=y1;m++){
					Q(max_child,m)=1;	
					b+=req.B[m];
				}
				if(Table[max_child].Type==1)//it's a switch
					subAssign(req,x1,y1,max_child,Q,map);
			}
			y1=y1-max_up;
			exclude[max_child]=true;
			int uplink=max_child;
			
			if(uplink>=0&&uplink<Ne&&b>0){
				map.Bandwidth[uplink]=b;
			}
			if(vm_added==y-x+1)
				return true;
		}
		else
			return false;		
	}	
	return false;
}
bool Tree::GreedyARAllocation(Cluster &req,Solution &map)
{
	sort(req.B.begin(),req.B.end());	// sort B1,...BN in acescending order		
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	//intiate the map variables:
	map.Clear();
	
	//quick fail 
	float res_port_B[nServer]={0},sumB=0;
	if (QuickFail(req,sumB,res_port_B)==false)
		return false;
	Matrix<bool>Q(Nv,req.N);
	for(int x=0;x<req.N;x++)
			for(int y=x;y<req.N;y++)
				split(x,y)=getSplit(req,sumB,x,y);
	//from down to root
	for (int v=0;v<Nv;v++)
		for(int x=0;x<req.N;x++)
			for(int y=x;y<req.N;y++){
			calculateAR(req,sumB,x,y,v);
			
			//cout<<" "<<v<<" "<<x<<" "<<y<<" ";
			//cout<<" "<<AR[v].split[x][y];
		}
	//from root to down
	bool accept=subAssign(req,0,req.N-1,root,Q,map);
	if(accept){
		for(int i=0;i<req.N;i++)for(int j=0;j<nServer;j++)
			if(Q(j,i)){
				assignment[i]=j;
				map.Slot[j]++;
			}
	}
	return accept;
			
}
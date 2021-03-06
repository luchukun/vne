#include <ilcplex/ilocplex.h>
#include<ctime>
#include "frandom.h"
#include "elements.h"
#include "Graphs.h"
ILOSTLBEGIN
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVarCubic;
using namespace std;
Graph::Graph()
{
	for( int i = 0; i < Nv; i++ )
	{	
		Table[i].Degree=0;
		Table[i].Dist = Infinity;
		Table[i].Prev = -1;//pointer to its self
		Table[i].Visited = 0;
		Table[i].Type=1;
		Table[i].Slot=-1;
	}

	for (int s=0;s<nServer;s++)
		for(int d=0;d<nServer;d++)
		{
			n_path[s][d]=0;
			for (int k=0;k<K;k++)
				for(int h=0;h<Nv;h++)
					Kpaths[s][d][k][h]=0;
			for (int i=0;i<Ne;i++)
			{
				f[s][d][i]=0;
				b[s][d][i]=false;
			}
		}
		Kmax=Kwidest;
		routingOption=default_routing;
}
Graph::~Graph()
{
	for( int i = 0; i < Nv; i++ )
		delete []Table[i].Adj;
}
void Graph::drawVL2()
{	
#ifdef _VL2
	float Cost=1;
	int Addr=0;
	//allocate *Adj according to node degree
	for (int i=0;i<nServer;i++)
		Table[i].Adj=new Edge[n_server_port]; //two adjcent edges for each server;
	for (int i=0;i<nToR;i++)
		Table[i+nServer].Adj=new Edge[nServerInRack+2];
	for (int i=0;i<nAggregate;i++)
		Table[i+nServer+nToR].Adj=new Edge[DA];
	for (int i=0;i<nIntermediate;i++)
		Table[i+nServer+nToR+nAggregate].Adj=new Edge[DI];

	//identify the servers
	for (int i=0;i<nServer;i++){
		Table[i].Type=0;	Table[i].Slot=maxSlot;
	}
	//layer 1
	for (int i=0;i<nToR;i++)
		for (int j=0;j<nServerInRack;j++)
		{	int Server_index=nServerInRack*i+j;
			int ToR_index=nServer+i;
			AddEdge(Server_index,ToR_index,Cost,gbpsServer2ToR,Addr++);
			Table[Server_index].Type=0;	Table[Server_index].Slot=maxSlot;
		}
	//layer 2
	for (int i=0;i<nToR;i++)
	{
		int ToR_index=nServer+i;
		int Groop_index=i/(DA/2);
		int Agg_index0=nServer+nToR+Groop_index*2;
		int Agg_index1=nServer+nToR+Groop_index*2+1;
		AddEdge(ToR_index,Agg_index0,Cost,gbpsToR2Agg,Addr++);
		AddEdge(ToR_index,Agg_index1,Cost,gbpsToR2Agg,Addr++);
		Table[ToR_index].Type=1;
	}
	//layer 3
	for (int i=0;i<nAggregate;i++)
		for (int j=0;j<nIntermediate;j++)
	{
		int Agg_index=nServer+nToR+i;
		int Int_index=nServer+nToR+nAggregate+j;
		AddEdge(Agg_index,Int_index,Cost,gbpsAgg2Int,Addr++);
		Table[Agg_index].Type=1;
		Table[Int_index].Type=1;
	}
#endif
}
void Graph::drawL2VL2()
{	
#ifdef _L2VL2
	float Cost=1;
	int Addr=0;
	//allocate *Adj according to node degree
	for (int i=0;i<nServer;i++){
		Table[i].Adj=new Edge[n_server_port]; //two adjcent edges for each server;
		Table[i].Type=0;	Table[i].Slot=maxSlot;	//identify the servers
	}
	for (int i=0;i<nToR;i++){
		Table[i+nServer].Adj=new Edge[nServerInRack+2];
		Table[i+nServer].Type=1;	// identify the switches
	}
	for (int i=0;i<nEdge;i++){
		Table[i+nServer+nToR].Adj=new Edge[nEdgePort];
		Table[i+nServer+nToR].Type=1;// identify the switches
	}

	//layer 1
	for (int i=0;i<nToR;i++){
		for (int j=0;j<nServerInRack;j++)
		{	int Server_index=nServerInRack*i+j;
			int ToR_index=nServer+i;
			AddEdge(Server_index,ToR_index,Cost,gbpsServer2ToR,Addr++);
		}
	}
	//layer 2
	for (int i=0;i<nEdge;i++)
	{	
		int Edge_index=nServer+nToR+i;
		for (int j=0;j<nToR;j++)
		{
			int ToR_index=nServer+j;		
			AddEdge(ToR_index,Edge_index,Cost,gbpsToR2Edge,Addr++);
		}
	}
#endif
}
void Graph::drawFattree()
{
#ifdef _FatTree
	float Cost=1;
	int Addr=0;

	//allocate *Adj according to node degree
	for (int i=0;i<nServer;i++)
		Table[i].Adj=new Edge[n_server_port]; //two adjcent edges for each server;
	for (int i=0;i<nToR+nAggregate+nCore;i++)
		Table[i+nServer].Adj=new Edge[H];
	//identify the servers
	for (int i=0;i<nServer;i++){
		Table[i].Type=0;	Table[i].Slot=maxSlot;
	}
	//layer 1
	for (int i=0;i<nToR;i++)
		for (int j=0;j<nServerInRack;j++)
		{	int Server_index=nServerInRack*i+j;
			int ToR_index=nServer+i;
			AddEdge(Server_index,ToR_index,Cost,gbpsCommodity*oversubstription,Addr++);
		}
	//layer 2
	for (int i=0;i<nToR;i++)
		for (int j=0;j<H/2;j++)
	{
		int ToR_index=nServer+i;//edge switch
		int Pod_index=i/(H/2);
		int Agg_index=nServer+nToR+Pod_index*H/2+j;
		AddEdge(ToR_index,Agg_index,Cost,gbpsCommodity*oversubstription,Addr++);
		Table[ToR_index].Type=1;
	}
	//layer 3
	for (int i=0;i<nCore;i++)
		for (int j=0;j<H;j++)// connect to H pod
	{   
		int Core_index=nServer+nToR+nAggregate+i;
		int Agg_index=nServer+nToR+j*(H/2)+i/(H/2);
		AddEdge(Core_index,Agg_index,Cost,gbpsCommodity*oversubstription,Addr++);
		Table[Core_index].Type=1;
		Table[Agg_index].Type=1;
	}
#endif

}
void Graph::drawBcube()
{
#ifdef _Bcube
	float Cost=1;
	int Addr=0;

	//allocate *Adj according to node degree
	for (int i=0;i<nServer;i++){
		Table[i].Adj=new Edge[n_server_port]; //two adjcent edges for each server;
		Table[i].Type=0;	Table[i].Slot=maxSlot;//identify the servers
	}
	for (int i=0;i<nToR+nAggregate;i++){
		Table[i+nServer].Adj=new Edge[H];
		Table[i+nServer].Type=1;Table[i+nServer].Slot=-1;
	}

	//layer 1
	for (int i=0;i<nToR;i++)
		for (int j=0;j<nServerInRack;j++)
		{	int Server_index=nServerInRack*i+j;
			int ToR_index=nServer+i;
			AddEdge(Server_index,ToR_index,Cost,gbpsCommodity,Addr++);
		}
	//layer 2
	for (int j=0;j<H;j++) // the jth group
		for (int i=0;i<nAggregate;i++)
		{
			int Server_index=nServerInRack*j+i;
			int Agg_index=nServer+nToR+i;
			AddEdge(Server_index,Agg_index,Cost,gbpsCommodity,Addr++);
		}	
#endif

}
int Graph::drawGragh(string &file_name)
{
	//const string file_name("Text1.txt");
	const int MaxLineLength = 256;
    static char OneLine[ MaxLineLength + 1 ];
    int Vs, Vd, Slot,Addr=0;
    float Bandwidth,Cost=1;
	if( file_name !="" )
    {             		
		ifstream Gs(file_name.c_str( ));
		if( Gs )
		{
			while( Gs.getline( OneLine, MaxLineLength ) )
				{
					istrstream Ls( OneLine, MaxLineLength );
					if (Ls>>Vs&&Ls>>Slot)
					{
						Table[Vs-1].Type=0;// the node is a server
						Table[Vs-1].Slot=Slot;
					}
					else if( Ls >> Vs && Ls >> Vd &&Ls>>Cost &&Ls>>Bandwidth)
					{	
						//std::cout<<Vs<<"-"<<Vd<<"-"<<Cost<<Bandwidth<<endl;
						AddEdge( Vs-1, Vd-1, Cost, Bandwidth,Addr++);
					}
					else
						 
						 std::cerr << "Bad line: "<<endl ;
				}
			Gs.close();
		}
        else
			std::cerr << "Error opening " << file_name<<endl ;
    } 
		
	return 1;
		
}

// Add the edge ( Source, Dest, Cost, Bandwidth ) to the Graph

void Graph::AddEdge( const int & s, const int & d, const float &cost,const float &bandw,const int&addr)
{	
	int i;
	//	e(s,d)
	i=Table[s].Degree;	
	(Table[s].Adj)[i].Assign(d,cost,addr);
	Bandwidth[addr]=bandw;	resBandwidth[addr]=bandw;
	Table[s].Degree++;	
	//	e(d,s)
	
	i=Table[d].Degree;
	(Table[d].Adj)[i].Assign(s,cost,addr+Ne/2);
	Bandwidth[addr+Ne/2]=bandw;	resBandwidth[addr+Ne/2]=bandw;
	Table[d].Degree++;

	opposite[addr]=addr+Ne/2;	opposite[addr+Ne/2]=addr;
}

void Graph::ClearTable()
{
	//for routing
	for( int i = 0; i < Nv; i++ )
	{	
		Table[i].Dist = Infinity;
		Table[i].Prev = -1;//pointer to its self
		Table[i].Visited = 0;
	}

}


void Graph::ClearNetwork()
{
	ClearTable();
	//renew server space
	for (int i=0;i<nServer;i++){
		Table[i].Slot=maxSlot;
	}
	// renew link capacity
	for (int e=0;e<Ne;e++){	
		resBandwidth[e]=Bandwidth[e];
	}		
}
void Graph::RandomNetwork(float p,float p_minResBw,float p_maxResBw)
{
	//renew server space
	for (int i=0;i<nServer;i++){
		if (rand_b01(p))
			Table[i].Slot=unif_int(0,maxSlot);
		else
			Table[i].Slot=maxSlot;
	}
	// renew link capacity
	for (int e=0;e<Ne/2;e++)
	{	
		if (rand_b01(p))
			resBandwidth[e]=unif_int(Bandwidth[e]*p_minResBw,Bandwidth[e]*p_maxResBw)/100*100;
		else
			resBandwidth[e]=Bandwidth[e];
		resBandwidth[e+Ne/2]=resBandwidth[e];
	}
}
void Graph::PrintPathRec( int DestNode ) const
{
    if( Table[DestNode].Prev !=-1 )
    {
        PrintPathRec( Table[DestNode].Prev );
		std::cout << " to ";
    }
    std::cout << DestNode+1;
}

// Driver routine to handle unreachables and print total cost
// It calls recursive routine to print shortest path to
// DestNode after a hortest path algorithm has run

void Graph::PrintPath(int s,int d) const
{
	if( Table[d].Dist ==Infinity)
		std::cout << d+1 << " is unreachable";
	else
	{	std::cout <<"The route to Node "<<d+1<<" is: ";
		PrintPathRec(d);
		std::cout << " (cost is " << Table[d].Dist << ")";
	}
	std::cout << endl;
		
}



int Graph::Dijkstra(int s)
{
	int V=-1,W;
    BinaryHeap<Comparable> PQ(Nv);
    ClearTable();
    Table[s].Dist = 0;
	
    PQ.insert(Comparable(s, 0));
	for( int NodesSeen = 0; NodesSeen < Nv; NodesSeen++  )
    {    
        do
        {	if( PQ.isEmpty( ))
			{	cout<<"PQ is Empty";
				return 0;
			}
            V = PQ.deleteMin().index;
        } while(Table[V].Visited);
		//std::cout<<" "<<V+1<<" ";
		Table[V].Visited = 1;     // Visited vertex as being seen
		
		for (int P =0;P<Table[V].Degree; P=P+1)
		{
			W = (Table[V].Adj)[P].Dest;
			if (Table[W].Visited==false) {
				float Cvw = (Table[V].Adj)[P].Cost;
				if( Cvw < 0 )
				{
					std::cerr << "Graph has negative edges" << endl;
					return 0;
				}
				if( Table[W].Dist > Table[V].Dist + Cvw )
				{
					Table[W].Dist = Table[V].Dist + Cvw;
					Table[W].Prev = V;
					PQ.insert(Comparable( W, Table[W].Dist ) );
				}
			}
		}
    }
	
	return 1;
}
float Graph::Dijkstra(int s, int d)
{
	int V=-1,W;
    BinaryHeap<Comparable> PQ(Nv);
    Comparable next;    // Stores the result of a DeleteMin
	int P;
    ClearTable();
    Table[s].Dist = 0;
	
    PQ.insert(Comparable(s, 0));
	while(Table[d].Visited ==false)
    {    
        do
        {	if( PQ.isEmpty( ))
			{	//cout<<"PQ is Empty";
				return 0;
			}
            V = PQ.deleteMin().index;
        } while(Table[V].Visited);

		//std::cout<<" "<<V+1<<" ";
		Table[V].Visited = 1;     // Visited vertex as being seen
		for (P =0;P<Table[V].Degree; P=P+1)
		{
			W = (Table[V].Adj)[P].Dest;
			if (Table[W].Visited==false) {
				float Cvw = (Table[V].Adj)[P].Cost;
				if( Cvw < 0 )
				{
					std::cerr << "Graph has negative edges" << endl;
					return 0;
				}
				if( Table[W].Dist > Table[V].Dist + Cvw )
				{
					Table[W].Dist = Table[V].Dist + Cvw;
					Table[W].Prev = V;
					PQ.insert(Comparable( W, Table[W].Dist ) );
				}
			}
		}
    }
	return Table[d].Dist;
}
// compute K shortest paths from s to d for each s=1:n,d=1:n
void Graph::YenKSP(int s, int d)
{	
	Dijkstra(s,d);//A[0] =Dijkstra(s,d);
	if (Table[s].Dist==Infinity){
		n_path[s][d]=0;
		return;
	}
	n_path[s][d]=1;
    
	Kpaths[s][d][0][0]=(unsigned short)Table[d].Dist;
	unsigned short relay=d;
	for (unsigned short h=Kpaths[s][d][0][0]+1;;h--)
    {
	   Kpaths[s][d][0][h]=relay;
	   if (relay==s)
		   break;
	   relay=Table[relay].Prev;					
    }

   // Determine the shortest path from the source to the sink.
	unsigned short (*A)[Nv+1]=Kpaths[s][d];	
	unsigned short rootPath[Nv+1]={0};   

   // Initialize the heap to store the potential kth shortest path.
	vector<vector<unsigned short>>B(K*Nv,vector<unsigned short>(Nv+1));
	BinaryHeap<Comparable> PQ(K*Nv);
	for(int k=1;k<=K-1;k++) 
	{
	   // The spur node ranges from the first node to the next to last node in the previous k-shortest path.
	   for (unsigned short i=1;i< A[k-1][0] ;i++)
	   {
		   // Spur node is retrieved from the previous k-shortest path, k − 1.
		   unsigned short spurNode = A[k-1][i];
		   // The sequence of nodes from the source to the spur node of the previous k-shortest path.
		   vector<unsigned short> rootPath(A[k-1],A[k-1]+i+1);	//rootPath = A[k-1].nodes(0, i);
		   rootPath[0]=i-1;// length of rootpath
		   for (int p=0;p<k;p++)	// for each path p in A:
		   {
			   // Remove the links that are part of the previous shortest paths which share the same root path.
			   bool flag=true;	
			   for(int j=1;j<=i;j++)
				   if (rootPath[j]!=A[p][j]){	flag=false;break;	} 
			   
			   if (flag==true)	//if rootPath == p.nodes(0, i):
				{	unsigned short next2spur=A[p][i+1];	//remove p.edge(i, i + 1) from Graph;
					for (int j =0;j<Table[spurNode].Degree; j=j+1)
					{
						if((Table[spurNode].Adj)[j].Dest==next2spur)
						{
							(Table[spurNode].Adj)[j].Cost=Infinity;break;
						}
					}
			   }
		   }
           for(int x=1;x<i;x++)//for each node in rootPath except spurNode:
		   {
			   unsigned short removeNode=rootPath[x];	//remove it from Graph;
			   for (int j =0;j<Table[removeNode].Degree; j=j+1)
				   (Table[removeNode].Adj)[j].Cost=Infinity;
		   }
		           			  
		   
		   // Calculate the spur path from the spur node to the sink.
		   Dijkstra(spurNode,d);
		   // Add back the edges and nodes that were removed from the graph.
		   for (int p=0;p<k;p++)	// for each path p in A:
		   {
			   	unsigned short next2spur=A[p][i+1];	//add p.edge(i, i + 1) to Graph;
				for (int j =0;j<Table[spurNode].Degree; j=j+1)
				{
					if((Table[spurNode].Adj)[j].Dest==next2spur)
					{
						(Table[spurNode].Adj)[j].Cost=1;break; //restore edges to Graph;
					}
				}
			   
		   }
		   //restore nodes in rootPath to Graph;
		   for(int x=1;x<i;x++)//for each node in rootPath except spurNode:
		   {
			   unsigned short removeNode=rootPath[x];	//add to Graph;
			   for (int j =0;j<Table[removeNode].Degree; j=j+1)
				   (Table[removeNode].Adj)[j].Cost=1;
		   }

		   if (Table[d].Dist==Infinity)
			   continue;
		   // totalPath = rootPath + spurPath;
		   unsigned short len_fullpath=rootPath[0]+(unsigned short)Table[d].Dist;
           vector<unsigned short>spurPath(len_fullpath+2);
		   int relay=d;
		   for (int j=len_fullpath+1;;j--)
		   {	
			   spurPath[j]=relay;
			   relay=Table[relay].Prev;			   
			   if (relay==spurNode)
				   break;						
		   } 
		   for(int j=1;j<=rootPath[0]+1;j++)
			   spurPath[j]=rootPath[j];
		   spurPath[0]=len_fullpath;

		   B.push_back(spurPath);
		   // Add the potential k-shortest path to the heap.
		   PQ.insert(Comparable(B.size()-1,(int)len_fullpath));
           
		   
	   }
	   // This handles the case of there being no spur paths, or no spur paths left.
	   if (PQ.isEmpty()||B.empty())
		   break;
	   // Sort the potential k-shortest paths by cost.
	   
	   int idx=PQ.deleteMin().index;
	   // Add the lowest cost path becomes the k-shortest path.
	   for (int j=0;j<= B[idx][0]+1;j++)
		   A[k][j] = B[idx][j];
	   // it's equal shortest path
	   if ((routingOption==_ECMP&&A[k][0]==A[0][0])||routingOption==_KshortestLB) 
		   n_path[s][d]++;

	 }	//for each k


	for(int k=0;k<min((int)n_path[s][d],K);k++)  //每条path
	{	
		
		for(int h=1;h<=A[k][0];h++)  //path上每条link
		{
			unsigned short p=A[k][h],q=A[k][h+1]; //link e(p,q)
			for (int x =0;x<Table[p].Degree; x=x+1)
			{
				if((Table[p].Adj)[x].Dest==(int)q)
				{
					int e=(Table[p].Adj)[x].Addr;
					b[s][d][e]= true; //
				}
			}
		}		
	}

}

float Graph::LoadBalance(int s,int d)
{	
	unsigned short (*A)[Nv+1]=Kpaths[s][d];
	float capacity[K]={0};// the capacity of k-th shortest path
	float sum_capacity = 0;
		
	//--- 计算path capacity ---	
	for(int k=0;k<min((int)n_path[s][d],K);k++) 
	{
		capacity[k] = Infinity;
		if (A[k][0]==0) {
			capacity[k] =0; 
			continue;
		}
		for(int h=1;h<=A[k][0];h++)
		{
			unsigned short p=A[k][h],q=A[k][h+1]; //link e(p,q)
			
			for (int x =0;x<Table[p].Degree; x=x+1)
			{
				if((Table[p].Adj)[x].Dest==(int)q)
				{	int e=(Table[p].Adj)[x].Addr;
					capacity[k]=min(capacity[k],resBandwidth[e]);
				}
			}
		}
		//sum_capacity += capacity[k];
	}
	
	//--- route with load balance ---
	for (int e=0;e<Ne;e++) 
		f[s][d][e]=0;
	
	bool path_selected[K]={0};
	
	for(int k=0;k<min((int)n_path[s][d],Kmax);k++)
	{	
		float max_cap=-1;
		int the_path=-1;
		for(int i=0;i<min((int)n_path[s][d],K);i++)
		{ 
			if((max_cap<capacity[i])&&(!path_selected[i]))
			{
				max_cap=capacity[i];
				the_path=i;
			}
		}
		path_selected[the_path]=true;
		sum_capacity += capacity[the_path];
	}

	if (sum_capacity==0)
		return sum_capacity;

	float split_ratio[K]={0};
	//++++++++++可交叉多路径主要体现在这里+++++++++	
	for(int k=0;k<min((int)n_path[s][d],K);k++)  //每条path
	{	
		if (path_selected[k]){
			for(int h=1;h<=A[k][0];h++)  //path上每条link
			{
				unsigned short p=A[k][h],q=A[k][h+1]; //link e(p,q)
				for (int x =0;x<Table[p].Degree; x=x+1)
				{
					if((Table[p].Adj)[x].Dest==(int)q)
					{
						int e=(Table[p].Adj)[x].Addr;
						split_ratio[k]=capacity[k]/sum_capacity;
						f[s][d][e]+= split_ratio[k] ; //f-sd在e(p,q)上的分配比例
						 // which path is routed through e
					}
				}
			}
		}
		//cout<<"  f on path "<<k+1<<": "<<capacity[k]/sum_capacity<<endl;
	}
	/*------------ the real capacity----------------------- */
	float real_capacity[K]={0},real_sum_capacity=0;

	for(int k=0;k<min((int)n_path[s][d],K);k++)  //每条path
	{	real_capacity[k] = Infinity;
		if (path_selected[k]){
			for(int h=1;h<=A[k][0];h++)  //path上每条link
			{
				unsigned short p=A[k][h],q=A[k][h+1]; //link e(p,q)
				for (int x =0;x<Table[p].Degree; x=x+1)
				{
					if((Table[p].Adj)[x].Dest==(int)q)
					{
						int e=(Table[p].Adj)[x].Addr;
						real_capacity[k]=min(real_capacity[k],resBandwidth[e]*split_ratio[k]/f[s][d][e]);
						
					}
				}
			}
			real_sum_capacity+=real_capacity[k];
		}
		//cout<<"  f on path "<<k+1<<": "<<capacity[k]/sum_capacity<<endl;
	}

	return real_sum_capacity;
}


float Graph::ECMP(int s,int d)
{	
	unsigned short (*A)[Nv+1]=Kpaths[s][d];
	float capacity[K]={0};// the capacity of k-th shortest path
	float sum_capacity = 0;
		
	//--- 计算path capacity ---	
	for(int k=0;k<min((int)n_path[s][d],K);k++) 
	{
		capacity[k] = Infinity;
		if (A[k][0]==0) {
			capacity[k] =0; 
			continue;
		}
		for(int h=1;h<=A[k][0];h++)
		{
			unsigned short p=A[k][h],q=A[k][h+1]; //link e(p,q)
			
			for (int x =0;x<Table[p].Degree; x=x+1)
			{
				if((Table[p].Adj)[x].Dest==(int)q)
				{	int e=(Table[p].Adj)[x].Addr;
					capacity[k]=min(capacity[k],resBandwidth[e]);
				}
			}
		}
		//sum_capacity += capacity[k];
	}
	
	//--- route with load balance ---
	for (int e=0;e<Ne;e++) 
		f[s][d][e]=0;
	
	bool path_selected[K]={0};
	
	for(int k=0;k<min((int)n_path[s][d],Kmax);k++)
	{	
		float max_cap=-1;
		int the_path=-1;
		for(int i=0;i<min((int)n_path[s][d],K);i++)
		{ 
			if((max_cap<capacity[i])&&(!path_selected[i]))
			{
				max_cap=capacity[i];
				the_path=i;
			}
		}
		path_selected[the_path]=true;
		sum_capacity += capacity[the_path];
	}

	if (sum_capacity==0)
		return sum_capacity;

	float split_ratio[K]={0};
	//++++++++++可交叉多路径主要体现在这里+++++++++	
	for(int k=0;k<min((int)n_path[s][d],K);k++)  //每条path
	{	
		if (path_selected[k]){
			split_ratio[k]=1/(float)n_path[s][d];//capacity[k]/sum_capacity;

			for(int h=1;h<=A[k][0];h++)  //path上每条link
			{
				unsigned short p=A[k][h],q=A[k][h+1]; //link e(p,q)
				for (int x =0;x<Table[p].Degree; x=x+1)
				{
					if((Table[p].Adj)[x].Dest==(int)q)
					{
						int e=(Table[p].Adj)[x].Addr;
						
						f[s][d][e]+= split_ratio[k] ; //f-sd在e(p,q)上的分配比例
						 // which path is routed through e
					}
				}
			}
		}
		//cout<<"  f on path "<<k+1<<": "<<capacity[k]/sum_capacity<<endl;
	}
	/*------------ the real capacity----------------------- */
	float real_capacity[K]={0},real_sum_capacity=0;

	for(int k=0;k<min((int)n_path[s][d],K);k++)  //每条path
	{	real_capacity[k] = Infinity;
		if (path_selected[k]){
			for(int h=1;h<=A[k][0];h++)  //path上每条link
			{
				unsigned short p=A[k][h],q=A[k][h+1]; //link e(p,q)
				for (int x =0;x<Table[p].Degree; x=x+1)
				{
					if((Table[p].Adj)[x].Dest==(int)q)
					{
						int e=(Table[p].Adj)[x].Addr;
						real_capacity[k]=min(real_capacity[k],resBandwidth[e]*split_ratio[k]/f[s][d][e]);
						
					}
				}
			}
			real_sum_capacity+=real_capacity[k];
		}
		//cout<<"  f on path "<<k+1<<": "<<capacity[k]/sum_capacity<<endl;
	}

	return real_sum_capacity;
}


float Graph::MaxFlow(int s,int d)
{	
	if(s==d)
		   return 0;
	//--------基本环境设置--------
	IloEnv env; //环境environment
	IloModel model(env); //model
		
	//-------常量变量设置--------

	//-------待定变量设置------ f[s][d][e]
	IloNumVarArray fx(env, Ne, 0, IloInfinity );	

//----------cplex列方程---------	
	//----优化目标----
	IloNumVar fmax (env);
	model.add(IloMaximize(env,fmax));	
	/* force to routing on equal cost paths*/
	for (int e=0;e<Ne;e++)
		if (b[s][d][e]==0)
			model.add(fx[e]==0);

	//--flow constraints--//
	for (int e=0;e<Ne;e++)
		model.add(fx[e]<=resBandwidth[e]);
	//flow conservation		

	   //C4+: sum{f[s][d][e(s,+)]}=f
	   IloExpr flow_s_pos(env);
	   
	   for (int p =0;p<Table[s].Degree; p=p+1)
	   {
		   int e=(Table[s].Adj)[p].Addr;
		   flow_s_pos+=fx[e];
	   }
	   //model.add( flow_s_pos==fmax);
	   //C5+: sum{f[s][d][e(d,-)]}=f
	   IloExpr flow_d_neg(env);
	   for(int v=0;v<Nv;v++)for(int p =0;p<Table[v].Degree;p=p+1)
	   {	
		   if ((Table[v].Adj)[p].Dest==d)
		   {
				int e=(Table[v].Adj)[p].Addr;
				flow_d_neg+=fx[e];
		   }
	   }
	   //model.add( flow_d_neg==fmax);

	   //C4-: sum{f[s][d][e(d,+)]}=0
	   IloExpr flow_d_pos(env);
	   
	   for (int p =0;p<Table[d].Degree; p=p+1)
	   {
		   int e=(Table[d].Adj)[p].Addr;
		   flow_d_pos+=fx[e];
	   }
	   //model.add( flow_d_pos==0);
	   
	   //C5-: sum{f[s][d][e(s,-)]}=0
	   IloExpr flow_s_neg(env);
	   for(int v=0;v<Nv;v++)for(int p =0;p<Table[v].Degree;p=p+1)
	   {	
		   if ((Table[v].Adj)[p].Dest==s)
		   {
				int e=(Table[v].Adj)[p].Addr;
				flow_s_neg+=fx[e];
		   }
	   }
	  // model.add( flow_s_neg==0);
	   model.add(flow_s_pos-flow_s_neg==fmax);
	   model.add(flow_d_pos-flow_d_neg+fmax==0);
	   //C6: sum{f[s][d][e(u,+)]}=sum{f[s][d][e(u,-)]}
	   for (int v=0;v<Nv;v++)
	   {
		   if (v==s||v==d)
			   continue;
		   IloExpr flow_v_pos(env),flow_v_neg(env);
		   for (int p =0;p<Table[v].Degree; p=p+1)
		   {	
			   int e1=(Table[v].Adj)[p].Addr;
			   flow_v_pos+=fx[e1];
		   }
		   for (int w=0;w<Nv;w++)for (int q =0;q<Table[w].Degree; q=q+1)
		   {	
			   if ((Table[w].Adj)[q].Dest==v)
			   {
					int e2=(Table[w].Adj)[q].Addr;
					flow_v_neg+=fx[e2];
			   }
		   }
		   model.add(flow_v_pos==flow_v_neg);
	   }				   

	//---------Cplex解方程---------	   
	IloCplex cplex(model); //依model建的cplex
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::RootAlg, IloCplex::Network);
	if(cplex.solve())
	{
		float tar =(float) cplex.getObjValue();
		//cplex.exportModel("max_flow.lp");		
		for(int e=0;e<Ne;e++) 
			f[s][d][e]=cplex.getValue(fx[e])/tar;		
		env.end();
		return tar;

	}
	else
	{
		cout<<"No Solution!\n\n";
		cplex.exportModel("max_flow.lp");
		env.end();
		return -1;
	}
}
int Graph::Sort_fsd(int e)
{	
	int num=0;
	for (int s=0;s<nServer;s++)
		for(int d=0;d<nServer;d++)
			if(f[s][d][e]>0)
			{
				pass[e][num]=Pair(s,d);
				num++;
			}

	//--------对pass[e]内的pair排序
	
	for(int a=0;a<num-1;a++)
		for(int b=a;b<num;b++)
		{
			int sa=pass[e][a].s,da=pass[e][a].d;
			int sb=pass[e][b].s,db=pass[e][b].d;

			if(f[sa][da][e] < f[sb][db][e]) //f-sd-e(i,j)
			{
				//交换a,b位置
				Pair temp(pass[e][a]); 
				pass[e][a] = pass[e][b];
				pass[e][b] = temp;

			}
		}
	return num; 
}

void Graph::updateBottlenecks(float *hostBw,const Solution&map,float *TS,float *TD)
{
	float mlu=0;	int the_link=0;
	for (int e=0;e<Ne;e++){
		if (mlu<map.Bandwidth[e]/resBandwidth[e])
		{
			mlu=map.Bandwidth[e]/resBandwidth[e];
			the_link=e; //find the most congested link
		}
	}

		//new method to find bottleneck 
	//float TS[nServer]={0},TD[nServer]={0};
	int s,d;
	for(int p=0;p<n_pair[the_link];p++)
	{
		s = pass[the_link][p].s; d = pass[the_link][p].d;
		if (hostBw[s]>0&&hostBw[d]>0){
#ifdef _BottleckXHostBw
			TS[s]+=f[s][d][the_link]*min(hostBw[s],hostBw[d]);
			TD[d]+=f[s][d][the_link]*min(hostBw[s],hostBw[d]);
#else
			TS[s]+=f[s][d][the_link];
			TD[d]+=f[s][d][the_link];

#endif		
		
		}
	}


}

int Graph::findBottleneck(float *hostBw,const Solution&map)
{
	float mlu=0;	int the_link=0;
	for (int e=0;e<Ne;e++){
		if (mlu<map.Bandwidth[e]/resBandwidth[e])
		{
			mlu=map.Bandwidth[e]/resBandwidth[e];
			the_link=e; //find the most congested link
		}
	}
	int s,d,the_server=-1;
	Pair bottleneck(-1,-1);
	//寻找bottleneck server pair
	for (int p=0;p<n_pair[the_link];p++)
	{
		s = pass[the_link][p].s; d = pass[the_link][p].d;
		if (hostBw[s]>0&&hostBw[d]>0)
			bottleneck=pass[the_link][p];
	}
	//new method to find bottleneck 
	float TS[nServer]={0},TD[nServer]={0};
	for(int p=0;p<n_pair[the_link];p++)
	{
		s = pass[the_link][p].s; d = pass[the_link][p].d;
		if (hostBw[s]>0&&hostBw[d]>0){
		//TS[s]+=f[s][d][the_link]*min(hostBw[s],hostBw[d]);
		//TD[d]+=f[s][d][the_link]*min(hostBw[s],hostBw[d]);
		TS[s]+=f[s][d][the_link];
		TD[d]+=f[s][d][the_link];
		}
	}
#ifdef _BottleckXHostBw
	for(int x=0;x<nServer;x++){
		TS[x]=TS[x]*hostBw[x];
		TD[x]=TD[x]*hostBw[x];
	}
#endif

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

	if (rand_b01(0.5))
		the_server=bottleneck.s;
	else	
		the_server=bottleneck.d;
	return the_server;

}

int Graph::findBottleneck(float *hostBw,float *TS,float *TD)
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

int Graph::findBottleneck(int the_link,float *hostBw)
{
	Pair bottleneck(-1,-1);
	int s,d,the_server=-1;
	//寻找bottleneck server pair
	for (int p=0;p<n_pair[the_link];p++)
	{
		s = pass[the_link][p].s; d = pass[the_link][p].d;
		if (hostBw[s]>0&&hostBw[d]>0)
			bottleneck=pass[the_link][p];
	}
	//new method to find bottleneck 
	float TS[nServer]={0},TD[nServer]={0};
	for(int p=0;p<n_pair[the_link];p++)
	{
		s = pass[the_link][p].s; d = pass[the_link][p].d;
		if (hostBw[s]>0&&hostBw[d]>0){
		//TS[s]+=f[s][d][the_link]*min(hostBw[s],hostBw[d]);
		//TD[d]+=f[s][d][the_link]*min(hostBw[s],hostBw[d]);
		TS[s]+=f[s][d][the_link];
		TD[d]+=f[s][d][the_link];
		}
	}
	
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
	if (rand_b01(0.5))
		the_server=bottleneck.s;
	else	
		the_server=bottleneck.d;
	return the_server;

}

float Graph::maxTraffic(int e,float* hostBw)
{ 
	float TS[nServer]={0},TD[nServer]={0};
	float TS_max[nServer]={0},TD_max[nServer]={0};
	float T_left=0,T_add=0;
	float sum_t=0; // sum of traffic on link e(i,j)
	if(n_pair[e]==0)
		return 0;

	//---- Initialize max concurrent traffic T_max[s] for each server s
	for(int p=0;p<n_pair[e];p++)
	{
		int s = pass[e][p].s; int d = pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			//从s出发的总流量限制
			TS_max[s] = hostBw[s];
			//到达d的总流量限制
			TD_max[d] = hostBw[d];
		}
	}
	/*for(int p=0;p<n_pair[e];p++)
	{
		int s = pass[e][p].s; int d = pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			//从s出发的总流量限制
			TS_max[s] += hostBw[d];
			//到达d的总流量限制
			TD_max[d] += hostBw[s];
		}
	}
	for(int p=0;p<n_pair[e];p++)
	{
		int s = pass[e][p].s; int d = pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			TS_max[s]=min(TS_max[s],hostBw[s]);//TS_max = min(hostBw[s],sum(hostBw[d]))
			TD_max[d]=min(TD_max[d],hostBw[d]);//TD_max = min(hostBw[d],sum(hostBw[s]))
		}
	}	
	*/
	//------done-对pass[e]内的pair排序
	
	//--------计算traffic,并寻找bn_s
	for(int p=0;p<n_pair[e];p++)
	{
		int s = pass[e][p].s; int d = pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			T_left = min(TS_max[s]-TS[s], TD_max[d]-TD[d]);
			T_add =max((float)0, min( min(hostBw[s],hostBw[d]), T_left));
			sum_t += T_add*f[s][d][e]; //max traffic
			TS[s] += T_add;	
			TD[d] += T_add;
		}
	}
	return sum_t;
}



float Graph::LPmaxTraffic(int e,float* hostBw)
{	
	float TS_max[nServer]={0},TD_max[nServer]={0};
	if(n_pair[e]==0)
		return 0;
	//---- Initialize max concurrent traffic T_max[s] for each server s 
	int numMe=0,num_TS_max=0,num_TD_max=0;
	for(int p=0;p<n_pair[e];p++)
	{
		int s = pass[e][p].s; int d = pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			//从s出发的总流量限制
			TS_max[s] += hostBw[d];
			//到达d的总流量限制
			TD_max[d] += hostBw[s];
			numMe++;
			
		}
	}
	for (int i=0;i<nServer;i++){
		TS_max[i] = min(TS_max[i],hostBw[i]);
		TD_max[i] = min(TD_max[i],hostBw[i]);
		if (TS_max[i]>0)	num_TS_max++;
		if (TD_max[i]>0)	num_TD_max++;
	}
	/* ----------greedy method----------------- */
	float TS[nServer]={0},TD[nServer]={0};
	float T_left=0,T_add=0;
	float sum_t=0; // sum of traffic on link e(i,j)

	for(int p=0;p<n_pair[e];p++)
	{
		int s = pass[e][p].s; int d = pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			T_left = min(TS_max[s]-TS[s], TD_max[d]-TD[d]);
			T_add =max((float)0, min( min(hostBw[s],hostBw[d]), T_left));
			sum_t += T_add*f[s][d][e]; //max traffic
			TS[s] += T_add;	
			TD[d] += T_add;
		}
	}

	if (num_TS_max==1||num_TD_max==1||sum_t>resBandwidth[e])
		return sum_t;
	/*------------ end of greedy------------ */

	vector<Pair> Me(numMe);
	
	for(int p=0,i=0;p<n_pair[e];p++)
	{
		int s = pass[e][p].s; int d = pass[e][p].d;
		if(hostBw[s]>0&&hostBw[d]>0)
		{
			Me[i++]=pass[e][p];
		}
	}
	

	//--------基本环境设置--------
	float tar=0;
	IloEnv env; //环境environment
	try{
		IloModel model(env); //model
					
		//-------常量变量设置--------
		IloNumVarArray x(env, numMe, 0, IloInfinity );
		for(int p=0;p<numMe;p++)
		{
			int s = Me[p].s; int d = Me[p].d;
			model.add(x[p]<=min(hostBw[s],hostBw[d]));
		}
		
		for (int s=0;s<nServer;s++)
		{
			if(TS_max[s]>0)
			{
				IloExpr Ts(env);
				for(int p=0;p<numMe;p++)
				{
					int i = Me[p].s; int j = Me[p].d;
					if(i==s)
						Ts+=x[p];
				}
				model.add(Ts<=TS_max[s]);
			}
				
			if(TD_max[s]>0)
			{
				IloExpr Td(env);
				for(int p=0;p<numMe;p++)
				{
					int i = Me[p].s; int j = Me[p].d;
					if(j==s)
						Td+=x[p];
				}
				model.add(Td<=TD_max[s]);
			}
		}
		IloExpr load_expr(env);
		for(int p=0;p<numMe;p++)
		{
			int s = Me[p].s; int d = Me[p].d;
			load_expr+=f[s][d][e]*x[p];
		}

		//model.add( load_expr>=0);
		//----优化目标----	
		model.add(IloMaximize(env, load_expr));
		//---------Cplex解方程---------	   
		IloCplex cplex(model); //依model建的cplex
		
		//cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
		cplex.setParam(IloCplex::RootAlg, IloCplex::Network);
		cplex.setOut(env.getNullStream());
		//cplex.setParam(IloCplex::RootAlg, IloCplex::Concurrent);
		if(cplex.solve())
		{
			tar =(float) cplex.getObjValue();	
			//if(tar>resBandwidth[e])
			//	cplex.exportModel("maxload.lp");
			/*if (abs(tar-greedy_load)>0.01)
					cout<<"\n"<<"the lp and max_traffic is not consistent!"<<"\n";		
			cout<<"\n"<<"load(greedy)="<<greedy_load<<";load(lp)="<<tar<<"\n";	*/
		}
		else
		{
			cout<<"No Solution!\n\n";
		}

	}

	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
	  cerr << "Error" << endl;
	}
	env.end();
	return tar;
}
void Graph::reserveResource(const Solution&map)
{	
	// bandwidth reservation
    for (int e=0;e<Ne;e++)	
		resBandwidth[e]-=map.Bandwidth[e];
	// VM slot reservation
	for (int i=0;i<nServer;i++)
		Table[i].Slot-=map.Slot[i];
}

void Graph::releaseResource(const Solution&map)
{
	// bandwidth release
    for (int e=0;e<Ne;e++)	
		resBandwidth[e]+=map.Bandwidth[e];

	// VM slot release
	for (int i=0;i<nServer;i++)
		Table[i].Slot+=map.Slot[i];
}

float Graph::ProcessRequest(char algorithm,bool enLProuting,int numOfreq,float load,
						float &max_utilization,float &success_rate,float &bandwidth_cost,float &RC)
{
	ClearNetwork();
	init_genrand(0);
	//init_genrand( (unsigned long)time( NULL ));				
	float t=0;
	for (int s=0;s<nServer;s++)	for (int d=0;d<nServer;d++)
	{
		if (s==d)
			continue;
		else
			YenKSP(s,d);
	}
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
			accept=Pertubation(req,enLProuting,map);
			break;
#ifdef _local_search
		case 'L':
			accept=LocalSearch(req,enLProuting,map);
			break;
#endif
		//case 'E':
		//	accept=ExhaustiveSearch(req,enLProuting,map);
		//	break;
		case 'F':
			accept=FirstFit(req,enLProuting,map);
			break;
		case 'G':
			accept=BestFit(req,enLProuting,map);
			break;
		case 'N':
			accept=NextFit(req,enLProuting,map);
			break;
		case 'B':
			accept=BackTracking(req,enLProuting,map);
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
			max_utilization+=max_ut;
			bandwidth_cost+=bw_used;
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

float Graph::SingleRequest(char algorithm,bool enLProuting,int numOfreq,
		float prob,float p_minResBw,float p_maxResBw,int  numOfVM,float minBw,float maxBw,
		float &max_utilization,float &success_rate,float &bandwidth_cost)
{
	ClearNetwork();
	init_genrand(0);
	//init_genrand( (unsigned long)time( NULL ));				
	float t=0;
	for (int s=0;s<nServer;s++)	for (int d=0;d<nServer;d++)
	{
		if (s==d)
			continue;
		else
			YenKSP(s,d);
	}
	Cluster req;	// embedding input
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
			accept=Pertubation(req,enLProuting,map);
			break;
		case 'R':
			accept=randomDrop(req,enLProuting,map);
			break;
		case 'B':
			accept=BackTracking(req,enLProuting,map);
			break;
		case 'F':
			accept=FirstFit(req,enLProuting,map);
			break;
		case 'G':
			accept=BestFit(req,enLProuting,map);
			break;
		case 'N':
			accept=NextFit(req,enLProuting,map);
			break;
		//case 'E':
		//	accept=ExhaustiveSearch(req,enLProuting,map);
		//	break;
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
				if (max_ut-1>1e-3)
					cout<<"/n"<<"max_ut="<<max_ut<<"/n";
			}			
			
			for (int e=0;e<Ne;e++){
				bw_used+=map.Bandwidth[e];
			}
			max_utilization+=max_ut;
			bandwidth_cost+=bw_used;//bandwidth_cost+=bw_used/bw_total;
		}	
	}
	success_rate=(float)n_success/(float)n_req;
	if(n_success>0){
		max_utilization=max_utilization/(float)n_success;
		bandwidth_cost=bandwidth_cost/(float)n_success;
	}
    return success_rate;
}


float Graph::LPRouting(const Cluster& req,float* hostBw,Solution&map)
{
	
	//-------常量变量设置--------
	float sumB=0;
	for (int i=0;i<req.N;i++)
		sumB+=req.B[i];
	
	//--------基本环境设置--------
	IloEnv env; //环境environment
	IloModel model(env); //model
		

	IloInt n_location=0; // number of used server
	for (int i=0;i<nServer;i++)
	   if (hostBw[i]>0)n_location++;
	if (n_location<=1){
		env.end();
		return 0;
	}
	IloNumArray m (env,n_location);
	IloIntArray the_server (env,n_location); //VM allocation
	for (int i=0,j=0;i<nServer;i++)
		if (hostBw[i]>0){
		   m[j]=hostBw[i]; the_server[j]=i;	j++;
		}
#ifdef _debug_LPRouting 
	float lb_mlu=0;
	// the mlu under load-balanced routing
	bool caculated[Ne]={0};
	for (int e=0;e<Ne;e++)
	{	
		if (caculated[e])
			continue;
		map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
		map.Bandwidth[opposite[e]]=map.Bandwidth[e];
		caculated[e]=true; caculated[opposite[e]]=true;

		lb_mlu=max(lb_mlu,1-(resBandwidth[e]-map.Bandwidth[e])/Bandwidth[e]);
	}
	/* display placement and routing information*/
	cout<<"\n"<<"load balance routing"<<"\n";
	for (int i=0;i<nServer;i++)
		cout<<hostBw[i]<<" ";
	for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
	{	if (s==d) continue;
		cout<<"\n"<<"f["<<the_server[s]<<"]["<<the_server[d]<<"]= ";
		for(int e=0;e<Ne;e++)
			if (f[the_server[s]][the_server[d]][e]>0)
			cout<<e<<" <"<<f[the_server[s]][the_server[d]][e]<<"> ";
	}

	if (lb_mlu<=0)
		return lb_mlu;
#endif
	//-------待定变量设置------ f[s][d][e]
	NumVarCubic fx(env, n_location); //
	for(int i=0;i<n_location;i++)
	{
	   fx[i] = NumVarMatrix(env,n_location);
	   for(int j=0;j<n_location;j++)
		   fx[i][j] = IloNumVarArray(env, Ne, 0, 1 );
	}
	// prime-dual variables
	NumVarCubic x(env, n_location); //x[s][d][e]
	for(int s=0;s<n_location;s++)
	{
	   x[s] = NumVarMatrix(env,n_location);  
	   for(int d=0;d<n_location;d++)
		   x[s][d] = IloNumVarArray(env, Ne, 0, IloInfinity );
	}
	NumVarMatrix y(env,n_location),z(env,n_location); //y[s][e],z[d][e]
	for(int s=0;s<n_location;s++)
	{
	   y[s] = IloNumVarArray(env, Ne, 0, IloInfinity );
	   z[s] = IloNumVarArray(env, Ne, 0, IloInfinity );
	}

//----------cplex列方程---------	
	//----优化目标----
	IloNumVar mlu (env);
	model.add(IloMinimize(env,mlu ));	
	// ---constraints-----//
	model.add( mlu <= 1);

	//--flow constraints--//

	//C1:	f[s][s][e] =0
	for(int e=0;e<Ne;e++)  
		for(int s=0;s<n_location;s++){
		   model.add( fx[s][s][e] == 0);
		}
	/* force to routing on equal cost paths*/
	for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)for(int e=0;e<Ne;e++)  
	{
		int Vs=the_server[s];int Vd=the_server[d];
		if(b[Vs][Vd][e]==0)
			model.add(fx[s][d][e]==0);
	}

	//C2:	f[s][d][e(v,w)]=f[d][s][e(w,v)]
	for(int s=0;s<n_location;s++)for(int d=s+1;d<n_location;d++)for(int e=0;e<Ne;e++) 
	{	
		model.add(fx[s][d][e]=fx[d][s][opposite[e]]); //e1(v,w) and e2(w,v)
	}
			
	//flow conservation		
	for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
	{
	   if(s==d)
		   continue;

	   int Vs=the_server[s];int Vd=the_server[d];

	   //C4+: sum{f[s][d][e(s,+)]}=1
	   IloExpr flow_s_pos(env);
	   
	   for (int p =0;p<Table[Vs].Degree; p=p+1)
	   {
		   int e=(Table[Vs].Adj)[p].Addr;
		   flow_s_pos+=fx[s][d][e];
	   }
	   model.add( flow_s_pos==1);

	   //C4-: sum{f[s][d][e(d,+)]}=0
	   IloExpr flow_d_pos(env);
	   
	   for (int p =0;p<Table[Vd].Degree; p=p+1)
	   {
		   int e=(Table[Vd].Adj)[p].Addr;
		   flow_d_pos+=fx[s][d][e];
	   }
	   model.add( flow_d_pos==0);


	   //C5+: sum{f[s][d][e(d,-)]}=1
	   IloExpr flow_d_neg(env);
	   for(int v=0;v<Nv;v++)for(int p =0;p<Table[v].Degree;p=p+1)
	   {	
		   if ((Table[v].Adj)[p].Dest==Vd)
		   {
				int e=(Table[v].Adj)[p].Addr;
				flow_d_neg+=fx[s][d][e];
		   }
	   }
	   model.add( flow_d_neg==1);

	   //C5-: sum{f[s][d][e(s,-)]}=0
	   IloExpr flow_s_neg(env);
	   for(int v=0;v<Nv;v++)for(int p =0;p<Table[v].Degree;p=p+1)
	   {	
		   if ((Table[v].Adj)[p].Dest==Vs)
		   {
				int e=(Table[v].Adj)[p].Addr;
				flow_s_neg+=fx[s][d][e];
		   }
	   }
	   model.add( flow_s_neg==0);

	   //C6: sum{f[s][d][e(u,+)]}=sum{f[s][d][e(u,-)]}
	   for (int v=0;v<Nv;v++)
	   {
		   if (v==Vs||v==Vd)
			   continue;
		   IloExpr flow_v_pos(env),flow_v_neg(env);
		   for (int p =0;p<Table[v].Degree; p=p+1)
		   {	
			   int e1=(Table[v].Adj)[p].Addr;
			   flow_v_pos+=fx[s][d][e1];
		   }
		   for (int w=0;w<Nv;w++)for (int q =0;q<Table[w].Degree; q=q+1)
		   {	
			   if ((Table[w].Adj)[q].Dest==v)
			   {
					int e2=(Table[w].Adj)[q].Addr;
					flow_v_neg+=fx[s][d][e2];
			   }
		   }
		   model.add(flow_v_pos==flow_v_neg);
	   }
	}
				   
	//	C7:prime-dual constraints on max utilization of links
	//		SUM sd{min(m[s],m[d])*x[e][s][d]}
	//			+SUM s{min(m[s],sumB-m[s])*y[e][s]}
	//			+SUM d{min(m[d],sumB-m[d])*z[e][s]} < u*c(e)
	//IloExpr bandwidth_cost(env);
	for(int e=0;e<Ne;e++)
	{	IloExpr mx(env);
		for (int s=0;s<n_location;s++)
		{	
			mx+=min(m[s],sumB-m[s])*(y[s][e]+z[s][e]);
			for (int d=0;d<n_location;d++)
				if (s!=d)
					mx+=min(m[s],m[d])*x[s][d][e];
		}
		//model.add(Bandwidth[e]-resBandwidth[e]+mx<=Bandwidth[e]*mlu);
		model.add(mx<=resBandwidth[e]*mlu);
		//bandwidth_cost+=mx;
		
	}
	//model.add(IloMinimize(env,bandwidth_cost));
	//	C8:prime-dual constraints:
	//		x[s][d][e]+y[s][e]+z[d][e] >= f[s][d][e]
	for(int e=0;e<Ne;e++)
		for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
			if(s!=d)
			   model.add( x[s][d][e]+y[s][e]+z[d][e] >= fx[s][d][e] );
			   //model.add(y[s][e]+z[d][e] >= fx[s][d][e] );

	//---------Cplex解方程---------	   
	IloCplex cplex(model); //依model建的cplex
	cplex.setOut(env.getNullStream());
	//cplex.setParam(IloCplex::RootAlg, IloCplex::Concurrent);
	if(cplex.solve())
	{	//cplex.exportModel("optRouting.lp");
		float tar =(float) cplex.getObjValue();
		float lp_mlu=0;

#ifndef _debug_LPRouting 

		for(int e=0;e<Ne;e++) for(int s=0;s<n_location;s++) for(int d=0;d<n_location;d++)
			f[the_server[s]][the_server[d]][e]=cplex.getValue(fx[s][d][e]);
		bool caculated[Ne]={0};
		for (int e=0;e<Ne;e++)
		{	
			if (caculated[e])
				continue;
			map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
			map.Bandwidth[opposite[e]]=map.Bandwidth[e];
			caculated[e]=true; caculated[opposite[e]]=true;
			//lp_mlu=max(lp_mlu,1-(resBandwidth[e]-map.Bandwidth[e])/Bandwidth[e]);
			lp_mlu=max(lp_mlu,map.Bandwidth[e]/resBandwidth[e]);
		}
			
#else
		if (tar<lb_mlu){
			for(int e=0;e<Ne;e++) for(int s=0;s<n_location;s++) for(int d=0;d<n_location;d++)
				f[the_server[s]][the_server[d]][e]=cplex.getValue(fx[s][d][e]);
			bool caculated[Ne]={0};
			for (int e=0;e<Ne;e++)
			{	
				if (caculated[e])
					continue;
				map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
				map.Bandwidth[opposite[e]]=map.Bandwidth[e];
				caculated[e]=true; caculated[opposite[e]]=true;

				lp_mlu=max(lp_mlu,map.Bandwidth[e]/resBandwidth[e]);
			}
			if (abs(tar-lp_mlu)>0.01)
				cout<<"\n"<<"the lp and max_traffic is not consistent!"<<"\n";

			cout<<"\n"<<"LP routing"<<"\n";
			for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
			{	if (s==d) continue;
				cout<<"\n"<<"f["<<the_server[s]<<"]["<<the_server[d]<<"]= ";
				for(int e=0;e<Ne;e++)
					if (f[the_server[s]][the_server[d]][e]>0)
					cout<<e<<" <"<<f[the_server[s]][the_server[d]][e]<<"> ";
			}

			
			float total_bandwidth_cost1=0
			for (int e=0;e<Ne;e++)
				total_bandwidth_cost1+=map.Bandwidth[e];		
		}
		else {
			cout<<"the solution is not optimal!\n";	
			cout<<"\n"<<"mlu(lb)="<<lb_mlu<<"; mlu(lp)="<<tar<<"\n";
		}
			
#endif
		env.end();
		return tar;

	}
	else
	{
		cout<<"No Solution!\n\n";
		//cplex.exportModel("optRouting.lp");
		env.end();
		return 2;
	}
}
float Graph::OptimalRouting(const Cluster& req,float* hostBw,Solution&map)
{
	
	//-------常量变量设置--------
	float sumB=0;
	for (int i=0;i<nServer;i++)
		sumB+=hostBw[i];
	
	//--------基本环境设置--------
	IloEnv env; //环境environment
	IloModel model(env); //model
		

	IloInt n_location=0; // number of used server
	for (int i=0;i<nServer;i++)
	   if (hostBw[i]>0)n_location++;
	if (n_location<=1){
		env.end();
		return 0;
	}
	IloNumArray m (env,n_location);
	IloIntArray the_server (env,n_location); //VM allocation
	for (int i=0,j=0;i<nServer;i++)
		if (hostBw[i]>0){
		   m[j]=hostBw[i]; the_server[j]=i;	j++;
		}
#ifdef _debug_LPRouting 
	float lb_mlu=0;
	// the mlu under load-balanced routing
	bool caculated[Ne]={0};
	for (int e=0;e<Ne;e++)
	{	
		if (caculated[e])
			continue;
		map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
		map.Bandwidth[opposite[e]]=map.Bandwidth[e];
		caculated[e]=true; caculated[opposite[e]]=true;

		lb_mlu=max(lb_mlu,1-(resBandwidth[e]-map.Bandwidth[e])/Bandwidth[e]);
	}
	/* display placement and routing information*/
	cout<<"\n"<<"load balance routing"<<"\n";
	for (int i=0;i<nServer;i++)
		cout<<hostBw[i]<<" ";
	for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
	{	if (s==d) continue;
		cout<<"\n"<<"f["<<the_server[s]<<"]["<<the_server[d]<<"]= ";
		for(int e=0;e<Ne;e++)
			if (f[the_server[s]][the_server[d]][e]>0)
			cout<<e<<" <"<<f[the_server[s]][the_server[d]][e]<<"> ";
	}

	if (lb_mlu<=0)
		return lb_mlu;
#endif
	//-------待定变量设置------ f[s][d][e]
	NumVarCubic fx(env, n_location); //
	for(int i=0;i<n_location;i++)
	{
	   fx[i] = NumVarMatrix(env,n_location);
	   for(int j=0;j<n_location;j++)
		   fx[i][j] = IloNumVarArray(env, Ne, 0, 1 );
	}
	// prime-dual variables
	/*NumVarCubic x(env, n_location); //x[s][d][e]
	for(int s=0;s<n_location;s++)
	{
	   x[s] = NumVarMatrix(env,n_location);  
	   for(int d=0;d<n_location;d++)
		   x[s][d] = IloNumVarArray(env, Ne, 0, IloInfinity );
	}*/
	NumVarMatrix y(env,n_location),z(env,n_location); //y[s][e],z[d][e]
	for(int s=0;s<n_location;s++)
	{
	   y[s] = IloNumVarArray(env, Ne, 0, IloInfinity );
	   z[s] = IloNumVarArray(env, Ne, 0, IloInfinity );
	}

//----------cplex列方程---------	
	//----优化目标----
	IloNumVar mlu (env);
	model.add(IloMinimize(env,mlu ));	
	// ---constraints-----//
	//model.add( mlu <= 1);

	//--flow constraints--//

	//C1:	f[s][s][e] =0
	for(int e=0;e<Ne;e++)  
		for(int s=0;s<n_location;s++){
		   model.add( fx[s][s][e] == 0);
		}
	/* force to routing on equal cost paths*/
	for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)for(int e=0;e<Ne;e++)  
	{
		int Vs=the_server[s];int Vd=the_server[d];
		if(b[Vs][Vd][e]==0)
			model.add(fx[s][d][e]==0);
	}

	//C2:	f[s][d][e(v,w)]=f[d][s][e(w,v)]
	for(int s=0;s<n_location;s++)for(int d=s+1;d<n_location;d++)for(int e=0;e<Ne;e++) 
	{	
		model.add(fx[s][d][e]=fx[d][s][opposite[e]]); //e1(v,w) and e2(w,v)
	}
			
	//flow conservation		
	for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
	{
	   if(s==d)
		   continue;

	   int Vs=the_server[s];int Vd=the_server[d];

	   //C4+: sum{f[s][d][e(s,+)]}=1
	   IloExpr flow_s_pos(env);
	   
	   for (int p =0;p<Table[Vs].Degree; p=p+1)
	   {
		   int e=(Table[Vs].Adj)[p].Addr;
		   flow_s_pos+=fx[s][d][e];
	   }
	   model.add( flow_s_pos==1);

	   //C4-: sum{f[s][d][e(d,+)]}=0
	   IloExpr flow_d_pos(env);
	   
	   for (int p =0;p<Table[Vd].Degree; p=p+1)
	   {
		   int e=(Table[Vd].Adj)[p].Addr;
		   flow_d_pos+=fx[s][d][e];
	   }
	   model.add( flow_d_pos==0);


	   //C5+: sum{f[s][d][e(d,-)]}=1
	   IloExpr flow_d_neg(env);
	   for(int v=0;v<Nv;v++)for(int p =0;p<Table[v].Degree;p=p+1)
	   {	
		   if ((Table[v].Adj)[p].Dest==Vd)
		   {
				int e=(Table[v].Adj)[p].Addr;
				flow_d_neg+=fx[s][d][e];
		   }
	   }
	   model.add( flow_d_neg==1);

	   //C5-: sum{f[s][d][e(s,-)]}=0
	   IloExpr flow_s_neg(env);
	   for(int v=0;v<Nv;v++)for(int p =0;p<Table[v].Degree;p=p+1)
	   {	
		   if ((Table[v].Adj)[p].Dest==Vs)
		   {
				int e=(Table[v].Adj)[p].Addr;
				flow_s_neg+=fx[s][d][e];
		   }
	   }
	   model.add( flow_s_neg==0);

	   //C6: sum{f[s][d][e(u,+)]}=sum{f[s][d][e(u,-)]}
	   for (int v=0;v<Nv;v++)
	   {
		   if (v==Vs||v==Vd)
			   continue;
		   IloExpr flow_v_pos(env),flow_v_neg(env);
		   for (int p =0;p<Table[v].Degree; p=p+1)
		   {	
			   int e1=(Table[v].Adj)[p].Addr;
			   flow_v_pos+=fx[s][d][e1];
		   }
		   for (int w=0;w<Nv;w++)for (int q =0;q<Table[w].Degree; q=q+1)
		   {	
			   if ((Table[w].Adj)[q].Dest==v)
			   {
					int e2=(Table[w].Adj)[q].Addr;
					flow_v_neg+=fx[s][d][e2];
			   }
		   }
		   model.add(flow_v_pos==flow_v_neg);
	   }
	}
				   
	//	C7:prime-dual constraints on max utilization of links
	//		SUM sd{min(m[s],m[d])*x[e][s][d]}
	//			+SUM s{min(m[s],sumB-m[s])*y[e][s]}
	//			+SUM d{min(m[d],sumB-m[d])*z[e][s]} < u*c(e)
	//IloExpr bandwidth_cost(env);
	for(int e=0;e<Ne;e++)
	{	IloExpr mx(env);
		for (int s=0;s<n_location;s++)
		{	
			mx+=min(m[s],sumB-m[s])*(y[s][e]+z[s][e]);
			/*for (int d=0;d<n_location;d++)
				if (s!=d)
					mx+=min(m[s],m[d])*x[s][d][e];*/
		}
		//model.add(Bandwidth[e]-resBandwidth[e]+mx<=Bandwidth[e]*mlu);
		model.add(mx<=resBandwidth[e]*mlu);
		//bandwidth_cost+=mx;
		
	}
	//model.add(IloMinimize(env,bandwidth_cost));
	//	C8:prime-dual constraints:
	//		x[s][d][e]+y[s][e]+z[d][e] >= f[s][d][e]
	for(int e=0;e<Ne;e++)
		for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
			if(s!=d)
				model.add(y[s][e]+z[d][e] >= fx[s][d][e] );
			//if(s!=d)
			   //model.add( x[s][d][e]+y[s][e]+z[d][e] >= fx[s][d][e] );
			   //model.add(y[s][e]+z[d][e] >= fx[s][d][e] );

	//---------Cplex解方程---------	   
	IloCplex cplex(model); //依model建的cplex
	cplex.setOut(env.getNullStream());
	//cplex.setParam(IloCplex::RootAlg, IloCplex::Concurrent);
	if(cplex.solve())
	{	//cplex.exportModel("optRouting.lp");
		float tar =(float) cplex.getObjValue();
		
#ifdef _debug_LPRouting 
		float lp_mlu=0;			
		if (tar<lb_mlu){
			for(int e=0;e<Ne;e++) for(int s=0;s<n_location;s++) for(int d=0;d<n_location;d++)
				f[the_server[s]][the_server[d]][e]=cplex.getValue(fx[s][d][e]);
			bool caculated[Ne]={0};
			for (int e=0;e<Ne;e++)
			{	
				if (caculated[e])
					continue;
				map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
				map.Bandwidth[opposite[e]]=map.Bandwidth[e];
				caculated[e]=true; caculated[opposite[e]]=true;

				lp_mlu=max(lp_mlu,map.Bandwidth[e]/resBandwidth[e]);
			}
			if (abs(tar-lp_mlu)>0.01)
				cout<<"\n"<<"the lp and max_traffic is not consistent!"<<"\n";

			cout<<"\n"<<"LP routing"<<"\n";
			for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
			{	if (s==d) continue;
				cout<<"\n"<<"f["<<the_server[s]<<"]["<<the_server[d]<<"]= ";
				for(int e=0;e<Ne;e++)
					if (f[the_server[s]][the_server[d]][e]>0)
					cout<<e<<" <"<<f[the_server[s]][the_server[d]][e]<<"> ";
			}

			
			float total_bandwidth_cost1=0
			for (int e=0;e<Ne;e++)
				total_bandwidth_cost1+=map.Bandwidth[e];		
		}
		else {
			cout<<"the solution is not optimal!\n";	
			cout<<"\n"<<"mlu(lb)="<<lb_mlu<<"; mlu(lp)="<<tar<<"\n";
		}
			
#endif
//caculate the max load for each link
		//enforce the optimal routing
		for(int e=0;e<Ne;e++) for(int s=0;s<n_location;s++) for(int d=0;d<n_location;d++)
				f[the_server[s]][the_server[d]][e]=cplex.getValue(fx[s][d][e]);

		for(int e=0;e<Ne;e++)
		{	map.Bandwidth[e]=0;
			
			for (int s=0;s<n_location;s++)
			{	
				map.Bandwidth[e]+=min(m[s],sumB-m[s])*(cplex.getValue(y[s][e]+z[s][e]));
			}
			if(map.Bandwidth[e]/resBandwidth[e]-tar>min_error)
				cout<<"found a link with larger mlu\n";
			
		}

		env.end();
		return tar;

	}
	else
	{
		cout<<"No Solution!\n\n";
		//cplex.exportModel("optRouting.lp");
		env.end();
		return 2;
	}
}

float Graph::CalcMaxLinkUt(const Cluster& req,float* hostBw,Solution&map)
{
	
	//-------常量变量设置--------
	float sumB=0;
	for (int i=0;i<nServer;i++)
		sumB+=hostBw[i];
	
	//--------基本环境设置--------
	IloEnv env; //环境environment
	IloModel model(env); //model
		

	IloInt n_location=0; // number of used server
	for (int i=0;i<nServer;i++)
	   if (hostBw[i]>0)n_location++;
	if (n_location<=1){
		env.end();
		return 0;
	}
	IloNumArray m (env,n_location);
	IloIntArray the_server (env,n_location); //VM allocation
	for (int i=0,j=0;i<nServer;i++)
		if (hostBw[i]>0){
		   m[j]=hostBw[i]; the_server[j]=i;	j++;
		}

	// prime-dual variables
	/*NumVarCubic x(env, n_location); //x[s][d][e]
	for(int s=0;s<n_location;s++)
	{
	   x[s] = NumVarMatrix(env,n_location);  
	   for(int d=0;d<n_location;d++)
		   x[s][d] = IloNumVarArray(env, Ne, 0, IloInfinity );
	}*/
	NumVarMatrix y(env,n_location),z(env,n_location); //y[s][e],z[d][e]
	for(int s=0;s<n_location;s++)
	{
	   y[s] = IloNumVarArray(env, Ne, 0, IloInfinity );
	   z[s] = IloNumVarArray(env, Ne, 0, IloInfinity );
	}

//----------cplex列方程---------	
	//----优化目标----
	IloNumVar mlu (env);
	model.add(IloMinimize(env,mlu ));		
				   
	//	C7:prime-dual constraints on max utilization of links
	//		SUM sd{min(m[s],m[d])*x[e][s][d]}
	//			+SUM s{min(m[s],sumB-m[s])*y[e][s]}
	//			+SUM d{min(m[d],sumB-m[d])*z[e][s]} < u*c(e)
	//IloExpr bandwidth_cost(env);
	for(int e=0;e<Ne;e++)
	{	
		
		IloExpr mx(env);
		for (int s=0;s<n_location;s++)
		{	
			mx+=min(m[s],sumB-m[s])*(y[s][e]+z[s][e]);
			/*for (int d=0;d<n_location;d++)
				if (s!=d)
					mx+=min(m[s],m[d])*x[s][d][e];*/
		}
		//model.add(Bandwidth[e]-resBandwidth[e]+mx<=Bandwidth[e]*mlu);
		model.add(mx<=resBandwidth[e]*mlu);
		//bandwidth_cost+=mx;
		
	}
	//model.add(IloMinimize(env,bandwidth_cost));
	//	C8:prime-dual constraints:
	//		x[s][d][e]+y[s][e]+z[d][e] >= f[s][d][e]
	for(int e=0;e<Ne;e++)
		for(int s=0;s<n_location;s++)for(int d=0;d<n_location;d++)
			if(s!=d)
				model.add(y[s][e]+z[d][e] >= f[the_server[s]][the_server[d]][e] );
			//if(s!=d)
			   //model.add( x[s][d][e]+y[s][e]+z[d][e] >= f[s][d][e] );
			   //model.add(y[s][e]+z[d][e] >= f[s][d][e] );

	//---------Cplex解方程---------	   
	IloCplex cplex(model); //依model建的cplex
	cplex.setOut(env.getNullStream());
	//cplex.setParam(IloCplex::RootAlg, IloCplex::Concurrent);
	if(cplex.solve())
	{	//cplex.exportModel("optRouting.lp");
		float tar =(float) cplex.getObjValue();

		//enforce the optimal routing

		for(int e=0;e<Ne;e++)
		{	map.Bandwidth[e]=0;
			//caculate the max load for each link
			for (int s=0;s<n_location;s++)
			{	
				map.Bandwidth[e]+=min(m[s],sumB-m[s])*(cplex.getValue(y[s][e]+z[s][e]));
			}
			if(map.Bandwidth[e]/resBandwidth[e]-tar>min_error)
				cout<<"found a link with larger mlu\n";
			
		}

		env.end();
		return tar;

	}
	else
	{
		cout<<"No Solution!\n\n";
		//cplex.exportModel("optRouting.lp");
		env.end();
		return 2;
	}
}

bool Graph::QuickFail(Cluster& req,Solution&map,float(*sum_capacity)[nServer],float& sumB,float* res_port_B)
{
	//sort B1,...BN in descending order
	sort(req.B.begin(),req.B.end(),greater<float>());
	//intiate the map variables:
	map.Clear();	
	//load balance routing 
	for (int s=0;s<nServer;s++)
		for (int d=s+1;d<nServer;d++)
		{   
			if(routingOption==_KshortestLB)
				sum_capacity[s][d]=LoadBalance(s,d);
			else if(routingOption==_ECMP)
				sum_capacity[s][d]=ECMP(s,d);
			sum_capacity[d][s]=sum_capacity[s][d];
			for (int e=0;e<Ne;e++)
				f[d][s][opposite[e]]=f[s][d][e];
		}
	// pre-sorting the pairs(s,d) through e
	for (int e=0;e<Ne;e++)
		n_pair[e]=Sort_fsd(e);

	// check residual slot and port bandwidths	
	bool resource_satisfied=true;
	float res_sumB=0;
	int res_slot=0;
	for (int i=0;i<req.N;i++)// total suscribed bandwidth
		sumB+=req.B[i];
	for(int i=0;i<nServer;i++){
		for (int p =0;p<Table[i].Degree; p=p+1)
		{
			int e=(Table[i].Adj)[p].Addr;
			res_port_B[i]+=resBandwidth[e];
		}
				
		if (Table[i].Slot>0&&res_port_B[i]>0){
			res_slot+=Table[i].Slot;
			res_sumB+=res_port_B[i];
		}
	}
	if (res_slot<req.N||res_sumB<sumB)
		resource_satisfied=false;
	
	return resource_satisfied;
	//end of quick fail !!	
	
}
bool Graph::CongestDetect(int x,const Cluster& req,float* hostBw,float (*sum_capacity)[nServer],Solution& map)
{		
	bool congest=false;
	bool caculated[Ne]={1};
	
	//step1:congestion test for server-to-server traffic
#ifdef _server2server
	for (int y=0;y<nServer;y++)
	{
		if (map.Slot[y]<=0||x==y)// no VM placed on server y
			continue; 
		// which edges has increased congestion?
		for (int e=0;e<Ne;e++) 
			if (b[x][y][e]||b[y][x][e])
				caculated[e]=0;
		if (sum_capacity[x][y]<min(hostBw[x],hostBw[y]))//capacity<load on paths between server x and  
		{	congest=true;
			break;// can not place VM(Bi) to server x,try x+1!
		}
	}
	if (congest)	//routing failure for server to server traffic
		return true;	//try the next server x+1 for VM(Bi)!!
#endif
	//step2: congestion test for concurrent traffic
	//Solution oldmap(map);
	//edge ranking according to load
	
#ifdef _edge_rank			
	int _rank[Ne];	float _load[Ne];
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
		map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
		map.Bandwidth[opposite[e]]=map.Bandwidth[e];
		caculated[e]=true; caculated[opposite[e]]=true;
		if (map.Bandwidth[e]>resBandwidth[e])
		{	//congestion detected!
			congest=true;
			//map=oldmap;
			break;
		}		
		
	}	//end for: congestion test for concurrent traffic
#else
	//whether e(v,w) is congested?
	for (int e=0;e<Ne;e++)
	{	
		map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
		if (map.Bandwidth[e]>resBandwidth[e])
		{	//congestion detected!
			congest=true;	
			map=oldmap;
			break;
		}		
	}	//end for: congestion test for concurrent traffic

#endif		
	return congest;

}

float Graph::assignBandwidth(Cluster& req,float* hostBw,Solution& map)
{		
	float mlu=0;
	bool caculated[Ne]={0};

			
	int _rank[Ne];	float _load[Ne];
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
		map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
		map.Bandwidth[opposite[e]]=map.Bandwidth[e];
		caculated[e]=true; caculated[opposite[e]]=true;
		mlu=max(mlu,map.Bandwidth[e]/resBandwidth[e]);
		if (mlu-1>min_error)	//congestion detected!
			break;
		
	}	//end for: congestion test for concurrent traffic

	return mlu;

}


// pertubation only when a VM can not placed
bool Graph::Pertubation(Cluster& req,bool enLProuting,Solution& map)
{	
	float sum_capacity[nServer][nServer];
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	vector<vector<bool>> tabu(req.N,vector<bool>(nServer,false));
	float res_port_B[nServer]={0},sumB=0;	//quick fail variables
	if (QuickFail(req,map,sum_capacity,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!	
	
	int numVMembedded=0,maxLoop=req.N*nServer;
	//int n_bottle[nServer]={0};
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
			
			bool congest=false;
			if(enLProuting)
				congest=OptimalRouting(req,hostBw,map)>1;
			else
				congest=CongestDetect(x,req,hostBw,sum_capacity,map);
			
			if (congest){
				// fail due to conestion,try pertubation
				map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;	
				updateBottlenecks(hostBw,map,TS,TD);
			}
			else	
				break;
			//else, try server x+1
		}	//end for server x
		
		if (numVMembedded==req.N)
		{	
			//if (enLProuting)
			//	LPRouting(req,hostBw,map);
			if (enLProuting){
				float mlu=assignBandwidth(req,hostBw,map);
				if(mlu-1>min_error)
					cout<<"error in Optimal Routing/n";
			}
			return true;
		}
		else if (x>=nServer)
		{
			if(n_usefull_server==0)	 // running out of VM slots!
				return false;
			// fail due to conestion,try pertubation						
			int the_server=-1,the_vm;
			// find the bottleneck server to lead to congestion
			//for(int j=0;j<nServer;j++)
			//	if (n_bottle[j]>max_bottle){
			//		max_bottle=n_bottle[j];
			//		the_server=j;
			//	}
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


// random drop a VM
bool Graph::randomDrop(Cluster& req,bool enLProuting,Solution& map)
{	
	float sum_capacity[nServer][nServer];
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	vector<vector<bool>> tabu(req.N,vector<bool>(nServer,false));
	float res_port_B[nServer]={0},sumB=0;	//quick fail variables
	if (QuickFail(req,map,sum_capacity,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!	
	
	int numVMembedded=0,maxLoop=req.N*nServer;
	//int n_bottle[nServer]={0};
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
			if (CongestDetect(x,req,hostBw,sum_capacity,map)){
				// fail due to conestion,try pertubation
				map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;	
				//int the_server;
				//// find the bottleneck server to lead to congestion
				//the_server=findBottleneck(hostBw,map);
				//if (the_server<0||the_server>=nServer) 
				//	return false; //fail due to bisection traffic, rather than concurrent congestion.
				//else
				//	n_bottle[the_server]++;
				//updateBottlenecks(hostBw,map,TS,TD);
			}
			else	
				break;
			//else, try server x+1
		}	//end for server x

		if (numVMembedded==req.N)
		{	
			if (enLProuting)
			LPRouting(req,hostBw,map);
			return true;
		}
		if (x>=nServer)
		{
			if(n_usefull_server==0)	 // running out of VM slots!
				return false;
			// fail due to conestion,try pertubation						
			int the_server=-1,the_vm;
			// find a random server that have hosted a VM
			int first_server=unif_int(0,nServer-1);
			for(int j=0;j<nServer;j++)
			{
				the_server=(first_server+j)%nServer;
				if (hostBw[the_server]>0)
					break;
			}
			
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
bool Graph::FirstFit(Cluster& req,bool enLProuting,Solution& map)
{	
	float sum_capacity[nServer][nServer];
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
		 
	float res_port_B[nServer]={0},sumB=0;//quick fail variables
	if (QuickFail(req,map,sum_capacity,sumB,res_port_B)==false)
		return false;
	//i: visit the VMs, j: visit the servers
	int numVMembedded=0;
	for (int i=0;i<req.N;i++)
	{	int x;
		for (x=0;x<nServer;x++)
		{	// packing VMs( i=1:N )into server x
			
			if (Table[x].Slot<=map.Slot[x])// there is no space in server x
				continue;	//open the next server 
			if(res_port_B[x]<min(hostBw[x]+req.B[i],sumB-hostBw[x]-req.B[i]))
				continue;	 //less ingress bandwidth!!  			
			map.Slot[x]+=1;hostBw[x]+=req.B[i];assignment[i]=x;numVMembedded++;
			if (CongestDetect(x,req,hostBw,sum_capacity,map))
			{	//undo the placement and try the next server
				map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;
			}
			else{
				if (numVMembedded==req.N)	
				{	// already found the host for the last VM
					if(enLProuting)
						LPRouting(req,hostBw,map);
					return true;
				}
				else // already found the host for the current VM
					break;
			}			
		}	//end for server x
		if (x>=nServer) // fail to find a suitable host
			return false;
	}	// end for VM(Bi)

	
}


bool Graph::NextFit(Cluster& req,bool enLProuting,Solution& map)
{	
	float sum_capacity[nServer][nServer];
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
		 
	float res_port_B[nServer]={0},sumB=0;//quick fail variables
	if (QuickFail(req,map,sum_capacity,sumB,res_port_B)==false)
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
			if (CongestDetect(x,req,hostBw,sum_capacity,map))
			{	//undo the placement and try the next server
				map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;
			}
			else{
				if (numVMembedded==req.N)	
				{	// already found the host for the last VM
					if(enLProuting)
						LPRouting(req,hostBw,map);
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


bool Graph::BestFit(Cluster& req,bool enLProuting,Solution& map)
{	
	float sum_capacity[nServer][nServer];
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
		 
	float res_port_B[nServer]={0},sumB=0;//quick fail variables
	if (QuickFail(req,map,sum_capacity,sumB,res_port_B)==false)
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
			CongestDetect(x,req,hostBw,sum_capacity,map);
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
			if(enLProuting)
				LPRouting(req,hostBw,map);
			return true;
		}

	}	// end for VM(Bi)	
}


// exhaustive search with load-balanced routing
bool Graph::BackTracking(Cluster &req,bool enLProuting, Solution &map)
{		
	float sum_capacity[nServer][nServer];
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	
	float res_port_B[nServer]={0},sumB=0; //quick fail variables
	if (QuickFail(req,map,sum_capacity,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!	

	//i: visit the VMs, j: visit the servers
	int numVMembedded=0;
	if (recursivePlacement(req,enLProuting,map,assignment,
							sumB,res_port_B,sum_capacity,hostBw,numVMembedded))	
		return true;
	else
		return false;


}
// used in backtracking algorihtm
bool Graph::recursivePlacement(Cluster &req,bool enLProuting,Solution &map,vector<int>& assignment,
		float sumB,float* res_port_B,float (*sum_capacity)[nServer],float* hostBw,int &numVMembedded)
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
		bool congest=false;
		if (enLProuting)
			congest=LPRouting(req,hostBw,map)>1;
		else
			congest=CongestDetect(x,req,hostBw,sum_capacity,map);
			
		if (congest){
			//routing failure,try the next server x+1 for VM(Bi)!!
			map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;	
		}
		else //not congestion
			if (numVMembedded==req.N)
				return true;
			else  //place the next VMs recursively
				if (recursivePlacement(req,enLProuting,map,assignment,
										sumB,res_port_B,sum_capacity,hostBw,numVMembedded))
					return true;
				else // back tracking!!
				{	map.Slot[x]-=1;	hostBw[x]-=req.B[i];	assignment[i]=-1; numVMembedded--;
					continue;//try the next server x+1 for VM(Bi)!!
				}

	}//end for server x
	if (x>=nServer)	
	return false;	

}




//-------code histroy-----------------------//
#ifdef _enumeration
bool Graph::ExhaustiveSearch(Cluster &req,bool enLProuting, Solution &map)
{
	sort(req.B.begin(),req.B.end(),greater<float>());	// sort B1,...BN in descending order	
	float sum_capacity[nServer][nServer];
	
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	
	//intiate the map variables:
	map.Clear();
	for (int s=0;s<nServer;s++)
		for (int d=s+1;d<nServer;d++)
		{   
			sum_capacity[s][d]=LoadBalance(s,d);
			sum_capacity[d][s]=sum_capacity[s][d];
			for (int e=0;e<Ne;e++)
				f[d][s][opposite[e]]=f[s][d][e];
		}

	for (int e=0;e<Ne;e++)
		n_pair[e]=Sort_fsd(e);
		//quick fail 
	float res_port_B[nServer]={0},sumB=0;	
	if (QuickFail(req,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!	


	//i: visit the VMs, j: visit the servers
	int numVMembedded=0;
	if (Enumeration(req, map,assignment,sumB,res_port_B,
		hostBw,numVMembedded))
	{
		return true;
	}
	else
		return false;

}
bool Graph::Enumeration(Cluster &req, Solution &map,vector<int>& assignment,float sumB,float* res_port_B,							   float* hostBw,int &numVMembedded)
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

			if (numVMembedded==req.N){
				if (LPRouting(req,hostBw,map)<=1)			
					return true;
				else // back tracking!!
				{	
					map.Slot[x]-=1;	hostBw[x]-=req.B[i];	assignment[i]=-1; numVMembedded--;
					continue;//try the next server x+1 for VM(Bi)!!

				}
				
			}
			else  //place the next VMs recursively
				if (Enumeration(req, map,assignment,sumB,res_port_B,
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


#endif
#ifndef _newBottleckCriteria
// pertubation only when a VM can not placed
bool Graph::Pertubation(Cluster& req,bool enLProuting,Solution& map)
{	
	float sum_capacity[nServer][nServer];
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in
	vector<vector<bool>> tabu(req.N,vector<bool>(nServer,false));
	float res_port_B[nServer]={0},sumB=0;	//quick fail variables
	if (QuickFail(req,map,sum_capacity,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!	
	
	int numVMembedded=0,maxLoop=req.N*nServer;
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
			if (CongestDetect(x,req,hostBw,sum_capacity,map)){
				map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;	}
			else	
				break;
			//else, try server x+1
		}	//end for server x

		if (numVMembedded==req.N)
		{	
			if (enLProuting)
			LPRouting(req,hostBw,map);
			return true;
		}
		if (x>=nServer)
		{
			if(n_usefull_server==0)	 // running out of VM slots!
				return false;
			// fail due to conestion,try pertubation						
			int the_server,the_vm;
			// find the bottleneck server to lead to congestion
			the_server=findBottleneck(hostBw,map);
			if (the_server<0) 
				return false; //fail due to bisection traffic, rather than concurrent congestion.
			for(the_vm=req.N-1;the_vm>=0;the_vm--)
				if (assignment[the_vm]==the_server)
					break;	
			//unload a VM from the bottleneck server
			map.Slot[the_server]-=1;	hostBw[the_server]-=req.B[the_vm];
			assignment[the_vm]=-1;		numVMembedded--;
			tabu[the_vm][the_server]=1;
		}
		//else, Let's place the next VM(Bi+1)
	}
	return false;
}


#endif
#ifdef _local_search
bool Graph::LocalSearch(Cluster& req,bool enLProuting,Solution& map)
{	
	sort(req.B.begin(),req.B.end(),greater<float>());	// sort B1,...BN in descending order	
	float sum_capacity[nServer][nServer];
	
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in

	vector<vector<bool>> tabu(req.N,vector<bool>(nServer,false));
	//intiate the map variables:
	map.Clear();
	for (int s=0;s<nServer;s++)
		for (int d=s+1;d<nServer;d++)
		{   
			sum_capacity[s][d]=LoadBalance(s,d);
			sum_capacity[d][s]=sum_capacity[s][d];
			for (int e=0;e<Ne;e++)
				f[d][s][opposite[e]]=f[s][d][e];
		}

	for (int e=0;e<Ne;e++)
		n_pair[e]=Sort_fsd(e);
	
		//quick fail 
	float res_port_B[nServer]={0},sumB=0;	
	if (QuickFail(req,map,sum_capacity,sumB,res_port_B)==false)
		return false;
	//end of quick fail !!	
	int numVMembedded=0,maxLoop=req.N*nServer;
	
	/*----------initial solution--------------*/	
	for (int i=0;i<req.N;i++)
	{
		int x=0;
		//to place VM(Bi)into the servers, {0,1,...tabu[i]} are not considered
		for (;x<nServer;x++)
		{	
			if (Table[x].Slot<=map.Slot[x])// there is no space in server m
				continue;	//try the next server for VM(Bi)
			if(res_port_B[x]<min(hostBw[x]+req.B[i],sumB-hostBw[x]-req.B[i]))
				continue;	 //less ingress bandwidth!! 
			map.Slot[x]+=1;hostBw[x]+=req.B[i];assignment[i]=x;numVMembedded++;
			break; // next VM
		}
		if(x>=nServer)
			return false;
	}
	// is the initial solution feasible?
	if(CongestDetect(req,hostBw,map)==false)//function code have changed
	{
		if (enLProuting)
		LPRouting(req,hostBw,map);
		return true; 
	}
	else while(1){
		float mlu=0;	int the_link=0;
		for (int e=0;e<Ne;e++){
			if (mlu<map.Bandwidth[e]){
				mlu=map.Bandwidth[e];
				the_link=e; //find the most congested link
			}
		}
		if (map.Bandwidth[the_link]<=resBandwidth[the_link])
			break; // congestion elimated!!
		int the_server=-1,the_vm=-1;
		//find the bottleneck server to lead to congestion
		the_server=findBottleneck(the_link,hostBw);
		if (the_server<0) 
			return false; //fail due to bisection traffic, rather than concurrent congestion.
		for(the_vm=req.N-1;the_vm>=0;the_vm--)
			if (assignment[the_vm]==the_server)
				break;	
		//unload a VM from the bottleneck server
		map.Slot[the_server]-=1;	hostBw[the_server]-=req.B[the_vm];
		assignment[the_vm]=-1;		numVMembedded--;	tabu[the_vm][the_server]=1;
		// caculate the new max_load
		bool re_caculated[Ne]={0};
		for (int e=0;e<Ne;e++)
			if(re_caculated[e]==false&&map.Bandwidth[e]>resBandwidth[e]){
				map.Bandwidth[e]=LPmaxTraffic(e,hostBw);
				map.Bandwidth[opposite[e]]=map.Bandwidth[e];
				re_caculated[e]=true; re_caculated[opposite[e]]=true;
			}
	}
	// continue embedding using pertubation algorithm
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
			if (CongestDetect(x,req,hostBw,sum_capacity,map)){
				map.Slot[x]-=1; hostBw[x]-=req.B[i]; assignment[i]=-1; numVMembedded--;	}
			else	
				break;
			//else, try server x+1
		}	//end for server x

		if (numVMembedded==req.N)
		{	
			if (enLProuting)
			LPRouting(req,hostBw,map);
			return true;
		}
		if (x>=nServer)
		{
			if(n_usefull_server==0)	 // running out of VM slots!
				return false;
			// fail due to conestion,try pertubation
			int the_server=-1,the_vm=-1;
			//find the bottleneck server to lead to congestion
			the_server=findBottleneck(hostBw,map);
			if (the_server<0) 
				return false; //fail due to bisection traffic, rather than concurrent congestion.
			for(the_vm=req.N-1;the_vm>=0;the_vm--)
				if (assignment[the_vm]==the_server)
					break;	
			//unload a VM from the bottleneck server
			map.Slot[the_server]-=1;	hostBw[the_server]-=req.B[the_vm];
			assignment[the_vm]=-1;		numVMembedded--;
			tabu[the_vm][the_server]=1;
		}
		//else, Let's place the next VM(Bi+1)
	}
	return false;
	
}

#endif
#ifdef _cut_check
bool Graph::CutCheck(const Cluster &req, float &sumB, float *hostBw)
{
	bool congest=false;
	for (int i=0;i<nToR;i++)
	{	int the_ToR=nServer+i,e,x;
		float traffic=0,capacity=0;
		for (int j=0;j<nServerInRack;j++){
			int the_server=i*nServerInRack+j;
			traffic+=hostBw[the_server];
		}
		traffic=min(traffic,sumB-traffic);

		for(int p=0;p<Table[the_ToR].Degree;p++){
			x=(Table[the_ToR].Adj)[p].Dest;
			e=(Table[the_ToR].Adj)[p].Addr;
			if (x>=nServer+nToR) // it's a up-layer switch
				capacity+=resBandwidth[e];
		}
		if(traffic>capacity)
		{
			congest=true; 
			break;
		}
	}
	return congest;
}
#endif
// pertubation whenever congestion is detected
#ifdef _pert_once_congest
bool Graph::Pertubation(Cluster& req,Solution& map)
{	
	sort(req.B.begin(),req.B.end(),greater<float>());	// sort B1,...BN in descending order	
	float sum_capacity[nServer][nServer];
	
	float hostBw[nServer]={0};	//total bandwidths of a server assigned to the current VDC 
	vector<int> assignment(req.N,-1); //the server index  B1,...BN located in

	vector<vector<bool>> tabu(req.N,vector<bool>(nServer,false));
	//intiate the map variables:
	map.Clear();
	for (int s=0;s<nServer;s++)
		for (int d=0;d<nServer;d++)
			sum_capacity[s][d]=LoadBalance(s,d);
	for (int e=0;e<Ne;e++)
		n_pair[e]=Sort_fsd(e);
		
	//i: visit the VMs, j: visit the servers
	int numVMembedded=0,maxLoop=req.N*nServer;
	
	for(int t=0;t<maxLoop+1;t++)
	{	
		// find VM(Bi) with max B, and try to place to server x
		int i,x;
		for (i=0;i<req.N;i++)
			if (assignment[i]==-1)
				break;
		//to place VM(Bi)into the servers, {0,1,...tabu[i]} are not considered
		for (x=0;x<nServer;x++)
		{	
			if (tabu[i][x])// this server x is in the tabu list of VM(Bi)
				continue;
			if (Table[x].Slot<=map.Slot[x])// there is no space in server m
				continue;	//try the next server for VM(Bi)
			bool congest=false;
			//step1:congestion test for server-to-server traffic
#ifdef server2server
				for (int y=0;y<nServer;y++)
				{
					if (map.Slot[y]<=0||x==y)// no VM placed on server y
						continue;
					if (sum_capacity[x][y]<min(hostBw[x]+req.B[i],hostBw[y]))//capacity<load on paths between server x and  
					{	congest=true;
						break;// can not place VM(Bi) to server x,try x+1!
					}
				}
				if (congest)	//routing failure for server to server traffic
					continue;	//try the next server x+1 for VM(Bi)!!
#endif
			//step2: congestion test for concurrent traffic
			map.Slot[x]+=1;hostBw[x]+=req.B[i];assignment[i]=x;numVMembedded++;
			Pair bottleneck; 
			for (int e=0;e<Ne;e++)
				map.Bandwidth[e]=0;
			for (int v=0;v<Nv;v++)	for (int p =0;p<Table[v].Degree; p=p+1)
			{	//whether e(v,w) is congested?
				int e=(Table[v].Adj)[p].Addr;
				do{
					map.Bandwidth[e]=maxTraffic(e,hostBw);
					if (map.Bandwidth[e]>resBandwidth[e])
					{	//congestion detected!
						congest=true;
						//try pertubation
						int the_server=-1,the_vm=-1;
						//find the bottleneck server to lead to congestion
						the_server=findBottleneck(e,hostBw);
						for(the_vm=req.N-1;the_vm>=0;the_vm--)
							if (assignment[the_vm]==the_server)
								break;								
						map.Slot[the_server]-=1;	hostBw[the_server]-=req.B[the_vm];															
						assignment[the_vm]=-1;tabu[the_vm][the_server]=1;
						numVMembedded--;
					}
					else
						congest=false;
				}while(congest);
			}	//end for: congestion test for concurrent traffic
			break;
		}	//end for server x

		if (numVMembedded==req.N)
		{
			//LocalSearch(req,tabu,map,assignment,hostBw);
			return true;
		}
		if (x>=nServer)
			return false;	
		//else, Let's place the next VM(Bi+1)
	}

}
#endif
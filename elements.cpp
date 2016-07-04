#include "elements.h"
#include <limits>      // Used for infinity (for integers)
#include <iostream>



// to define the Virutal data center
Cluster::Cluster(const Cluster& Rhs)
{
	N=Rhs.N;
	Arrivaltime=Rhs.Arrivaltime;
	Holdtime=Rhs.Holdtime;
	for(int i=0;i<N;i++)
		B[i]=Rhs.B[i];
}
void Cluster::random()
{
	N=unif_int(minN,maxN);
	B.resize(N);
	for (int i=0;i<N;i++)
		B[i]=unif_int(minB,maxB);
	Arrivaltime=exprnd(muArrivaltime);
	Holdtime=exprnd(muHoldtime);
}
#ifndef _Tree
void Cluster::random(float load)
{
	N=unif_int(minN,maxN);
	B.resize(N);
	for (int i=0;i<N;i++)
		B[i]=unif_int(minB,maxB)/100*100;
	//muArrivalTime=0.5*(minN+maxN)*muHoldtime/(load*nServer*maxSlot);
	Holdtime=exprnd(muHoldtime);
	float vm_load=0.5*(minN+maxN)/(nServer*maxSlot),
		bw_load=0.25*(minN+maxN)*(minB+maxB)/(nServer*gbpsCommodity*n_server_port);
	Arrivaltime=exprnd(max(vm_load,bw_load)*(muHoldtime/load));
}
#else
void Cluster::random(float load)
{
	N=unif_int(minN,maxN);
	B.resize(N);
	for (int i=0;i<N;i++)
		B[i]=unif_int(minB,maxB)/10*10;
	//muArrivalTime=0.5*(minN+maxN)*muHoldtime/(load*nServer*maxSlot);
	Holdtime=exprnd(muHoldtime);
	Arrivaltime=exprnd(0.5*(minN+maxN)*muHoldtime/(load*nServer*maxSlot));
}
#endif
void Cluster::random(int numOfVM,float minBw,float maxBw)
{
	N=numOfVM;
	B.resize(N);
	for (int i=0;i<N;i++)
		B[i]=unif_int(minBw,maxBw)/10*10;
}
void Cluster::random(int numOfVM1,int numOfVM2,float minBw,float maxBw)
{
	N=unif_int(numOfVM1,numOfVM2);
	B.resize(N);
	for (int i=0;i<N;i++)
		B[i]=unif_int(minBw,maxBw)/10*10;
}
// to store the embedding solution
Solution::Solution()
{
	Departtime=0;
	for (int i=0;i<Ne;i++)
		Bandwidth[i]=0;
	for (int i=0;i<nServer;i++)
		Slot[i]=0;
}
Solution::Solution(const Solution&Rhs)
{
	for (int i=0;i<Ne;i++)
		Bandwidth[i]=Rhs.Bandwidth[i];
	for (int i=0;i<nServer;i++)
		Slot[i]=Rhs.Slot[i];
	Departtime=Rhs.Departtime;
}
void Solution::Clear()
{
	Departtime=0;
	for (int i=0;i<Ne;i++)
		Bandwidth[i]=0;
	for (int i=0;i<nServer;i++)
		Slot[i]=0;
	
}
const Solution& Solution::operator =(const Solution &Rhs) 
{
	Departtime=Rhs.Departtime;
	for (int i=0;i<Ne;i++)
		Bandwidth[i]=Rhs.Bandwidth[i];
	for (int i=0;i<nServer;i++)
		Slot[i]=Rhs.Slot[i];
	
	return *this;
	
}
void Performance::clear()
{
	sucess_rate=0;max_utilization=0;RC=0;
}

void allocateRange::resize(int nVM,int nSlot)
{
	n=nVM;
	a=nSlot;

	p=new bool **[n];
	//split=new int *[n]; 
	for(int x=0;x<n;x++)
	{	//x=0,1,...N-1
		p[x]=new bool *[n-x]; //y=x,...N-1
		//split[x]=new int[n-x];
		
		for(int j=0;j<n-x;j++)
		{
			int y=x+j;
			if(a==-1)//it's a switch
				p[x][j]=new bool[y-x+2];
			else
				p[x][j]=new bool[min(a+1,y-x+2)];
		}
	}
	
}
allocateRange::~allocateRange()
{	if(n>0)
	release();
}
void allocateRange::release()
{
	//p=new int **[n];
	for(int x=0;x<n;x++)
	{	//x=0,1,...N-1
		
		for(int j=0;j<n-x;j++)
		{
			delete []p[x][j];
			
		}
		delete []p[x];
		//delete []split[x];
		
		//p[x]=new int *[n-x]; //y=x,...N-1
	}
	delete []p;
	n=-1;
	//delete []split;
	
}
bool & allocateRange::operator ()(int x, int y, int i)
{
	return p[x][y-x][i];
	
}
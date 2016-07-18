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
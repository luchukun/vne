#include "wavepath.h"
#include <limits>      // Used for infinity (for integers)
#include <iostream>
// to maintain the wavelengths in the links

Link::Link()
{
	for (int i=0;i<Nw;i++)
	{
		HoldTime[i]=0;
	}
}
Link::Link(const Link &x)
{
	for (int i=0;i<Nw;i++)
		HoldTime[i]=x.HoldTime[i];
}


Clink::Clink(const Clink &xin)
{
	for (int i=0;i<Nw;i++)
		Avai[i]=xin.Avai[i];
}
void Clink::operator =(Link &xin)
{
	for (int i=0;i<Nw;i++)
		if (xin.HoldTime[i]==0)
		Avai[i]=1;		
		else if (xin.HoldTime[i]>0)
		Avai[i]=0;
}
void Clink::operator =(bool c)
{
	for (int i=0;i<Nw;i++)
		Avai[i]=c;
}
void Clink::operator =(const Clink &xin)
{
	for (int i=0;i<Nw;i++)
		Avai[i]=xin.Avai[i];		
		
}
Clink Clink::operator &(const Clink &xin)
{
	Clink y;
	for (int i=0;i<Nw;i++)
		y.Avai[i]=Avai[i]&xin.Avai[i];		
	return y;
		
}

int Clink::MinWave()
{
	
	for (int i=0;i<Nw;i++)
	{
		if (Avai[i]==1)
			{return i; }
	}
	return -1;

}
/*
bool Clink::CommonWave(WCR * WL,bool *S,int & SameWave,int &owave)
{
	int iwave;	
	for (int i=0;i<Nw;i++)
	{	iwave=WL[i].InWave;
		if (WL[i].OutWave)
			if((WL[i].copy)&&(WL[i].InWave!=i)&&(S[iwave]==0))
				{SameWave=iwave;owave=i; return true;}
	}
	return false;

}
bool Clink::CommonWave(WCR * WL,int & SameWave)
{
	
	for (int i=0;i<Nw;i++)
	{
		if (WL[i].OutWave&&WL[i].copy&&WL[i].InWave==i)
		{SameWave=i; return true;}
	}
	return false;

}
*/
bool Clink::IsFull()
{
	for (int i=0;i<Nw;i++)
	{
		if (Avai[i]==1)
			return false;
	}
	return true;

}

// to describe the Edge//

void Edge::Init(const int &d, const int &a,float c)
{
	Dest=d;Addr=a;Cost=c;
}

Edge::Edge(const Edge & E) 
{	Dest=E.Dest;	
	Addr=E.Addr;
	Cost=E.Cost;
}
/*
class WCR
{
public:
	int InWave;
	bool OutWave[Nw];
	bool copy;
	WCR(){}
	Init(const int & i,const bool & o,const bool & c);
};
WCR::Init(const int & i,const bool & o,const bool & c)
{
	InWave=i,OutWave=o;copy=c;
}
*/
void Node::Init(const int &node_name)
{
	Name=node_name;
	MC=true;		
	Degree=0;Dist=Inf;Mark=0;Prev=-1;//OnTree=0;
	//OnWave=0;OnForest=0;
	WaveList=0;
	Wavelength=-1;
	nBranch=0;
}

request::request()
{
	for (int i=0;i<Nv;i++)
	{Dnodes[i]=0;Member[Nv]=0;}
	Type=0;
	Snode=-1;
	MemNumber=0;
	HoldTime=0;
	ConvTimes=0;
}
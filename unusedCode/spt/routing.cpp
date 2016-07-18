#include "routing.h"
#include "frandom.h"
#include "time.h"
#include "wavepath.h"
Graph::Graph()
{	
	int i;
	for (i=0;i<Nv;i++)
	{
		for (i=0;i<Nv;i++)
		{	Table[i].Init(i);
		}
	}	

	/*AllLink=new Link [Max_link];
	if(AllLink == 0 )
	{
		std::cerr << "Insufficient memory for AllLink" << endl;
		//return -1;
	}	
	*/
}
/*Graph::~Graph()
{
	int i;
	for (i=0;i<Nv;i++)
	{	delete [] (Table[i]).Adj;
		delete [] (Table[i]).WaveList;
	}
	delete [] AllLink;
	delete [] Table;

	for (i=0;i<Max_link;i++)
	{	delete [] (AllLink[i]).HoldTime;
	}	
}*/
int Graph::Gen_gragh()
{
	const string file_name("Text1.txt");
	const int MaxLineLength = 256;
    static char OneLine[ MaxLineLength + 1 ];
    int Vs, Vd, N_link=0;;
    float Cost;
	if( file_name !="" )
    {             		
		ifstream Gs(file_name.c_str( ));
		if( Gs )
		{
			while( Gs.getline( OneLine, MaxLineLength ) )
				{
					istrstream Ls( OneLine, MaxLineLength );
					if( Ls >> Vs &&	Ls >> Vd &&Ls >> Cost)
					{	
						//std::cout<<Vs<<"-"<<Vd<<"-"<<Cost<<endl;
						AddEdge( Vs-1, Vd-1, Cost, N_link);
						N_link=N_link+1;
					}
					 else
						 std::cerr << "Bad line: "<<endl ;
				}
			Gs.close();
		}
        else
			std::cerr << "Error opening " << file_name<<endl ;
    } 
	for( int i = 0; i < Nv; i++ )
	{
		for(int j=0;j<Nw;j++)
			Table[i].WaveRange[j]=Table[i].Degree;
	}
	// computed the number of total nodes;
	
	return 1;
		
}

// Add the edge ( Source, Dest, Cost ) to the Graph
void Graph::AddEdge( const int & S1, const int & D1, const float &C1,const int &N_link )
{	
	int i;
	i=Table[S1].Degree;	
	(Table[S1].Adj)[i].Init(D1,N_link,C1);	
	Table[S1].Degree++;	
	i=Table[D1].Degree;
	(Table[D1].Adj)[i].Init(S1,N_link,C1);
	Table[D1].Degree++;
	AllLink[N_link].Node1=S1;
	AllLink[N_link].Node2=D1;

}

void Graph::ClearTable(int type )
{
	if (type==0) //for routing
	{	for( int i = 0; i < Nv; i++ )
		{
			Table[i].Dist = Inf;
			Table[i].Prev = -1;//pointer to its self
			Table[i].Mark = 0;
			Table[i].OnTree=0;	
			Table[i].nBranch=0;
		}
	} 
	
	else if (type==1) // for WAA
	{	
		for( int i = 0; i < Nv; i++ )
		{	Table[i].OnWave=0;
			Table[i].WaveList=0;
			Table[i].Wavelength=-1;			
		}		
	}
	else if(type==2)
	{
		for( int i = 0; i < Nv; i++ )
		{
			for(int j=0;j<Nw;j++)
				Table[i].WaveRange[j]=Table[i].Degree;
		}
	}
//type==1, first time to initate Table	
}

bool Graph::AssignWave(int &Vm,int &Vs,queue<MidNode> &RevNode)
{	
	int V0=Vm;
	int V1,V2,w1,LinkIndex,j;
	int P;
	//std::cout<<"\n"<<"The member "<<Vm+1<<"is accepted on";
	do
	{	if(RevNode.empty()) 
			std::cout<<"\n"<<"RevNode is empty!";
		else
		{
			MidNode x(RevNode.front());
			V1=x.inode;
			w1=x.iwave;
			RevNode.pop();			
			for (j=V0;j!=V1;j=V2)
			{	V2=Table[j].Prev;					
				for (P =0;P<Table[j].Degree&&(Table[j].Adj)[P].Dest!=V2; P=P+1);
				LinkIndex=(Table[j].Adj)[P].Addr;
				float *a=AllLink[LinkIndex].HoldTime;
				a[w1]=Req.HoldTime;
				Table[j].WaveList.Avai[w1]=1;
				Table[j].OnWave=true;
			//	std::cout<<"+ Node "<<j+1<<" Wave "<<w1+1;
				if	(Req.Dnodes[j])	Req.Member[j]=1;
			}
			//Req.ConvTimes=Req.ConvTimes+1;
			
			V0=V1;
		}
	}while(V0!=Vs);

	//Req.ConvTimes=Req.ConvTimes-1;
	while (!RevNode.empty())   {RevNode.pop();}
	return 1;
}

void Graph::WaveBlock(int &Vm,queue<MidNode> &RevNode)
{
	Req.Member[Vm]=false;// fail to reserve waves for node i
	while (!RevNode.empty())   {RevNode.pop();} //clear the reserve results
	//std::cout<<"\n"<<"The member "<<Vm+1<<"is blocked";
}

float Graph::WAA_down2up_MWCn()
{
	int i,Pv,Cv,V,Vm;
	Clink S0,S1,S2;S0=0;
	int P,LinkIndex;
	int mwave, cwave;	
	queue<MidNode> RevNode;//to save the reserve results (middle node and waves)
	Clink E1,E2;
	bool Leaf[Nv];
	for (i=0; i<Nv;i++)	{Leaf[i]=false;Req.Member[i]=0;}
	ClearTable(1);

	for (i=0; i<Nv;i++)
	{
		if ((Req.Dnodes[i]==1)&&(!Req.Member[i]))
		{	
			S0=0;
			bool RevEnd=false;// end of the reserve process
			V=i;Vm=i;Cv=i;		
			while(!RevEnd)	// to assign wave for a member
			{	
				if (V==Req.Snode||Table[Cv].OnWave==1) //source or leaf node?
				{	
					if (S0.IsFull())
					{	WaveBlock(i,RevNode);
						RevEnd=true;	
					}
					else
					{	mwave=S0.MinWave();		//busy waves in c(v)-v
						MidNode nd1(V,mwave);
						RevNode.push(nd1);// put a intermidate node
						RevEnd=true;
						Req.Member[i]=true;
						Leaf[i]=true;
						if (Leaf[V]==true) Leaf[V]=false; // to connect to a original leaf node
						AssignWave(i,V,RevNode); //assign waves from mem[i] to V
					}
				}//if (V==Req.Snode||Leaf[V]==true)
				else if(Leaf[V]==true)
				{	
					if (S0.IsFull())
					{	WaveBlock(i,RevNode);
						RevEnd=true;	
					}
					else 
					{	
						cwave=(Table[V].WaveList&S0).MinWave();
						if(cwave>=0&&cwave<Nw)//cwave is the input wave of node v
						{	mwave=cwave;
							MidNode nd1(V,mwave);
							RevNode.push(nd1);
						}
						else  
						{
							Req.ConvTimes++;
							mwave=S0.MinWave();		//busy waves in c(v)-v
							MidNode nd1(V,mwave);
							RevNode.push(nd1);// put a intermidate node
						}
						RevEnd=true;
						Req.Member[i]=true;
						Leaf[i]=true;
						if (Leaf[V]==true) Leaf[V]=false; // to connect to a original leaf node
						AssignWave(i,V,RevNode);						
						//Table[V].WaveList[mwave].Init(cwave,true,false);
					
					}					
				}
				else
				{	
					Pv=Table[V].Prev;					
					for (P =0;P<Table[V].Degree&&(Table[V].Adj)[P].Dest!=Pv; P=P+1);
					LinkIndex=(Table[V].Adj)[P].Addr;
					S1=AllLink[LinkIndex];//S1= busy waves in v-P(V)
					S2=S0&S1;				
				    if (S1.IsFull()==true)
					{	
						WaveBlock(i,RevNode);
						RevEnd=true;						
					}
					else
					{											
						if (S2.IsFull()==true)
						{	mwave=S0.MinWave();
							MidNode nd1(V,mwave);
							RevNode.push(nd1);					
							S0=S1;
							//Vm=V;						
						}
						else	S0=S2;						
					}																	
				} //if (V==Req.Snode||Leaf[V]==true)
				Cv=V;
				V=Pv;	
			}//while(!RevEnd)
		}//if (Req.Dnodes[i]==1)
	}//for (i=0; i<Nv;i++)	

	int Mem_included=0;float block_rate=0;
	for (i=0;i<Nv;i++)
	{
		if (Req.Member[i])
		Mem_included++;				
	}
	
	//block_rate=(1-float(Mem_included)/float(Req.MemNumber));
	block_rate=float(Req.MemNumber-Mem_included);
	return block_rate;

}

float Graph::WAA_down2up_Sparse_MWCn()
{
	int i,Pv,Cv,V,Vm;
	Clink S0,S1,S2;
	S0=1;
	int P,LinkIndex;
	int mwave, cwave;	
	queue<MidNode> RevNode;//to save the reserve results (middle node and waves)
	Clink E1,E2;
	bool Leaf[Nv];
	for (i=0; i<Nv;i++)	{Leaf[i]=false;Req.Member[i]=0;}
	ClearTable(1);
	for (i=0; i<Nv;i++)
	{
		if ((Req.Dnodes[i]==1)&&(!Req.Member[i]))
		{	
			S0=1;
			bool RevEnd=false;// end of the reserve process
			V=i;Vm=i;Cv=i;Pv=i;			
			while(!RevEnd)	// to assign wave for a member
			{	
				if (V==Req.Snode||Table[Cv].OnWave||(Table[V].OnWave&&Table[V].MC)) //source or leaf node?
				{	
					if (S0.IsFull())
					{	WaveBlock(i,RevNode);
						RevEnd=true;	
					}
					else
					{	mwave=S0.MinWave();		//busy waves in c(v)-v
						MidNode nd1(V,mwave);
						RevNode.push(nd1);// put a intermidate node
						RevEnd=true;
						Req.Member[i]=true;
						Leaf[i]=true;
						if (Leaf[V]==true) Leaf[V]=false; // to connect to a original leaf node
						AssignWave(i,V,RevNode); //assign waves from mem[i] to V
					}
				}//if (V==Req.Snode||Leaf[V]==true)
				else if(Leaf[V]==true)
				{	
					if (S0.IsFull())
					{	WaveBlock(i,RevNode);
						RevEnd=true;	
					}
					else 
					{	
						cwave=(Table[V].WaveList&S0).MinWave();
						if(cwave>=0&&cwave<Nw)//cwave is the input wave of node v
						{	mwave=cwave;
							MidNode nd1(V,mwave);
							RevNode.push(nd1);
						}
						else  
						{
							Req.ConvTimes++;
							mwave=S0.MinWave();		//busy waves in c(v)-v
							MidNode nd1(V,mwave);
							RevNode.push(nd1);// put a intermidate node
						}
						RevEnd=true;
						Req.Member[i]=true;
						Leaf[i]=true;
						if (Leaf[V]==true) Leaf[V]=false; // to connect to a original leaf node
						AssignWave(i,V,RevNode);						
						//Table[V].WaveList[mwave].Init(cwave,true,false);
					
					}					
				}
				else
				{	
					Pv=Table[V].Prev;					
					for (P =0;P<Table[V].Degree&&(Table[V].Adj)[P].Dest!=Pv; P=P+1);
					LinkIndex=(Table[V].Adj)[P].Addr;
					S1=AllLink[LinkIndex];//S1= busy waves in v-P(V)
					S2=S0&S1;				
				    if (S1.IsFull()==true)
					{	
						WaveBlock(i,RevNode);
						RevEnd=true;						
					}
					else
					{											
						if (S2.IsFull()==true)
						{	mwave=S0.MinWave();
							MidNode nd1(V,mwave);
							RevNode.push(nd1);					
							S0=S1;
							//Vm=V;						
						}
						else	S0=S2;						
					}																	
				} //if (V==Req.Snode||Leaf[V]==true)
				Cv=V;
				V=Pv;	
			}//while(!RevEnd)
		}//if (Req.Dnodes[i]==1)
	}//for (i=0; i<Nv;i++)	

	int Mem_included=0;float block_rate=0;
	for (i=0;i<Nv;i++)
	{
		if (Req.Member[i])
		Mem_included++;				
	}
	block_rate=float(Req.MemNumber-Mem_included);
	//block_rate=(1-float(Mem_included)/float(Req.MemNumber));
	return block_rate;
}

float Graph::WAA_down2up()
{
	int i,Pv,Cv,V,Vm;
	Clink S0,S1,S2;
	S0=1;
	int P,LinkIndex;
	int mwave, cwave;	
	queue<MidNode> RevNode;//to save the reserve results (middle node and waves)
	Clink E1,E2;
	bool Leaf[Nv];
	for (i=0; i<Nv;i++)	{Leaf[i]=false;Req.Member[i]=0;}
	ClearTable(1);
	for (i=0; i<Nv;i++)
	{
		if ((Req.Dnodes[i]==1)&&(!Req.Member[i]))
		{	
			S0=1;
			bool RevEnd=false;// end of the reserve process
			V=i;Vm=i;Cv=i;Pv=i;			
			while(!RevEnd)	// to assign wave for a member
			{	
				if (V==Req.Snode) //source or leaf node?
				{	
					if (S0.IsFull())
					{	WaveBlock(i,RevNode);
						RevEnd=true;	
					}
					else
					{	mwave=S0.MinWave();		//busy waves in c(v)-v
						MidNode nd1(V,mwave);
						RevNode.push(nd1);// put a intermidate node
						RevEnd=true;
						Req.Member[i]=true;
						Leaf[i]=true;
						if (Leaf[V]==true) Leaf[V]=false; // to connect to a original leaf node
						AssignWave(i,V,RevNode); //assign waves from mem[i] to V
					}
				}//if (V==Req.Snode||Leaf[V]==true)
				else if(Leaf[V]==true||Table[V].OnWave)
				{	
					if (S0.IsFull())
					{	WaveBlock(i,RevNode);
						RevEnd=true;	
					}
					else 
					{	
						cwave=(Table[V].WaveList&S0).MinWave();
						if(cwave>=0&&cwave<Nw)//cwave is the input wave of node v
						{	mwave=cwave;
							MidNode nd1(V,mwave);
							RevNode.push(nd1);
						}
						else  
						{
							Req.ConvTimes++;
							mwave=S0.MinWave();		//busy waves in c(v)-v
							MidNode nd1(V,mwave);
							RevNode.push(nd1);// put a intermidate node
						}
						RevEnd=true;
						Req.Member[i]=true;
						Leaf[i]=true;
						if (Leaf[V]==true) Leaf[V]=false; // to connect to a original leaf node
						AssignWave(i,V,RevNode);						
						//Table[V].WaveList[mwave].Init(cwave,true,false);
					
					}					
				}
				else
				{	
					Pv=Table[V].Prev;					
					for (P =0;P<Table[V].Degree&&(Table[V].Adj)[P].Dest!=Pv; P=P+1);
					LinkIndex=(Table[V].Adj)[P].Addr;
					S1=AllLink[LinkIndex];//S1= busy waves in v-P(V)
					S2=S0&S1;				
				    if (S1.IsFull()==true)
					{	
						WaveBlock(i,RevNode);
						RevEnd=true;						
					}
					else
					{											
						if (S2.IsFull()==true)
						{	mwave=S0.MinWave();
							MidNode nd1(V,mwave);
							RevNode.push(nd1);					
							S0=S1;
							//Vm=V;						
						}
						else	S0=S2;						
					}																	
				} //if (V==Req.Snode||Leaf[V]==true)
				Cv=V;
				V=Pv;	
			}//while(!RevEnd)
		}//if (Req.Dnodes[i]==1)
	}//for (i=0; i<Nv;i++)	

	int Mem_included=0;float block_rate=0;
	for (i=0;i<Nv;i++)
	{
		if (Req.Member[i])
		Mem_included++;				
	}
	
	//block_rate=(1-float(Mem_included)/float(Req.MemNumber));
	block_rate=float(Req.MemNumber-Mem_included);
	return block_rate;
}
void Graph::anyWAA(int Snode)
{	int Pv,V;
	Clink S1;S1=1;
	int P,LinkIndex;
	V=Snode;
	Pv=Table[V].Prev;					
	int Pnext=-1;
	bool RevEnd=0;
	while(!RevEnd)
	{	if(Table[V].nBranch==0)	
			RevEnd=1;
		else if(Table[V].nBranch==1)
		{	
			for (P =0;P<Table[V].Degree; P=P+1)
			{	Pnext=(Table[V].Adj)[P].Dest;		
				if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
				{	//nBranch=nBranch+1;		
					LinkIndex=(Table[V].Adj)[P].Addr;					
					S1=AllLink[LinkIndex];
					if(S1.IsFull())
					{	RevEnd=1;
					}
					else
					{	Table[Pnext].Wavelength=S1.MinWave();						
						
					}
					P=Table[V].Degree;
				}	
				
			}
			Pv=V;
			V=Pnext;
		}
		else if(Table[V].nBranch>1)
		{	
			for (P =0;P<Table[V].Degree; P=P+1)
			{	Pnext=(Table[V].Adj)[P].Dest;		
				if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
				{	LinkIndex=(Table[V].Adj)[P].Addr;
					S1=AllLink[LinkIndex];
					if(!S1.IsFull())
					{									
						Table[Pnext].Wavelength=S1.MinWave();
						anyWAA(Pnext);					
					}
				}		
			}
			RevEnd=true;			
		}		
	}
}
void Graph::oneWAA(int Snode)
{	int Pv,V;
	Clink S1;S1=1;
	int P,LinkIndex;
	int mwave;	
	V=Snode;
	Pv=Table[V].Prev;					
	int Pnext=-1;
	int wtmp[Nw];
	bool RevEnd=0;
	while(!RevEnd)
	{	if(Table[V].nBranch==0)	
			RevEnd=1;
		else if(Table[V].nBranch==1)
		{	
			for (P =0;P<Table[V].Degree; P=P+1)
			{	Pnext=(Table[V].Adj)[P].Dest;		
				if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
				{	//nBranch=nBranch+1;
					
					LinkIndex=(Table[V].Adj)[P].Addr;					
					S1=AllLink[LinkIndex];
					if(S1.IsFull())
					{	RevEnd=1;
					}
					else
					{	Table[Pnext].Wavelength=S1.MinWave();						
						
					}
					P=Table[V].Degree;
				}	
				
			}
		Pv=V;
		V=Pnext;
		}
		else if(Table[V].nBranch>1)
		{	for( int i=0;i<Nw;i++)	wtmp[i]=0;			
			for (P =0;P<Table[V].Degree; P=P+1)
			{	Pnext=(Table[V].Adj)[P].Dest;		
				if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
				{	LinkIndex=(Table[V].Adj)[P].Addr;
					S1=AllLink[LinkIndex];
					for( int i=0;i<Nw;i++)	{if (S1.Avai[i]==1)	wtmp[i]++;}		
				}		
			}
			mwave=-1;int max=0;
			for( int i=0;i<Nw;i++)
			{	if(wtmp[i]==Table[V].nBranch)	
				{max=Table[V].nBranch; mwave=i; i=Nw;}
				else if(wtmp[i]>max)
				{max=wtmp[i];mwave=i;}
			}
			if(max==0)
			{	RevEnd=1;}
			else
			{	for (P =0;P<Table[V].Degree; P=P+1)
				{	Pnext=(Table[V].Adj)[P].Dest;		
					if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
					{	LinkIndex=(Table[V].Adj)[P].Addr;
						S1=AllLink[LinkIndex];
						if (S1.Avai[mwave]==1)			
						{	Table[Pnext].Wavelength=mwave;
							oneWAA(Pnext);
						}
					}		
				}
			}
			RevEnd=1;
		}		
	}
}
void Graph::oneWAA_v1(int Snode)
{	int Pv,V;
	Clink S1;S1=1;
	int P,LinkIndex;
	int mwave,tmp;	
	int maxa=0,minb=Nd;
	V=Snode;
	//cout<<"\n"<<"Snode="<<Snode+1;
	Pv=Table[V].Prev;					
	int Pnext=-1;
	int wtmp[Nw];
	bool RevEnd=0;
	while(!RevEnd)
	{	if(Table[V].nBranch==0)	
		{	RevEnd=1;
			//cout<<"Node "<<V+1<<" is a leaf node";
			//if(Table[V].Wavelength==-1) 
			//cout<<"\n"<<"Node "<<V+1<<" is blocked";
		}
		else if(Table[V].nBranch==1)
		{	//cout<<"\n"<<"Node "<<V+1<<" is a Mid node";
			for (P =0;P<Table[V].Degree; P=P+1)
			{	Pnext=(Table[V].Adj)[P].Dest;		
				if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
				{	//nBranch=nBranch+1;
					LinkIndex=(Table[V].Adj)[P].Addr;					
					S1=AllLink[LinkIndex];
					if(S1.IsFull())
					{	RevEnd=1;
					}
					else
					{	//Table[Pnext].Wavelength=S1.MinWave();	
						mwave=-1;minb=Table[V].Degree+1;
						for (int i=0;i<Nw;i++)
						{	tmp=Table[V].WaveRange[i];
							if (S1.Avai[i]==1&&tmp<minb)
							{mwave=i;minb=tmp; }
						}
						if(mwave==-1)
							RevEnd=1;
						else
						{	//cout<<"\n"<<"Node "<<Pnext+1<<" is accepted on Wave:"<<mwave+1;
							Table[Pnext].Wavelength=mwave;						
						}
					}
					P=Table[V].Degree;
				}
				
			}
			Pv=V;
			V=Pnext;
		}
		else if(Table[V].nBranch>1)
		{	//cout<<"\n"<<"Node "<<V+1<<" is a MC node :";
			for( int i=0;i<Nw;i++)	wtmp[i]=0;			
			for (P =0;P<Table[V].Degree; P=P+1)
			{	Pnext=(Table[V].Adj)[P].Dest;		
			if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
				{	LinkIndex=(Table[V].Adj)[P].Addr;
					S1=AllLink[LinkIndex];
					for( int i=0;i<Nw;i++)	{if (S1.Avai[i]==1)	wtmp[i]++;}		
				}		
			}
			mwave=-1;maxa=0;minb=Table[V].Degree+1;
			for( int i=0;i<Nw;i++)
			{	tmp=Table[V].WaveRange[i];			
				if(wtmp[i]>Table[V].nBranch)
				cout<<"error in wtmp of oneWAA_v1";
				else if(wtmp[i]>maxa)
				{	maxa=wtmp[i];minb=tmp;mwave=i;}
				else if((wtmp[i]==maxa)&&(tmp<minb))
				{	maxa=wtmp[i];minb=tmp;mwave=i;}
				
			}
			if(mwave==-1)
			{	RevEnd=1;
			}
			else
			{	for (P =0;P<Table[V].Degree; P=P+1)
				{	Pnext=(Table[V].Adj)[P].Dest;		
					if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
					{	LinkIndex=(Table[V].Adj)[P].Addr;
						S1=AllLink[LinkIndex];
						if (S1.Avai[mwave]==1)			
						{	//cout<<"\n"<<"Node "<<Pnext+1<<" is accepted on Wave:"<<mwave+1;
							Table[Pnext].Wavelength=mwave;
						if(Table[Pnext].nBranch>0)
							oneWAA_v1(Pnext);
						}
					}		
				}
				RevEnd=1;
			}
			
		}		
	}
}
float Graph::WAA_up2down(sim_type WAA_flag)
{
	int Pv,V;
	int P,LinkIndex;
	int mwave;	
	ClearTable(1);
	for (int i=0; i<Nv;i++)	Req.Member[i]=0;
	if (WAA_flag==fWC_up2down)	anyWAA(Req.Snode);
	if(WAA_flag==sWC_up2down)	oneWAA(Req.Snode);
	if(WAA_flag==sWC_up2down_v1)oneWAA_v1(Req.Snode);
	for (int i=0; i<Nv;i++)
	{	if(Req.Dnodes[i]&&Table[i].OnWave==false)
		{	//cout<<"\n";
			bool RevEnd=false;// end of the reserve process
			V=i;		
			while(!RevEnd)	// to assign wave for a member
			{	if(V==Req.Snode==true||Table[V].OnWave)	RevEnd=true;							
				else
				{	Pv=Table[V].Prev;					
					for (P =0;P<Table[V].Degree&&(Table[V].Adj)[P].Dest!=Pv; P=P+1);
					LinkIndex=(Table[V].Adj)[P].Addr;
					mwave=Table[V].Wavelength;
					if (mwave==-1)
					{	RevEnd=true;
						//cout<<"Node "<<i+1<<" is blocked"<<"\n";
					}
					else
					{	if(Req.Dnodes[V])	Req.Member[V]=true;
						float *a=AllLink[LinkIndex].HoldTime;
						a[mwave]=Req.HoldTime;
						Table[V].WaveRange[mwave]=Table[V].WaveRange[mwave]-1;
						Table[Pv].WaveRange[mwave]=Table[Pv].WaveRange[mwave]-1;
						Table[V].OnWave=true;
						//cout<<"+ Node "<<V+1<<" Wave "<<mwave+1;
						V=Pv;
						
					}			    
				} 					
			}
			
		}
		
	}	
	//cout<<"\n";
	int Mem_included=0;float block_rate=0;
	for (int i=0;i<Nv;i++)
	{
		if (Req.Member[i])
		Mem_included++;				
	}
	block_rate=float(Req.MemNumber-Mem_included);
	//block_rate=(1-float(Mem_included)/float(Req.MemNumber));
	return block_rate;
}

int Graph::Dijkstra()
{
	int V,U, W;
    BinaryHeap<Path> PQ(Path(Req.Snode, 0));
    Path Lnext;    // Stores the result of a DeleteMin
	int P;
    ClearTable(0);
    Table[Req.Snode].Dist = 0;
	
    PQ.Insert(Path(Req.Snode, 0));
	int Mem_included=0;
//	int Leaf[Nv];
	Table[Req.Snode].OnTree=1;
    while (Mem_included<Req.MemNumber)
    {
        do
        {	if( PQ.IsEmpty( ))
			{	cout<<"PQ is Empty";
				return 1;
			}
            PQ.DeleteMin(Lnext);
        } while(Table[Lnext.Vdest].Mark);
        V = Lnext.Vdest;
		//std::cout<<" "<<V+1<<" ";
        Table[V].Mark = 1;     // Mark vertex as being seen
		for (P =0;P<Table[V].Degree; P=P+1)
        {
            W = (Table[V].Adj)[P].Dest;
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
                PQ.Insert(Path( W, Table[W].Dist ) );
            }
        }
		if (Req.Dnodes[V] == 1) 
		{	Mem_included=Mem_included+1;
			U=V;
			while(U!=Req.Snode&&Table[U].OnTree!=1)	
			{	Table[U].OnTree=1;
				U=Table[U].Prev;
			}

			//std::cout<<"\n"<<V+1<<" is routed!\n ";
		}
    }
	int Pnext;
	for (V=0;V<Nv;V++)
	{
		for (int P =0;P<Table[V].Degree; P=P+1)
		{	
			Pnext=(Table[V].Adj)[P].Dest;		
			if(Table[Pnext].Prev==V&&Table[Pnext].OnTree==1)
			{	
				Table[V].nBranch=Table[V].nBranch+1;		
			}		
		}
		//cout<<"\n"<<V+1<<" 's Branch="<<Table[V].nBranch;
	}
	//cout<<"\n";
	return 1;
}
/*
int Graph::Member_only( )//muticast routing algorithm
{
	int V,U, W;
    BinaryHeap<Path> PQ(Path(Req.Snode, 0));
    Path Pnext;    // Stores the result of a DeleteMin
	int P;
    ClearTable(0);
    Table[Req.Snode].Dist = 0;
    PQ.Insert(Path(Req.Snode, 0));
	int Tree_number=0;
	int Mem_included=0;
    while (Mem_included<Req.MemNumber)
    {
        do
        {
            if( PQ.IsEmpty( ) )		
			{		
				Tree_number=Tree_number+1;
				ClearTable(1);
				Table[Req.Snode].Dist = 0;
				PQ.Insert(Path(Req.Snode, 0));
				// return 1;
			}
			else
               
            PQ.DeleteMin(Pnext);
        } while(Table[Pnext.Vdest].Mark);

        V = Pnext.Vdest;
        Table[V].Mark = 1;     // Mark vertex as being seen
		if (Req.Dnodes[V]==0)	// V is not a Member
			for (P =0;P<Table[V].Degree; P=P+1)
			{
				W = (Table[V].Adj)[P].Dest;
				int link_addr=(Table[V].Adj)[P].Addr;
				float Cvw = AllLink[link_addr].Cost;
				if( Cvw < 0 )
				{	std::cerr << "Graph has negative edges" << endl; 
					return 0;
				}
				// update the dist of v's neighbour v
				if( Table[W].Dist > Table[V].Dist + Cvw )
				{
					Table[W].Dist = Table[V].Dist + Cvw;
					Table[W].Prev = V;
					PQ.Insert( Path( W, Table[W].Dist ) );
				}
			}
		else //V is a Member
		{	Mem_included=Mem_included+1;
			U=V;
			PQ.MakeEmpty( );
			do
			{	Table[U].OnTree=1;
				Table[U].Onforest=1;
				Table[U].Degree=Table[U].Degree-1;//u->v then u--;
				U=Table[U].Prev;
				Table[U].Degree=Table[U].Degree-1;//U<-u; then v--
			}while(U!=Req.Snode);

			for( int i = 0; i < Nv; i++ )
			{
				if (Table[i].OnTree==0)
				{	Table[i].Dist = Inf;
					Table[i].Prev = -1;
					Table[i].Mark = 0;
				}
				else 
					if (Table[i].Degree>0)
						for (P =0;P<Table[V].Degree; P=P+1)
						{
							W = (Table[V].Adj)[P].Dest;							
							int link_addr=(Table[V].Adj)[P].Addr;
							float Cvw = AllLink[link_addr].Cost;
							if (Table[W].OnTree==0)
							{	if( Cvw < 0 )
							{	std::cerr << "Graph has negative edges" << endl;
									return 0;
								}
								if( Table[W].OnTree)
									PQ.Insert(Path(W, Cvw));									
							}
						}				
			}

		}
    }
	return 1;
}

int Graph::CDrouting( )//muticast routing algorithm
{
	int V,U, W;
    BinaryHeap<Path> PQ(Path(Req.Snode, 0));
    Path Pnext;    // Stores the result of a DeleteMin
	int P;
    ClearTable(0);
    Table[Req.Snode].Dist = 0;
    PQ.Insert(Path(Req.Snode, 0));
	int Tree_number=0;
	int Mem_included=0;
    while (Mem_included<Req.MemNumber)
    {
        do
        {
            if( PQ.IsEmpty( ) )		
			{		
				Tree_number=Tree_number+1;
				ClearTable(1);
				Table[Req.Snode].Dist = 0;
				PQ.Insert(Path(Req.Snode, 0));
			}
               
            PQ.DeleteMin(Pnext);
        } while(Table[Pnext.Vdest].Mark);

        V = Pnext.Vdest;
        Table[V].Mark = 1;     // Mark vertex as being seen
		if (Req.Dnodes[V]==0)	// V is not a Member
			for (P =0;P<Table[V].Degree; P=P+1)
			{
				W = (Table[V].Adj)[P].Dest;
				int link_addr=(Table[V].Adj)[P].Addr;
				float Cvw = AllLink[link_addr].Cost;
				if( Cvw < 0 )
				{	std::cerr << "Graph has negative edges" << endl; 
					return 0;
				}
				// update the dist of v's neighbour v
				if( Table[W].Dist > Table[V].Dist + Cvw )
				{
					Table[W].Dist = Table[V].Dist + Cvw;
					Table[W].Prev = V;
					PQ.Insert( Path( W, Table[W].Dist ) );
				}
			}
		else //V is not a Member
		{	Mem_included=Mem_included+1;
			U=V;
			PQ.MakeEmpty( );
			do
			{	Table[U].OnTree=1;
				Table[U].Onforest=1;
				Table[U].Degree=Table[U].Degree-1;//u->v then u--;
				U=Table[U].Prev;
				Table[U].Degree=Table[U].Degree-1;//U<-u; then v--
			}while(U!=Req.Snode);

			for( int i = 0; i < Nv; i++ )
			{
				if (Table[i].OnTree==0)
				{	Table[i].Dist = Inf;
					Table[i].Prev = -1;
					Table[i].Mark = 0;
				}
				else 
					if (Table[i].Degree>0)
						for (P =0;P<Table[V].Degree; P=P+1)
						{
							W = (Table[V].Adj)[P].Dest;							
							int link_addr=(Table[V].Adj)[P].Addr;
							float Cvw = AllLink[link_addr].Cost;
							if (Table[W].OnTree==0)
							{	if( Cvw < 0 )
							{	std::cerr << "Graph has negative edges" << endl;
									return 0;
								}
								if( Table[W].OnTree)
									PQ.Insert(Path(W, Cvw));									
							}
						}				
			}

		}
    }
	return 1;
}

*/
void Graph::PrintPathRec( int DestNode ) const
{
    if( Table[DestNode].Prev !=-1 )
    {
        PrintPathRec( Table[DestNode].Prev );
		std::cout << " to ";
		//bool *a=Table[DestNode].WaveAssigned;
		
		//for (int i=0;i<Nw;i++)
		//{	if (a[i]) 
			//std::cout<<" Wave["<<i+1<<"] ";
		//}	
    }
    std::cout << Table[DestNode].Name+1;
}

// Driver routine to handle unreachables and print total cost
// It calls recursive routine to print shortest path to
// DestNode after a hortest path algorithm has run

void Graph::PrintPath() const
{	for (int i=0;i<Nv;i++)
		if (Req.Dnodes[i])
		{	if( Table[i].Dist ==Inf)
				std::cout << Table[i].Name+1 << " is unreachable";
			else
			{	std::cout <<"The route to Node "<<i+1<<" is: ";
				PrintPathRec(i);
				std::cout << " (cost is " << Table[i].Dist << ")";
			}
			std::cout << endl;
		}
}

void Graph::LinkUpdate(float & arrive_time)
{	float * remain_time;
	int Nd1,Nd2;
	for (int i=0;i<Ne;i++)
	{	remain_time=AllLink[i].HoldTime;
		//std::cout<<"\n"<<"link "<<i+1<<" 's free waves are:";
		for (int j=0;j<Nw;j++)
		{	
			if (remain_time[j]>arrive_time)
				remain_time[j]=remain_time[j]-arrive_time;
			else
			{	remain_time[j]=0;
				Nd1=AllLink[i].Node1;
				Nd2=AllLink[i].Node2;
				Table[Nd1].WaveRange[j]=Table[Nd1].WaveRange[j]+1;
				Table[Nd2].WaveRange[j]=Table[Nd2].WaveRange[j]+1;
			}


		//	if	(remain_time[j]==0)
			//	std::cout<<j+1<<"-";
		}
	}
}
void Graph::LinkRandom(const float & link_ratio)
{	float * remain_time;
	for (int i=0;i<Ne;i++)
	{	remain_time=AllLink[i].HoldTime;
		std::cout<<"\n"<<"link "<<i+1<<" 's free waves are:";
		for (int j=0;j<Nw;j++)
		{				
			remain_time[j]=(float)100*rand_b01(link_ratio);
			if	(remain_time[j]==0)
				std::cout<<j+1<<"-";
		}
	}

}
void Graph::LinkClear()
{	
	for (int i=0;i<Ne;i++)
	{	
		for (int j=0;j<Nw;j++)
		{	AllLink[i].HoldTime[j]=0;
		}
	}
}
void Graph::fMC()
{
	int i;
	for (i=0;i<Nv;i++)
	{						
		Table[i].MC=1;			
	}
}
void Graph::MCRandom(const float & MC_ratio)
{
	int i;
	for (i=0;i<Nv;i++)
	{						
		Table[i].MC=rand_b01(MC_ratio);			
	}
}

void Graph::MCMax(int *NodePriority)
{
	int s,i,j,P,W;
	int *route_degree=new int [Nv],temp_degree=0;
	//int *NodePriority=new int [Nv];
	for (i=0;i<Nv;i++)
	{	
		route_degree[i]=0;
	}
	for (s=0;s<Nv;s++)
	{	
		Req.Snode=s;
		for (i=0;i<Nv;i++)
		{	if (i!=s)
				Req.Dnodes[i]=true;
			else 
				Req.Dnodes[i]=false;
		}		
		Req.MemNumber=Nv-1;
		Dijkstra();	
		for (i=0;i<Nv;i++)
		{	temp_degree=0;
			for (P =0;P<Table[i].Degree; P=P+1)
			{
				W = (Table[i].Adj)[P].Dest;
				if (Table[W].Prev==i)
					temp_degree++;
					
			}
			if (temp_degree>=2)
			route_degree[i]=route_degree[i]+(temp_degree-1);
		}		
	}
	int max_degree=0;
	int max_deg_node=-1;
	//for (i=0;i<Nv;i++)
		//std::cout<<"\n"<<"Node["<<i<<"]="<<route_degree[i];

	for (i=0;i<Nv;i++)
	{	max_degree=0;
		for (j=0;j<Nv;j++)
		{
			if (route_degree[j]>max_degree)
			{
				max_degree=route_degree[j];
				max_deg_node=j;
			}
		}	
		NodePriority[i]=max_deg_node;
		route_degree[max_deg_node]=-1;
	}
	
}




float Graph::ProcessRequest(int SimTime,float load,float mul_ratio,sim_type WAA_flag)
{
    float arri_time=0;
	int i,j,k,t;
//	float load=0.5,mem_ratio=0.5;
	init_genrand( (unsigned long)time( NULL ));
	int Mem_included=0;
	float block_rate=0;
	int SimT=SimTime;
	LinkClear();
	//cout<<"\n"<<"load="<<load<<" mul_ratio="<<mul_ratio;
	for (t=1;t<SimTime;t++)
	{	/*	define the arrivals by the mem_ratio
		LinkUpdate(arri_time);
		//Req.HoldTime=0.5;
		Req.HoldTime=exprnd(1);		
		Req.Snode=unif_int(Nv)-1;
		Req.MemNumber=0;
		for (i=0;i<Nv;i++)
		{			
			if 	(i==Req.Snode)	Req.Dnodes[i]=0;
			else	Req.Dnodes[i]=rand_b01(mem_ratio);
			Req.Member[i]=0;
			if (Req.Dnodes[i])	Req.MemNumber++;
		}
		*/
		LinkUpdate(arri_time);
		Req.HoldTime=exprnd(1);	
		Req.Type=rand_b01(mul_ratio);
		Req.Snode=unif_int(Nv)-1;
		if(Req.Type)
			Req.MemNumber=unif_int(Nv-2)+1;//[2,Nv-1]
		else
			Req.MemNumber=1;	
		
		for (i=0;i<Nv;i++)
		{	Req.Dnodes[i]=0;
			Req.Member[i]=0;
		}
		Req.Dnodes[Req.Snode]=1;
		for (i=1;i<=Req.MemNumber;i=i+1)
		{	int temp=0;
			temp=unif_int(Nv-i);
			for(j=0,k=0;j<temp&&k<Nv;)
			{	if (Req.Dnodes[k]==0)	j=j+1;
					k=k+1;
			}
			Req.Dnodes[k-1]=1;
		}
		Req.Dnodes[Req.Snode]=0;
		/*
			cout<<"\n"<<"Req:"<<t<<" HoldTime is: "<<Req.HoldTime
		<<" ArriveTime is: "<<arri_time;
		cout<<" MemNumber="<<Req.MemNumber<<" ";
		cout<<"\n"<<"Source is: "<<Req.Snode+1<<"; ";
		cout<<" Dest is: ";
		for(int i=0;i<Nv;i++)
			if(Req.Dnodes[i]) cout<<" "<<i+1;		
		std::cout<<"\n";
		*/
		Dijkstra();	
		//PrintPath();
		switch (WAA_flag)
		{	case MWC_down2up:
				block_rate+=WAA_down2up_MWCn();break;
			case fWC_down2up:
				block_rate+=WAA_down2up();break;
			case fWC_up2down:
				block_rate+=WAA_up2down(WAA_flag);break;
			case sWC_up2down:
				block_rate+=WAA_up2down(WAA_flag);break;
			case sWC_up2down_v1:
				ClearTable(2);
				block_rate+=WAA_up2down(WAA_flag);break;
			case MWC_sparse:
				block_rate+=WAA_down2up_Sparse_MWCn();break;
			default:break;
		}
		arri_time=exprnd(load);
		//arri_time=1;
	}

	block_rate=block_rate/float(SimT);
	//std::cout<<"\n"<<"Blocking probabiltiy is: "<<block_rate;
	//std::cout<<block_rate<<",";
    return block_rate;
}

int Graph::ProcessRequest(int SimTime,float load,float mul_ratio,float &unicast_bp,float &multicast_bp)
{
    float arri_time=0;
	int i,j,k,t;
//	float load=0.5,mem_ratio=0.5;
	init_genrand( (unsigned long)time( NULL ));
	//init_genrand(0);
	int Mem_included=0;	
		//cout<<"\n"<<"n="<<n<<"Multicast ratio="<<mul_ratio;
	int unicast_t=0,multicast_t=0;
	multicast_bp=0;unicast_bp=0;	
	for (t=1;t<SimTime;t++)
	{	
			//cout<<"\n"<<"t="<<t<<" ";	
		LinkUpdate(arri_time);
			//Req.HoldTime=0.5;
		Req.HoldTime=exprnd(1);	
		Req.Type=rand_b01(mul_ratio);
		Req.Snode=unif_int(Nv)-1;
			//cout<<"Snode="<<Req.Snode<<" ";
		if(Req.Type)
			Req.MemNumber=unif_int(Nv-1)+1;//[2,Nv]
		else
			Req.MemNumber=1;	
			//cout<<"MemNumber="<<Req.MemNumber<<" ";
			
		for (i=0;i<Nv;i++)
		{	Req.Dnodes[i]=0;
				Req.Member[i]=0;
		}
		Req.Dnodes[Req.Snode]=1;
			//cout<<"Dnodes are: ";
		for (i=1;i<=Req.MemNumber;i=i+1)
		{	int temp=0;
			temp=unif_int(Nv-i);
			for(j=0,k=0;j<temp&&k<Nv;)
			{	if (Req.Dnodes[k]==0)	j=j+1;
					k=k+1;
			}
			Req.Dnodes[k-1]=1;
		}
			
		Req.Dnodes[Req.Snode]=0;
		Dijkstra();	
			//PrintPath();	
		WAA_down2up(); 
				//WAA_down2up_MWCn();
		Mem_included=0;
		for (i=0;i<Nv;i++)
		{
			if (Req.Member[i])
				Mem_included++;				
		}	
		if (Req.Type)
		{	multicast_bp=multicast_bp+(1-(float) Mem_included/(float)(Req.MemNumber));
				multicast_t++;
		}
		else
		{	unicast_bp=unicast_bp+(1-(float)Mem_included/(float)Req.MemNumber);
			unicast_t++;
		}	
		arri_time=exprnd(load);
	}
	if (multicast_t>0)
		multicast_bp=multicast_bp/float(multicast_t);
	else 
		cerr<<"number of multicast is 0 ";
	if (unicast_t>0)
		unicast_bp=unicast_bp/float(unicast_t);
	else 
		cerr<<"number of unicast is 0 ";
			
    return 1;
}

int Graph::UnicastRequest(int SimTime,float load,float &unicast_bp)
{
    float arri_time=0;
	int i,j,k,t;
//	float load=0.5,mem_ratio=0.5;
	init_genrand( (unsigned long)time( NULL ));
	//init_genrand(0);
	int Mem_included=0;

	
		//cout<<"\n"<<"n="<<n<<"Multicast ratio="<<mul_ratio;
	int unicast_t=0;
	unicast_bp=0;	
	for (t=1;t<SimTime;t++)
	{	
		//cout<<"\n"<<"t="<<t<<" ";	
		LinkUpdate(arri_time);
		//Req.HoldTime=0.5;
		Req.HoldTime=exprnd(1);	
		Req.Type=0;
		Req.Snode=unif_int(Nv)-1;
		//cout<<"Snode="<<Req.Snode<<" ";			
		Req.MemNumber=1;	
		//cout<<"MemNumber="<<Req.MemNumber<<" ";
			
		for (i=0;i<Nv;i++)
		{	Req.Dnodes[i]=0;
				Req.Member[i]=0;
		}
		Req.Dnodes[Req.Snode]=1;
		//cout<<"Dnodes are: ";
		for (i=1;i<=Req.MemNumber;i=i+1)
		{	int temp=0;
			temp=unif_int(Nv-i);
			for(j=0,k=0;j<temp&&k<Nv;)
			{	if (Req.Dnodes[k]==0)	j=j+1;
					k=k+1;
			}
			Req.Dnodes[k-1]=1;
		}
			
		Req.Dnodes[Req.Snode]=0;
		Dijkstra();	
			//PrintPath();	
		WAA_down2up(); 
			//WAA_down2up_MWCn();
		Mem_included=0;
		for (i=0;i<Nv;i++)
		{
			if (Req.Member[i])
				Mem_included++;				
		}	

		unicast_bp=unicast_bp+(1-(float)Mem_included/(float)Req.MemNumber);
		unicast_t++;
				
		arri_time=exprnd(load);
	}

		if (unicast_t>0)
			unicast_bp=unicast_bp/float(unicast_t);
		else 
			cerr<<"number of unicast is 0 ";
			
    return 1;
}

int Graph::MulticastRequest(int SimTime,float load,float &multicast_bp)
{
    float arri_time=0;
	int i,j,k,t;
//	float load=0.5,mem_ratio=0.5;
	init_genrand( (unsigned long)time( NULL ));
	//init_genrand(0);
	int Mem_included=0;	
	//cout<<"\n"<<"n="<<n<<"Multicast ratio="<<mul_ratio;
	int multicast_t=0;
	multicast_bp=0;
	for (t=1;t<SimTime;t++)
	{	
			//cout<<"\n"<<"t="<<t<<" ";	
		LinkUpdate(arri_time);
			//Req.HoldTime=0.5;
		Req.HoldTime=exprnd(1);	
		Req.Type=1;
		Req.Snode=unif_int(Nv)-1;
			//cout<<"Snode="<<Req.Snode<<" ";
			
		Req.MemNumber=unif_int(Nv-1)+1;//[2,Nv]
			//cout<<"MemNumber="<<Req.MemNumber<<" ";
			
		for (i=0;i<Nv;i++)
		{	Req.Dnodes[i]=0;
			Req.Member[i]=0;
		}
		Req.Dnodes[Req.Snode]=1;
			//cout<<"Dnodes are: ";
		for (i=1;i<=Req.MemNumber;i=i+1)
		{	int temp=0;
			temp=unif_int(Nv-i);
			for(j=0,k=0;j<temp&&k<Nv;)
			{	if (Req.Dnodes[k]==0)	j=j+1;
					k=k+1;
			}
			Req.Dnodes[k-1]=1;
		}
			
		Req.Dnodes[Req.Snode]=0;
		Dijkstra();	
			//PrintPath();	
		WAA_down2up(); 
			//WAA_down2up_MWCn();
		Mem_included=0;
		for (i=0;i<Nv;i++)
		{
			if (Req.Member[i])
				Mem_included++;				
		}	
			
		multicast_bp=multicast_bp+(1-(float) Mem_included/(float)(Req.MemNumber));
		multicast_t++;
		arri_time=exprnd(load);
	}
		if (multicast_t>0)
			multicast_bp=multicast_bp/float(multicast_t);
		else 
			cerr<<"number of multicast is 0 ";
			
			
    return 1;
}

int Graph::ProcessRequest()
{
    float arri_time=0;
	int i;
///	float load=0.5,mem_ratio=0.5;
//	init_genrand( (unsigned long)time( NULL ));				
	float block_rate=0;
		LinkUpdate(arri_time);
		Req.HoldTime=1;
//		Req.HoldTime=exprnd(1);		
		Req.Snode=14-1;
		for (i=0;i<Nv;i++)	Req.Dnodes[i]=0;			
		Req.Dnodes[1-1]=1; 
		Req.Dnodes[2-1]=1; 
		Req.Dnodes[3-1]=1;
		Req.Dnodes[5-1]=1;
		Req.Dnodes[6-1]=1;
		Req.Dnodes[7-1]=1;
		//Req.Dnodes[8-1]=1;
		Req.Dnodes[9-1]=1; 
		Req.Dnodes[10-1]=1;
		Req.Dnodes[11-1]=1;
		Req.Dnodes[12-1]=1;
		Req.Dnodes[13-1]=1;
		//Req.Dnodes[14-1]=1;
		std::cout<<"\n"<<"\n"<<"Req: HoldTime is: "<<Req.HoldTime
			<<" ArriveTime is: "<<arri_time;
		std::cout<<"\n"<<"Source is: "<<Req.Snode+1<<"; ";
		std::cout<<"Dest is: ";
		Req.MemNumber=0;
		for (i=0;i<Nv;i++)
		{	Req.Member[i]=0;		
			if (Req.Dnodes[i])
			{	Req.MemNumber++;
				std::cout<<" "<<i+1<<" ";
			}
		}
		std::cout<<"\n";		
		Dijkstra();	
		PrintPath();	
		sim_type WAA_flag=sWC_up2down_v1;
		switch (WAA_flag)
		{	case MWC_down2up:
				block_rate+=WAA_down2up_MWCn();break;
			case fWC_down2up:
				block_rate+=WAA_down2up();break;
			case fWC_up2down:
				block_rate+=WAA_up2down(WAA_flag);break;
			case sWC_up2down:
				block_rate+=WAA_up2down(WAA_flag);break;
			case sWC_up2down_v1:
				ClearTable(2);
				block_rate+=WAA_up2down(WAA_flag);break;
			case MWC_sparse:
				block_rate+=WAA_down2up_Sparse_MWCn();break;
			default:break;
		}
		//arri_time=exprnd(load);
		arri_time=1;	
    return 1;
}




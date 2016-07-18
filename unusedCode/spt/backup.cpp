void Graph::ClearTable(int type )
{
			
    for( int i = 0; i < N_node; i++ )
    {
        Table[i].Dist = Inf;
        Table[i].Prev = -1;
        Table[i].Mark = 0;
		Table[i].Ontree=0;
		
		if (Gmap[i].MC==1)
			Table[i].degree=Gmap[i].degree;
		else
			Table[i].degree=2;
    }
//type==1, first time to initate Table
	if type==0
	{	for( int i = 0; i < N_node; i++ )
			Table[i].Onforest=0;
	}	
}


int Graph::Member_only( int Snode,int *Dnodes,int Mem_number)//muticast routing algorithm
{
	int V,U, W;
    BinaryHeap<link> PQ(link(Snode, 0));
    link Pnext;    // Stores the result of a DeleteMin
	list<Edge> ::iterator P;
    ClearTable(0);
    Table[Snode].Dist = 0;
    //PQ.Insert(link(Snode, 0));
	int Tree_number=0
	int Mem_included=0;
    while (Mem_included<Mem_number)
    {
        do
        {
            if( PQ.IsEmpty( ) )		
			{		
				Tree_number=Tree_number+1;
				ClearTable(1);
				Table[Snode].Dist = 0;
				PQ.Insert(link(Snode, 0));
				// return 1;
			}
               
            PQ.DeleteMin(Pnext);
        } while(Table[Pnext.Dest].Mark);

        V = Pnext.Dest;
        Table[V].Mark = 1;     // Mark vertex as being seen
		if (Dnodes[V]==0)	// V is not a Member
			for (P = Gmap[V].Adj.begin(); P != Gmap[V].Adj.end(); ++P)
			{
				W = (*P).Dest;
				float Cvw =(*P).Cost;
				if( Cvw < 0 )
				{	std::cerr << "Graph has negative edges" << endl; 
					return 0;
				}
				// update the dist of v's neighbour v
				if( Table[W].Dist > Table[V].Dist + Cvw )
				{
					Table[W].Dist = Table[V].Dist + Cvw;
					Table[W].Prev = V;
					PQ.Insert( link( W, Table[W].Dist ) );
				}
			}
		else //V is not a Member
		{	Mem_included=Mem_included+1;
			U=V;
			PQ.MakeEmpty( );
			do
			{	Table[U].Ontree=1;
				Table[U].Onforest=1;
				Table[U].degree=Table[U].degree-1;//u->v then u--;
				U=Table[U].Prev;
				Table[U].degree=Table[U].degree-1;//U<-u; then v--
			}while(U!=Snode);

			for( int i = 0; i < N_node; i++ )
			{
				if (Table[i].Ontree==0)
				{	Table[i].Dist = Inf;
					Table[i].Prev = -1;
					Table[i].Mark = 0;
				}
				else 
					if (Table[i].degree>0)
						for (P = Gmap[i].Adj.begin(); P != Gmap[i].Adj.end(); ++P)
						{
							W = (*P).Dest;
							float Cvw =(*P).Cost;
							if (Table[W].Ontree==0)
							{	if( Cvw < 0 )
							{	std::cerr << "Graph has negative edges" << endl;
									return 0;
								}
								if( Table[W].Ontree)
									PQ.Insert(link(W, Cvw));									
							}
						}				

			}

		}
    }
	return 1;
}
int Graph::Dijkstra( int SNode,int *Dnodes,int mem_number)
{

	int V, W;
    BinaryHeap<link> PQ(link(SNode, 0));
    link Vnext;    // Stores the result of a DeleteMin
	list<Edge> ::iterator P;
    ClearTable( );
    Table[SNode].Dist = 0;
    //PQ.Insert(link(SNode, 0));
	
	int mem_included=0;
    while (mem_included<mem_number)
    {
        do
        {
            if( PQ.IsEmpty( ) )
                return 1;
            PQ.DeleteMin(Vnext);
        } while(Table[Vnext.Dest].Mark);

        Vnext = Vnext.Dest;
        Table[Vnext].Mark = 1;     // Mark vertex as being seen
		for (P = Gmap[Vnext].Adj.begin(); P != Gmap[Vnext].Adj.end(); ++P)
        {
            W = *P.Dest;
            DistType Cvw = *P.Cost;
            if( Cvw < 0 )
            {
                cerr << "Graph has negative edges" << endl;
                return 0;
            }
            if( Table[W].Dist > Table[V].Dist + Cvw )
            {
                Table[W].Dist = Table[V].Dist + Cvw;
                Table[W].Prev = V;
                PQ.Insert( link<DistType>( W, Table[W].Dist ) );
            }
        }

		if Dnodes[W] = True;
			mem_included=mem_included+1;
    }
	return 1;
}

int Graph::MRA( int Snode,int *Dnodes,int mem_number)//muticast routing algorithm
{
	int V, W;
    BinaryHeap<link> PQ(link(SNode, 0));
    link Pnext;    // Stores the result of a DeleteMin
	list<Edge> ::iterator P;
    ClearTable( );
    Table[Snode].Dist = 0;
    //PQ.Insert(link(SNode, 0));
	
	int mem_included=0;
    while (mem_included<mem_number)
    {
        do
        {
            if( PQ.IsEmpty( ) )
                return 1;
            PQ.DeleteMin(Pnext);
        } while(Table[Pnext.Dest].Mark);

        Vnext = Pnext.Dest;
        Table[Vnext].Mark = 1;     // Mark vertex as being seen
		if Dnodes[Vnext]==0	// V is not a member
			for (P = Gmap[Vnext].Adj.begin(); P != Gmap[Vnext].Adj.end(); ++P)
			{
				W = *P.Dest;
				DistType Cvw = *P.Cost;
				if( Cvw < 0 )
				{	cerr << "Graph has negative edges" << endl; 
					return 0;
				}
				// update the dist of v's neighbour v
				if( Table[W].Dist > Table[Vnext].Dist + Cvw )ertex
				{
					Table[W].Dist = Table[Vnext].Dist + Cvw;
					Table[W].Prev = Vnext;
					PQ.Insert( link( W, Table[W].Dist ) );
				}
			}
		else //V is not a member
		{	mem_included=mem_included+1;
			U=Vnext;			
			do
			{	Table[U].Ontree=true;
				Table[U].deree=Table[U].degree-1;//u->v then u--;
				U=Table[U].prev;
				Table[U].deree=Table[U].degree-1;//U<-u; then v--
			}while(U!=Snode);

			for( int i = 0; i < N_node; i++ )
			{
				if Table[i].Ontree==false
				{	Table[i].Dist = Inf;
					Table[i].Prev = -1;
					Table[i].Mark = 0;
				}
				else 
					if (Table[i].deree>0)
						for (P = Gmap[i].Adj.begin(); P != Gmap[i].Adj.end(); ++P)
						{
							W = *P.Dest;
							DistType Cvw =*P.Cost;
							if Table[W].Ontree==0
							{	if( Cvw < 0 )
								{	cerr << "Graph has negative edges" << endl;
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

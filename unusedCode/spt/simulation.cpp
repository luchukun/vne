#include "routing.h"
#include "wavepath.h"
void output2txt(float* x,int size,char *filename)
{
	// open input file
	ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
        // copy file contents to cout

	int i=0;
    while (i<size) {
         	outFile<<x[i]<<" ";
		 i++;
    }
	outFile<<endl;
	outFile.close();
}

void output2txt(float x[11][15],int Hsize,int Vsize,char *filename)
{
	// open input file
	ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
        // copy file contents to cout

	int i=0,j=0;
	while(i<Vsize)
	{ j=0;
		while (j<Hsize) {
         		outFile<<x[i][j]<<" ";
			 j++;
		}
		outFile<<"\n";
		i++;
	}
	outFile<<endl;
	outFile.close();
}

int blocking_probability_sep()
{
	int SimTime=10;
	float load[9]={10,20,30,40,50,60,70,80,90};
	float multicast_ratio[9]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	float unicast_bp[9],multicast_bp[9],unicast_bp1[9],multicast_bp1[9];
	int size=9;
	Graph G; //参数16改成8是对的	
	G.Gen_gragh();
	//G.LinkRandom(link_ratio);
	G.fMC();
	output2txt(multicast_ratio,size,"mr.txt");
	for (int i=0;i<size;i++)
	{
		G.LinkClear();
		G.ProcessRequest(SimTime,100,multicast_ratio[i],unicast_bp[i],multicast_bp[i]);
		G.LinkClear();
		G.MulticastRequest(SimTime,load[i],multicast_bp1[i]);
		G.LinkClear();
		G.UnicastRequest(SimTime,100-load[i],unicast_bp1[i]);
	}
	output2txt(unicast_bp,size,"ubp.txt");
	output2txt(multicast_bp,size,"mbp.txt");
	output2txt(unicast_bp1,size,"ubp1.txt");
	output2txt(multicast_bp1,size,"mbp1.txt");
	cout<<"end of main";
	char c;
	cin>>c;
	return 0;
}

void block_probality(Graph &G)
{	
	int SimTime=10000;
	float load=256;
	float multicast_ratio[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
	//float multicast_ratio[11]={1,1,1,1,1,1,1,1,1,1};
	sim_type WAA_flag;//fWC_down2up,fWC_up2down,MWC_down2up,sWC_up2down,MWC_sparse//
	//G.LinkRandom(link_ratio);
	int Hsize=15,Vsize=11,size=11;
	//cout<<"\n"<<"Load is: "<<load;
/*
	float bp_MC_sparse[11][15];
	WAA_flag=MWC_sparse;
	int NodePriority[14];	

	G.MCMax(NodePriority);	
	for (int mr=0;mr<11;mr++)
	{	
		//std::cout<<"\n"<<"multicast_ratio ratio = "<<multicast_ratio_ratio[mr];
		//std::cout<<"\n"<<"Block probability= ";
		for (int j=0;j<=14;j++)
		{	
			for (int i=0;i<14;i++)
			{	
				if (i<j)
				G.Table[NodePriority[i]].MC=true;
				else G.Table[NodePriority[i]].MC=false;
				//std::cout<<"\n"<<"Node["<<NodePriority[i]+1<<"] is added to be MC.";
			}
			
			bp_MC_sparse[mr][j]=G.ProcessRequest(SimTime,load,multicast_ratio[mr],WAA_flag);
		}
	}
	output2txt(bp_MC_sparse,Hsize,Vsize,"bp_MC_sparse.txt");
*/	
	//fWC_down2up,fWC_up2down,MWC_down2up,sWC_up2down,MWC_sparse//
/*
	G.fMC();
	float bp_MWC_down2up[11];
	WAA_flag=MWC_down2up;
	for (int mr=0;mr<size;mr++)	
		bp_MWC_down2up[mr]=G.ProcessRequest(SimTime,load,multicast_ratio[mr],WAA_flag);		
	output2txt(bp_MWC_down2up,size,"bp_MWC_down2up.txt");
*/
	float bp_fWC_down2up[11];
	WAA_flag=fWC_down2up;
	for (int mr=0;mr<size;mr++)		
		bp_fWC_down2up[mr]=G.ProcessRequest(SimTime,load,multicast_ratio[mr],WAA_flag);		
	output2txt(bp_fWC_down2up,size,"bp_fWC_down2up.txt");
cout<<"\n"<<"end of fWC_down2up";
	float bp_fWC_up2down[11];
	WAA_flag=fWC_up2down;
	for (int mr=0;mr<size;mr++)	
		bp_fWC_up2down[mr]=G.ProcessRequest(SimTime,load,multicast_ratio[mr],WAA_flag);		
	output2txt(bp_fWC_up2down,size,"bp_fWC_up2down.txt");
cout<<"\n"<<"end of fWC_up2down";
	float bp_sWC_up2down[11];
	WAA_flag=sWC_up2down;
	for (int mr=0;mr<size;mr++)		
		bp_sWC_up2down[mr]=G.ProcessRequest(SimTime,load,multicast_ratio[mr],WAA_flag);		
	output2txt(bp_sWC_up2down,size,"bp_sWC_up2down.txt");
cout<<"\n"<<"end of sWC_up2down";
	float bp_sWC_up2down_v1[11];
	WAA_flag=sWC_up2down_v1;
	for (int mr=0;mr<size;mr++)		
		bp_sWC_up2down_v1[mr]=G.ProcessRequest(SimTime,load,multicast_ratio[mr],WAA_flag);	
		//G.ProcessRequest();
	output2txt(bp_sWC_up2down_v1,size,"bp_sWC_up2down_v1.txt");
cout<<"\n"<<"end of sWC_up2down_v1"	;
}


int main()
{	Graph G; //参数16改成8是对的
	G.Gen_gragh();
	block_probality(G);
	cout<<"\n"<<"end of Main"<<"\n"<<"Type any key to continue..";
	char c;
	cin.get(c);
	
}
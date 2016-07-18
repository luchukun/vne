#include "Graphs.h"
#include "Trees.h"
#include <ctime>
using namespace std;
void testKSP(Graph& G)
{
	int s,d;
	{
		for (int t=1;t<=8;t++)
		{
			cout<<"\n"<<"Input: source dest"<<"\n";
			if(cin>>s&&cin>>d&&s!=d)
			{	s--;d--;
				G.YenKSP(s,d);
				for (int k=0;k<G.n_path[s][d];k++){
					for (unsigned short i=1;i<=G.Kpaths[s][d][k][0]+1;i++)
						cout<<G.Kpaths[s][d][k][i]+1<<"-";
					cout<<"\n";
				}
			}
		}
	}

}


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
void output2txt(int size,float t,float load[],float max_utilization[],float success_rate[],float RC[],char *filename)
{
	// open input file
	ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
	outFile<<"\n"<<"N=( "<<minN<<" , "<<maxN<<" )";
	outFile<<"\n"<<"B=( "<<minB<<" , "<<maxB<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
#ifndef _Tree
	outFile<<"\n"<<"process time is "<<t<<"ms"; 
#else
	outFile<<"\n"<<"process time is "<<t<<"us"; 
#endif
	int i=0;
	outFile<<"\n"<<"load = ";
    while (i<size) {
         	outFile<<load[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         	outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         	outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"R/C = ";
    while (i<size) {
         	outFile<<RC[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}
void output2txt(int size,float t[],int numOfVM[],float max_utilization[],float success_rate[],float bw[],char *filename)
{
	// open input file
	ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
	outFile<<"\n"<<"N=( "<<minN<<" , "<<maxN<<" )";
	outFile<<"\n"<<"B=( "<<minB<<" , "<<maxB<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
	
	int i=0;
	
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t[i]<<"ms ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"VM= ";
    while (i<size) {
         outFile<<numOfVM[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}

void dynamic_embedding(char algorithm)
{
	Tree G;
	G.drawTree();
	//input
	int numOfreq=1e3;
	const int n=4,m=2;
	float load[n]={0.8,0.6,0.4,0.2};
	int n_test=4;
	//output
	float max_utilization[n]={0},success_rate[n]={0},RC[n]={0};
	// start simulation	
	time_t Atime,Btime; 
	time(&Atime); 
	float running_time=0;
	
	for (int i=0;i<n_test;i++)
		G.ProcessRequest(algorithm,numOfreq,load[i],max_utilization[i],success_rate[i],RC[i]);
	
	time(&Btime);
	running_time=1000*1000*(float)(Btime - Atime)/(float)(numOfreq*n_test);
	cout<<"process time is "<<running_time<<"us"<<endl;
	
	switch(algorithm)
	{
		case 'P':
			output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"tree_pertubation.txt");
			break;
		case 'B':
			output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"tree_backtracking.txt");
			break;
		case 'F':			 
			output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"tree_firstfit.txt");
			break;
		default:
			break;
	}	
}

void dynamic_embedding(char algorithm,bool enLProuting)
{
	Graph G;
#ifdef _FatTree
	G.drawFattree();
#endif
#ifdef _L2VL2
	G.drawL2VL2();
#endif
#ifdef _VL2
	G.drawVL2();
#endif

	//testKSP(G);
	//input
	int numOfreq=100;
	const int n=4,m=2;
	float load[n]={0.8,0.6,0.4,0.2};
	int n_test=4;
	//output
	float max_utilization[n]={0},success_rate[n]={0},RC[n]={0};
	// start simulation
	
	time_t Atime,Btime; 
	time(&Atime); 
	float running_time=0;
	
	for (int i=0;i<n_test;i++)
		G.ProcessRequest(algorithm,enLProuting,numOfreq,load[i],max_utilization[i],success_rate[i],RC[i]);
	
	time(&Btime);
	running_time=1000*(float)(Btime - Atime)/(float)(numOfreq*n_test);
	cout<<"process time is "<<running_time<<"ms"<<endl;
	
	switch(algorithm)
	{
		case 'P':
			if (enLProuting)
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"pertubation_lp.txt");
			else 
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"pertubation.txt");
			break;
		case 'L':
			if (enLProuting)
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"localsearch_lp.txt");
			else 
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"localsearch.txt");
			break;
		case 'E':
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"exhaustive.txt");
			break;
		case 'B':
			if (enLProuting)
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"backtracking_lp.txt");
			else
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"backtracking.txt");
			break;
		case 'F':
			if (enLProuting)
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"firstfit_lp.txt");
			else 
				output2txt(n_test,running_time,load,max_utilization,success_rate,RC,"firstfit.txt");
			break;
		default:
			break;
	}	
}

void static_embedding(char algorithm,bool enLProuting)
{
	Graph G;
#ifdef _FatTree
	G.drawFattree();
#endif
#ifdef _L2VL2
	G.drawL2VL2();
#endif
#ifdef _VL2
	G.drawVL2();
#endif

	//testKSP(G);
	//input
	int numOfreq=100;
	const int n=4,m=2;
	float p=1;
	int numOfVM[n]={2,4,6,8};
	int n_test=4;
	//output
	float max_utilization[n]={0},success_rate[n]={0},bandwidth_cost[n]={0};
	// start simulation
	
	time_t Atime,Btime; 
	
	float running_time[n]={0};
	
	for (int i=0;i<n_test;i++)
	{
		time(&Atime); 
		G.SingleRequest(algorithm,enLProuting,numOfreq,p,numOfVM[i],
						max_utilization[i],success_rate[i],bandwidth_cost[i]);	
		time(&Btime);
		running_time[i]=1000*(float)(Btime - Atime)/(float)(numOfreq*n_test);
		cout<<"process time is "<<running_time[i]<<"ms"<<endl;
	}

	
	switch(algorithm)
	{
		case 'P':
			if (enLProuting)
				output2txt(n_test,running_time,numOfVM,max_utilization,success_rate,bandwidth_cost,"s_pertubation_lp.txt");
			else 
				output2txt(n_test,running_time,numOfVM,max_utilization,success_rate,bandwidth_cost,"s_pertubation.txt");
			break;
		case 'B':
			if (enLProuting)
				output2txt(n_test,running_time,numOfVM,max_utilization,success_rate,bandwidth_cost,"s_backtracking_lp.txt");
			else 
				output2txt(n_test,running_time,numOfVM,max_utilization,success_rate,bandwidth_cost,"s_backtracking.txt");
			break;
		case 'E':
				output2txt(n_test,running_time,numOfVM,max_utilization,success_rate,bandwidth_cost,"s_exhaustive.txt");
			break;
		case 'F':
			if (enLProuting)
				output2txt(n_test,running_time,numOfVM,max_utilization,success_rate,bandwidth_cost,"s_firstfit_lp.txt");
			else 
				output2txt(n_test,running_time,numOfVM,max_utilization,success_rate,bandwidth_cost,"s_firstfit.txt");
			break;
		default:
			break;
	}	
}


int main()
{		
	//P: Pertubation; B:Backtracking; F:FirstFit
#ifndef _Tree	
	char algorithm='P';	bool enLProuting=0;
	//single request simulation	
	if(1){
		algorithm='F';	enLProuting=0;
		static_embedding(algorithm,enLProuting);
			algorithm='P';	enLProuting=0;
			static_embedding(algorithm,enLProuting);
				algorithm='B';	enLProuting=1;
				static_embedding(algorithm,enLProuting);
	}				
	// dynamic simulation
	else{
	algorithm='P';	enLProuting=0;
	dynamic_embedding(algorithm,enLProuting);
		algorithm='F';	
		dynamic_embedding(algorithm,enLProuting);
			algorithm='B';	enLProuting=1;
			dynamic_embedding(algorithm,enLProuting);
				}
#else
	char algorithm='P';	
	dynamic_embedding(algorithm);
		algorithm='F';
		dynamic_embedding(algorithm);
			algorithm='B';
			dynamic_embedding(algorithm);
#endif
	// finish
	char c;
	cin.get(c);
	cout<<"\n"<<"end of Main"<<"\n"<<"Type any key to continue..";
	
}
